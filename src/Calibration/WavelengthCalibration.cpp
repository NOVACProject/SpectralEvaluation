#include <algorithm>
#include <chrono>
#include <fstream>
#include <SpectralEvaluation/Calibration/WavelengthCalibration.h>
#include <SpectralEvaluation/Calibration/ReferenceSpectrumConvolution.h>
#include <SpectralEvaluation/Calibration/WavelengthCalibrationByRansac.h>
#include <SpectralEvaluation/Calibration/FraunhoferSpectrumGeneration.h>
#include <SpectralEvaluation/Calibration/InstrumentLineShapeEstimation.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Evaluation/WavelengthFit.h>
#include <SpectralEvaluation/Math/PolynomialFit.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Spectra/SpectrumUtils.h>
#include <SpectralEvaluation/VectorUtils.h>
#include <SpectralEvaluation/File/TXTFile.h>
#include <SpectralEvaluation/File/File.h>
#include <SpectralEvaluation/Interpolation.h>

namespace novac
{

// ---------------------- Free functions used in calibration -------------------------

std::vector<double> GetPixelToWavelengthMappingFromFile(const std::string& clbFile)
{
    CSpectrum initialWavelengthCalibrationSpectrum;
    CTXTFile::ReadSpectrum(initialWavelengthCalibrationSpectrum, clbFile);

    if (initialWavelengthCalibrationSpectrum.m_wavelength.size() > 0)
    {
        return initialWavelengthCalibrationSpectrum.m_wavelength;
    }
    else
    {
        return std::vector<double> { initialWavelengthCalibrationSpectrum.m_data, initialWavelengthCalibrationSpectrum.m_data + initialWavelengthCalibrationSpectrum.m_length };
    }
}

/// <summary>
/// Very special baseline removal where all points below the baseline are set to the baseline level.
/// </summary>
/// <param name="spectrum"></param>
void RemoveBaseline(CSpectrum& spectrum)
{
    std::vector<double> spectrumValues(spectrum.m_data, spectrum.m_data + spectrum.m_length);
    std::vector<double> lowestPoints;
    FindNLowest(spectrumValues, 20, lowestPoints);
    const double baseline = Average(lowestPoints);
    for (long ii = 0; ii < spectrum.m_length; ++ii)
    {
        if (spectrum.m_data[ii] < baseline)
        {
            spectrum.m_data[ii] = baseline;
        }
    }
}

double MaxIntensity(const std::vector<novac::SpectrumDataPoint>& measuredPeaks)
{
    if (measuredPeaks.size() == 0)
    {
        return 0.0;
    }

    double maxIntensity = measuredPeaks[0].intensity;
    for (size_t idx = 1; idx < measuredPeaks.size(); ++idx)
    {
        maxIntensity = std::max(maxIntensity, measuredPeaks[idx].intensity);
    }

    return maxIntensity;
}

bool IsSaturatedAround(const CSpectrum& measuredSpectrum, double fractionalPixel)
{
    const size_t pixel = static_cast<size_t>(std::round(fractionalPixel));

    if (pixel > 0 && (measuredSpectrum.m_data[pixel - 1] / measuredSpectrum.m_data[pixel - 1]) < 0.01)
    {
        return true;
    }
    if (pixel < static_cast<size_t>(measuredSpectrum.m_length - 1) && (measuredSpectrum.m_data[pixel] / measuredSpectrum.m_data[pixel + 1]) < 0.01)
    {
        return true;
    }

    return false;
}

/// <summary>
/// Returns true if any of the provided SpectrumDataPoint:s shows a flat-top
/// </summary>
bool AnyLineIsSaturated(const CSpectrum& measuredSpectrum, const std::vector<novac::SpectrumDataPoint>& measuredPeaks)
{
    if (measuredPeaks.size() == 0)
    {
        return false;
    }

    for (size_t idx = 0; idx < measuredPeaks.size(); ++idx)
    {
        if (IsSaturatedAround(measuredSpectrum, measuredPeaks[idx].pixel))
        {
            return true;
        }
    }

    return false;
}

// --------------------------------------- Performing a pixel-to-wavelength calibration from a measured mercury spectrum ---------------------------------------

bool MercuryCalibration(
    const CSpectrum& measuredMercurySpectrum,
    int polynomialOrder,
    double minimumWavelength,
    double maximumWavelength,
    SpectrometerCalibrationResult& result,
    MercurySpectrumCalibrationState* state)
{
    if (measuredMercurySpectrum.m_length < 50)
    {
        return false;
    }
    else if (polynomialOrder < 1 || polynomialOrder > 3)
    {
        return false;
    }
    else if (maximumWavelength < minimumWavelength)
    {
        return false;
    }

    // List of known, strong mercury lines, in nm(air)
    struct emissionLine
    {
        emissionLine(double w, double s) : wavelength(w), intensity(s) { }
        double wavelength;
        double intensity;
    };

    const std::vector<emissionLine> knownMercuryLines = {
        emissionLine(365.4836, 1.00),
        emissionLine(313.1839, 0.90),
        emissionLine(404.6563, 0.47),
        emissionLine(296.7280, 0.41),
        emissionLine(366.3279, 0.16),
        emissionLine(302.1498, 0.08),
        emissionLine(289.4,    0.04),
        emissionLine(334.148,  0.04),
        emissionLine(407.783,  0.03)
    };
    /* More known lines :
        253.652,  296.7280,
        302.1498, 312.5668, 313.1548, 313.1839, 334.148, 365.0153, 365.4836, 366.3279, 398.3931,
        404.6563, 407.783, 435.8328,
        546.074, 576.960, 579.066 }; */

        // Find the peaks of the spectrum
    std::vector<novac::SpectrumDataPoint> measuredPeaks;
    FindEmissionLines(measuredMercurySpectrum, measuredPeaks);

    if (measuredPeaks.size() < polynomialOrder + 1)
    {
        // TODO: Set the state here?
        return false;
    }

    // Figure out how many points are saturated (if any)
    std::vector<novac::SpectrumDataPoint> sortedPeaks{ begin(measuredPeaks), end(measuredPeaks) };
    std::sort(begin(sortedPeaks), end(sortedPeaks), [](const novac::SpectrumDataPoint& a, const novac::SpectrumDataPoint& b) {return a.intensity > b.intensity; });
    const double maximumMeasuredIntensity = sortedPeaks[0].intensity;

    int numberOfSaturatedPeaks = 0;
    for (size_t ii = 0; ii < sortedPeaks.size(); ++ii)
    {
        if (sortedPeaks[ii].flatTop && sortedPeaks[ii].intensity > 0.8 * maximumMeasuredIntensity)
        {
            ++numberOfSaturatedPeaks;
        }
    }

    // Construct a model by assuming that the N strongest peaks in the measured are the N strongest known mercury peaks.
    std::vector<double> initialCalibration;
    bool initialCalibrationSuccess = false;
    if (numberOfSaturatedPeaks == 0)
    {
        std::vector<double> pixels;
        std::vector<double> wavelengths;
        for (size_t ii = 0; ii < 5; ++ii)
        {
            pixels.push_back(sortedPeaks[ii].pixel);
            wavelengths.push_back(knownMercuryLines[ii].wavelength);
        }
        std::sort(begin(pixels), end(pixels));
        std::sort(begin(wavelengths), end(wavelengths));

        PolynomialFit polyFit{ static_cast<int>(polynomialOrder) };
        initialCalibrationSuccess = polyFit.FitPolynomial(pixels, wavelengths, initialCalibration);
    }

    // Try to figure out which known mercury line corresponds to which measured peak
    std::vector<Correspondence> correspondences;
    const double initialWavelengthTolerance = initialCalibrationSuccess ? 10.0 : 50.0;
    const double relativePositionTolerance = (numberOfSaturatedPeaks > 0) ? 0.5 : 0.25;
    for (size_t measIdx = 0; measIdx < sortedPeaks.size(); ++measIdx)
    {
        const double initialWavelengthOfPeak = initialCalibrationSuccess ? PolynomialValueAt(initialCalibration, sortedPeaks[measIdx].pixel) : minimumWavelength + sortedPeaks[measIdx].pixel * (maximumWavelength - minimumWavelength) / (double)(measuredMercurySpectrum.m_length - 1);
        const double positionInList = measIdx / (double)sortedPeaks.size();

        for (size_t lineIdx = 0; lineIdx < knownMercuryLines.size(); ++lineIdx)
        {
            const double linePositionInList = (lineIdx / (double)knownMercuryLines.size());

            if (std::abs(initialWavelengthOfPeak - knownMercuryLines[lineIdx].wavelength) < initialWavelengthTolerance &&
                std::abs(positionInList - linePositionInList) < relativePositionTolerance)
            {
                novac::Correspondence corr;
                corr.measuredIdx = measIdx;
                corr.measuredValue = sortedPeaks[measIdx].pixel;
                corr.theoreticalIdx = lineIdx;
                corr.theoreticalValue = knownMercuryLines[lineIdx].wavelength;
                corr.error = std::abs(positionInList - linePositionInList);

                correspondences.push_back(corr);
            }
        }
    }

    // Run the calibration using Ransac
    RansacWavelengthCalibrationSettings calibrationSettings;
    calibrationSettings.modelPolynomialOrder = polynomialOrder;
    calibrationSettings.numberOfRansacIterations = 50000;
#ifdef DEBUG
    calibrationSettings.numberOfThreads = 4; // auto
#else
    calibrationSettings.numberOfThreads = 1;
#endif // DEBUG

    RansacWavelengthCalibrationSetup calibrationSetup{ calibrationSettings };

    auto ransacResult = calibrationSetup.DoWavelengthCalibration(correspondences);

    if (ransacResult.highestNumberOfInliers == 0)
    {
        return false;
    }

    // Save the result
    result.pixelToWavelengthMappingCoefficients = ransacResult.bestFittingModelCoefficients;
    result.pixelToWavelengthMapping = std::vector<double>(measuredMercurySpectrum.m_length);
    for (size_t pixelIdx = 0; pixelIdx < static_cast<size_t>(measuredMercurySpectrum.m_length); ++pixelIdx)
    {
        result.pixelToWavelengthMapping[pixelIdx] = PolynomialValueAt(ransacResult.bestFittingModelCoefficients, static_cast<double>(pixelIdx));
    }

    if (nullptr != state)
    {
        // Remap the points from the (possible) wavelength calibration of the measured spectrum to
        //  our own (home brewn) pixel-to-wavelength calibration.   
        for (auto& peak : measuredPeaks)
        {
            novac::LinearInterpolation(result.pixelToWavelengthMapping, peak.pixel, peak.wavelength);
        }

        state->peaks = measuredPeaks;
    }

    return true;
}



/// --------------------------- Wavelength Calibration --------------------------------


WavelengthCalibrationSetup::WavelengthCalibrationSetup(const WavelengthCalibrationSettings& calibrationSettings)
    : settings(calibrationSettings)
{
}

SpectrometerCalibrationResult WavelengthCalibrationSetup::DoWavelengthCalibration(const CSpectrum& measuredSpectrum)
{
    // TODO: Validation of the setup and incoming parameters!

    // Magic parameters...
    const double minimumPeakIntensityInMeasuredSpectrum = 0.02; // in the normalized units, was 1000
    const double minimumPeakIntensityInFraunhoferReference = 0.01; // in the normalized units, was 600
    novac::RansacWavelengthCalibrationSettings ransacSettings; // Magic method parameters. These needs to be optimized...
    novac::CorrespondenceSelectionSettings correspondenceSelectionSettings; // Magic selection parameters...

    // Setup
    novac::RansacWavelengthCalibrationSetup ransacCalibrationSetup{ ransacSettings };

    // Start by removing any remaining baseline from the measuerd spectrum and normalizing the intensity of it, such that we can compare it to the fraunhofer spectrum.
    this->calibrationState.measuredSpectrum = std::make_unique<CSpectrum>(measuredSpectrum);
    RemoveBaseline(*calibrationState.measuredSpectrum);
    Normalize(*calibrationState.measuredSpectrum);

    // Get the envelope of the measured spectrum (used to correct the shape of the fraunhofer spectrum to the detector sensitivity + optics absorption of the spectrometer)
    novac::GetEnvelope(*calibrationState.measuredSpectrum, calibrationState.measuredSpectrumEnvelopePixels, calibrationState.measuredSpectrumEnvelopeIntensities);

    // Find the keypoints of the measured spectrum
    novac::FindKeypointsInSpectrum(*calibrationState.measuredSpectrum, minimumPeakIntensityInMeasuredSpectrum, calibrationState.measuredKeypoints);

    // Get the Fraunhofer spectrum
    novac::FraunhoferSpectrumGeneration fraunhoferSetup{ settings.highResSolarAtlas, settings.crossSections };
    calibrationState.originalFraunhoferSpectrum = fraunhoferSetup.GetFraunhoferSpectrum(settings.initialPixelToWavelengthMapping, settings.initialInstrumentLineShape);
    // calibrationState.originalFraunhoferSpectrum = fraunhoferSetup.GetFraunhoferSpectrumMatching(settings.initialPixelToWavelengthMapping, *calibrationState.measuredSpectrum, measuredInstrumentLineShape);
    calibrationState.fraunhoferSpectrum = std::make_unique<CSpectrum>(*calibrationState.originalFraunhoferSpectrum); // create a copy which we can modify

    SpectrometerCalibrationResult result;
    result.pixelToWavelengthMapping = settings.initialPixelToWavelengthMapping;

    // Run the calibration a number of times, each time adjusting the shape of the Fraunhofer spectrum to better match the measured intensity
    const int numberOfIterations = 2;
    for (int iterationIdx = 0; iterationIdx < numberOfIterations; ++iterationIdx)
    {
        // Get all the keypoints from the fraunhofer spectrum
        novac::FindKeypointsInSpectrum(*calibrationState.fraunhoferSpectrum, minimumPeakIntensityInFraunhoferReference, calibrationState.fraunhoferKeypoints);

        // List all possible correspondences (with some filtering applied).
        this->calibrationState.allCorrespondences = novac::ListPossibleCorrespondences(calibrationState.measuredKeypoints, *calibrationState.measuredSpectrum, calibrationState.fraunhoferKeypoints, *calibrationState.fraunhoferSpectrum, correspondenceSelectionSettings);

        // The actual wavelength calibration by ransac
        auto startTime = std::chrono::steady_clock::now();
        auto ransacResult = ransacCalibrationSetup.DoWavelengthCalibration(calibrationState.allCorrespondences);
        auto stopTime = std::chrono::steady_clock::now();

        // Save the result
        result.pixelToWavelengthMappingCoefficients = ransacResult.bestFittingModelCoefficients;
        result.pixelToWavelengthMapping = GetPixelToWavelengthMapping(ransacResult.bestFittingModelCoefficients, settings.initialPixelToWavelengthMapping.size());
        calibrationState.correspondenceIsInlier = ransacResult.correspondenceIsInlier;

        {
            // Output for debugging
            std::cout << "Calibration by Ransac took " << std::chrono::duration_cast<std::chrono::seconds>(stopTime - startTime).count() << " seconds" << std::endl;
            std::cout << "Best fitting model includes " << ransacResult.highestNumberOfInliers << " inliers out of the " << ransacResult.numberOfPossibleCorrelations << " possible correspondences" << std::endl;
            std::cout << "Best fitting model: " << std::endl;
            for (int orderIdx = 0; orderIdx <= (int)ransacResult.modelPolynomialOrder; ++orderIdx)
            {
                std::cout << "  c[" << orderIdx << "]: " << ransacResult.bestFittingModelCoefficients[orderIdx] << std::endl;
            }
        }

        std::cout << " -- Adjusting Fraunhofer Spectrum -- " << std::endl;

        // Adjust the shape of the Fraunhofer spectrum using what we now know:
        //  We have the sensitivity of the spectrometer (pixel -> intesity) from the calculated envelope
        //  We have the pixel-to-wavelength mapping of the spectrometer from the ransac calibration
        // This allows us to calculate the wavelength -> intensity mapping of the 
        if (iterationIdx < numberOfIterations - 1)
        {
            // Adjust the selection parameter for the maximum error in the wavelength calibration such that the search space decreases for each iteration
            correspondenceSelectionSettings.maximumPixelDistanceForPossibleCorrespondence = std::max(10, correspondenceSelectionSettings.maximumPixelDistanceForPossibleCorrespondence / 2);

            if (this->settings.estimateInstrumentLineShape == InstrumentLineshapeEstimationOption::Gaussian)
            {
                InstrumentLineShapeEstimation ilsEstimator{ result.pixelToWavelengthMapping };
                if (settings.initialInstrumentLineShape.GetSize() > 0)
                {
                    ilsEstimator.UpdateInitialLineShape(settings.initialInstrumentLineShape);
                }

                double lineShapeFwhm = 0.0;
                ilsEstimator.EstimateInstrumentLineShape(fraunhoferSetup, *calibrationState.measuredSpectrum, result.estimatedInstrumentLineShape, lineShapeFwhm);

                // TODO: Save the result !!

                // Re-convolve the Fraunhofer spectrum to get it on the new pixel grid and with the new instrument line shape
                calibrationState.fraunhoferSpectrum = fraunhoferSetup.GetFraunhoferSpectrum(result.pixelToWavelengthMapping, result.estimatedInstrumentLineShape);
            }
            else
            {
                // Re-convolve the Fraunhofer spectrum to get it on the new pixel grid
                calibrationState.fraunhoferSpectrum = fraunhoferSetup.GetFraunhoferSpectrum(result.pixelToWavelengthMapping, settings.initialInstrumentLineShape);
                // calibrationState.fraunhoferSpectrum = fraunhoferSetup.GetFraunhoferSpectrumMatching(result.pixelToWavelengthMapping, *calibrationState.measuredSpectrum, measuredInstrumentLineShape);
            }

            // Create a wavelength -> intensity spline using the new wavelength calibration and the envelope of the measured spectrum
            // TODO: Move to separate function
            {
                std::vector<double> measuredSpectrumWavelength(calibrationState.measuredSpectrumEnvelopePixels.size());
                for (size_t ii = 0; ii < calibrationState.measuredSpectrumEnvelopePixels.size(); ++ii)
                {
                    measuredSpectrumWavelength[ii] = novac::PolynomialValueAt(ransacResult.bestFittingModelCoefficients, calibrationState.measuredSpectrumEnvelopePixels[ii]);
                }
                MathFit::CVector modelInput(measuredSpectrumWavelength.data(), (int)measuredSpectrumWavelength.size(), 1, false);
                MathFit::CVector modelOutput(calibrationState.measuredSpectrumEnvelopeIntensities.data(), (int)calibrationState.measuredSpectrumEnvelopeIntensities.size(), 1, false);

                // Create a spline from the slit-function.
                MathFit::CCubicSplineFunction apparentSensitivitySpline(modelInput, modelOutput);

                // Adjust the shape of the Fraunhofer spectrum to match the measured
                for (size_t pixelIdx = 0; pixelIdx < calibrationState.fraunhoferSpectrum->m_length; ++pixelIdx)
                {
                    const double apparentSensitivity = apparentSensitivitySpline.GetValue((MathFit::TFitData)calibrationState.fraunhoferSpectrum->m_wavelength[pixelIdx]);
                    calibrationState.fraunhoferSpectrum->m_data[pixelIdx] *= apparentSensitivity;
                }
                Normalize(*calibrationState.fraunhoferSpectrum);
            }
        }
    }

    return result;
}

std::vector<double> WavelengthCalibrationSetup::GetPixelToWavelengthMapping(const std::vector<double>& polynomialCoefficients, size_t detectorSize)
{
    std::vector<double> result(detectorSize);

    for (size_t pixelIdx = 0; pixelIdx < detectorSize; ++pixelIdx)
    {
        result[pixelIdx] = novac::PolynomialValueAt(polynomialCoefficients, (double)pixelIdx);
    }

    return result;
}


}