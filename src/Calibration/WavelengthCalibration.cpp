#include <algorithm>
#include <chrono>
#include <SpectralEvaluation/Calibration/WavelengthCalibration.h>
#include <SpectralEvaluation/Calibration/ReferenceSpectrumConvolution.h>
#include <SpectralEvaluation/Calibration/WavelengthCalibrationByRansac.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Evaluation/WavelengthFit.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Spectra/SpectrumUtils.h>
#include <SpectralEvaluation/VectorUtils.h>
#include <SpectralEvaluation/File/TXTFile.h>
#include <SpectralEvaluation/File/File.h>

namespace novac
{

// ---------------------- Free functions used in calibration -------------------------

std::vector<double> GetPixelToWavelengthMappingFromFile(const std::string& clbFile)
{
    CSpectrum initialWavelengthCalibrationSpectrum;
    SpectrumIO::CTXTFile::ReadSpectrum(initialWavelengthCalibrationSpectrum, clbFile);
    std::vector<double> pixelToWavelengthMapping{ initialWavelengthCalibrationSpectrum.m_data, initialWavelengthCalibrationSpectrum.m_data + initialWavelengthCalibrationSpectrum.m_length };
    return pixelToWavelengthMapping;
}

// --------------------------- FraunhoferSpectrumGeneration ---------------------------

std::unique_ptr<CSpectrum> FraunhoferSpectrumGeneration::GetFraunhoferSpectrum(
    const std::vector<double>& pixelToWavelengthMapping,
    const ::Evaluation::CCrossSectionData& measuredInstrumentLineShape,
    const std::string& highResolutionOzoneCrossSection,
    double ozoneTotalColumn)
{
    std::pair<std::string, double> ozone{ highResolutionOzoneCrossSection, ozoneTotalColumn };
    std::vector< std::pair<std::string, double>> crossSections{ ozone };
    return GetFraunhoferSpectrum(pixelToWavelengthMapping, measuredInstrumentLineShape, crossSections);
}

std::unique_ptr<CSpectrum> FraunhoferSpectrumGeneration::GetFraunhoferSpectrum(
    const std::vector<double>& pixelToWavelengthMapping,
    const ::Evaluation::CCrossSectionData& measuredInstrumentLineShape,
    const std::vector<std::pair<std::string, double>>& highResolutionCrossSections)
{
    // Get the high res solar spectrum
    Evaluation::CCrossSectionData solarCrossSection;
    FileIo::ReadCrossSectionFile(this->solarAtlasFile, solarCrossSection);

    for each (auto crossSection in highResolutionCrossSections)
    {
        const std::string& path = crossSection.first;
        const double totalColumn = crossSection.second;

        // Turn the molecular absorption into an absorbance spectrum and multiply with the high res solar
        if (std::abs(totalColumn) > std::numeric_limits<double>::epsilon())
        {
            // Get the high res ozone spectrum
            Evaluation::CCrossSectionData ozoneCrossSection;
            FileIo::ReadCrossSectionFile(path, ozoneCrossSection);

            Mult(ozoneCrossSection.m_crossSection, -totalColumn);
            Exp(ozoneCrossSection.m_crossSection);
            std::vector<double> resampledOzoneCrossSection;
            Evaluation::Resample(ozoneCrossSection, solarCrossSection.m_waveLength, resampledOzoneCrossSection);
            Mult(resampledOzoneCrossSection, solarCrossSection.m_crossSection);
        }
    }

    // Generate a theoretical solar spectrum by convolving the high-res solar atlas with the measured slf
    auto startTime = std::chrono::steady_clock::now();
    std::vector<double> theoreticalFraunhoferSpectrumData;
    ::Evaluation::ConvolveReference(pixelToWavelengthMapping, measuredInstrumentLineShape, solarCrossSection, theoreticalFraunhoferSpectrumData, WavelengthConversion::None, ConvolutionMethod::Fft);
    auto stopTime = std::chrono::steady_clock::now();
    std::cout << "Convolution of Fraunhofer Reference took " << std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime).count() << " ms" << std::endl;

    std::unique_ptr<CSpectrum> theoreticalFraunhoferSpectrum = std::make_unique<CSpectrum>(pixelToWavelengthMapping, theoreticalFraunhoferSpectrumData);
    Normalize(*theoreticalFraunhoferSpectrum); // normalizes the intensity

    return theoreticalFraunhoferSpectrum;
}


bool MercuryCalibration(
    const CSpectrum& measuredMercurySpectrum,
    int polynomialOrder,
    double minimumWavelength,
    double maximumWavelength,
    std::vector<SpectrumDataPoint>& /*foundPeaks*/,
    std::vector<double> pixelToWavelengthPolynomial)
{
    if (measuredMercurySpectrum.m_length < 50)
    {
        return false;
    }
    else if (polynomialOrder < 0 || polynomialOrder > 5)
    {
        return false;
    }
    else if (maximumWavelength < minimumWavelength)
    {
        return false;
    }

    // List of known mercury lines, in nm(air)
    const std::vector<double> knownMercuryLines = { 284.7675, 302.1498, 312.5668, 313.1548, 313.1839, 365.0153, 365.4836, 366.3279, 398.3931, 404.6563, 435.8328 };



    return false;
}



/// --------------------------- Wavelength Calibration --------------------------------

// Helper structure for keeping track of the envelope
struct SpectrumEnvelope
{
    std::vector<double> pixel;
    std::vector<double> intensity;

    /// <summary>
    /// Adds one more point to this envelope by adding it to the end of 'pixel' and 'intensity'
    /// </summary>
    void Add(double x, double y)
    {
        this->pixel.push_back(x);
        this->intensity.push_back(y);
    }

    size_t Size() const
    {
        return this->pixel.size();
    }
};


WavelengthCalibrationSetup::WavelengthCalibrationSetup(const WavelengthCalibrationSettings& calibrationSettings)
    : settings(calibrationSettings)
{
}

SpectrometerCalibrationResult WavelengthCalibrationSetup::DoWavelengthCalibration(CSpectrum& measuredSpectrum, const Evaluation::CCrossSectionData& measuredInstrumentLineShape)
{
    // TODO: Validation of the setup and incoming parameters!

    // Magic parameters...
    const double minimumPeakIntensityInMeasuredSpectrum = 0.06; // in the normalized units, was 1000
    const double minimumPeakIntensityInFraunhoferReference = 0.01; // in the normalized units, was 600
    novac::RansacWavelengthCalibrationSettings ransacSettings; // Magic method parameters. These needs to be optimized...
    novac::CorrespondenceSelectionSettings correspondenceSelectionSettings; // Magic selection parameters...

    // Setup
    novac::RansacWavelengthCalibrationSetup ransacCalibrationSetup{ ransacSettings };

    // Start by normalizing the intensity of the measured spectrum, such that we can compare it to the fraunhofer spectrum.
    Normalize(measuredSpectrum);

    // Get the envelope of the measured spectrum (used to correct the shape of the fraunhofer spectrum to the detector sensitivity + optics absorption of the spectrometer)
    SpectrumEnvelope measuredSpectrumEnvelope;
    novac::GetEnvelope(measuredSpectrum, measuredSpectrumEnvelope.pixel, measuredSpectrumEnvelope.intensity);

    // Find the keypoints of the measured spectrum
    std::vector<novac::SpectrumDataPoint> measuredKeypoints;
    novac::FindKeypointsInSpectrum(measuredSpectrum, minimumPeakIntensityInMeasuredSpectrum, measuredKeypoints);

    // Get the Fraunhofer spectrum
    novac::FraunhoferSpectrumGeneration fraunhoferSetup{ settings.highResSolarAtlas };
    const auto originalFraunhoferSpectrum = fraunhoferSetup.GetFraunhoferSpectrum(settings.initialPixelToWavelengthMapping, measuredInstrumentLineShape, settings.crossSections);
    std::unique_ptr<CSpectrum> fraunhoferSpectrum = std::make_unique<CSpectrum>(*originalFraunhoferSpectrum); // create a copy which we can modify

    SpectrometerCalibrationResult result;
    result.pixelToWavelengthMapping = settings.initialPixelToWavelengthMapping;

    // Run the calibration a number of times, each time adjusting the shape of the Fraunhofer spectrum to better match the measured intensity
    const int numberOfIterations = 3;
    for (int iterationIdx = 0; iterationIdx < numberOfIterations; ++iterationIdx)
    {
        // Get all the keypoints from the fraunhofer spectrum
        std::vector<novac::SpectrumDataPoint> fraunhoferKeypoints;
        novac::FindKeypointsInSpectrum(*fraunhoferSpectrum, minimumPeakIntensityInFraunhoferReference, fraunhoferKeypoints);

        // List all possible correspondences (with some filtering applied).
        const auto possibleCorrespondences = novac::ListPossibleCorrespondences(measuredKeypoints, measuredSpectrum, fraunhoferKeypoints, *fraunhoferSpectrum, ransacSettings, correspondenceSelectionSettings);

        // The actual wavelength calibration by ransac
        auto startTime = std::chrono::steady_clock::now();
        auto ransacResult = ransacCalibrationSetup.DoWavelengthCalibration(possibleCorrespondences);
        auto stopTime = std::chrono::steady_clock::now();

        // Save the result
        result.pixelToWavelengthMappingCoefficients = ransacResult.bestFittingModelCoefficients;
        result.pixelToWavelengthMapping = GetPixelToWavelengthMapping(ransacResult.bestFittingModelCoefficients, settings.initialPixelToWavelengthMapping.size());

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
            // Re-convolve the Fraunhofer spectrum to get it on the new pixel grid
            fraunhoferSpectrum = fraunhoferSetup.GetFraunhoferSpectrum(result.pixelToWavelengthMapping, measuredInstrumentLineShape, settings.crossSections);

            // Create a wavelength -> intensity spline using the new wavelength calibration and the envelope of the measured spectrum
            // TODO: Move to separate function
            {
                std::vector<double> measuredSpectrumWavelength(measuredSpectrumEnvelope.Size());
                std::vector<double> measuredSpectrumIntensity(measuredSpectrumEnvelope.Size());
                for (size_t ii = 0; ii < measuredSpectrumEnvelope.Size(); ++ii)
                {
                    measuredSpectrumWavelength[ii] = novac::PolynomialValueAt(ransacResult.bestFittingModelCoefficients, measuredSpectrumEnvelope.pixel[ii]);
                    measuredSpectrumIntensity[ii] = measuredSpectrumEnvelope.intensity[ii];
                }
                MathFit::CVector modelInput(measuredSpectrumWavelength.data(), (int)measuredSpectrumWavelength.size(), 1, false);
                MathFit::CVector modelOutput(measuredSpectrumIntensity.data(), (int)measuredSpectrumIntensity.size(), 1, false);

                // Create a spline from the slit-function.
                MathFit::CCubicSplineFunction apparentSensitivitySpline(modelInput, modelOutput);

                // Adjust the shape of the Fraunhofer spectrum to match the measured
                for (size_t pixelIdx = 0; pixelIdx < fraunhoferSpectrum->m_length; ++pixelIdx)
                {
                    const double apparentSensitivity = apparentSensitivitySpline.GetValue((MathFit::TFitData)fraunhoferSpectrum->m_wavelength[pixelIdx]);
                    fraunhoferSpectrum->m_data[pixelIdx] *= apparentSensitivity;
                }
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