#include <algorithm>
#include <chrono>
#include <fstream>
#include <SpectralEvaluation/Calibration/Correspondence.h>
#include <SpectralEvaluation/Calibration/WavelengthCalibration.h>
#include <SpectralEvaluation/Calibration/ReferenceSpectrumConvolution.h>
#include <SpectralEvaluation/Calibration/WavelengthCalibrationByRansac.h>
#include <SpectralEvaluation/Calibration/FraunhoferSpectrumGeneration.h>
#include <SpectralEvaluation/Calibration/CrossSectionSpectrumGenerator.h>
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

    std::vector<double> GetPixelToWavelengthMapping(const std::vector<double>& polynomialCoefficients, size_t detectorSize)
    {
        std::vector<double> result(detectorSize);

        for (size_t pixelIdx = 0; pixelIdx < detectorSize; ++pixelIdx)
        {
            result[pixelIdx] = novac::PolynomialValueAt(polynomialCoefficients, (double)pixelIdx);
        }

        return result;
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

    bool IsInlier(size_t measuredIdx, const std::vector<Correspondence>& allCorrespondences, const std::vector<bool>& correspondenceIsInlier)
    {
        for (size_t correspondenceIdx = 0; correspondenceIdx < allCorrespondences.size(); ++correspondenceIdx)
        {
            if (correspondenceIsInlier[correspondenceIdx] && allCorrespondences[correspondenceIdx].measuredIdx == measuredIdx)
            {
                return true;
            }
        }

        return false; // not found
    }

    size_t FindPeakWithPixel(double pixel, const std::vector<novac::SpectrumDataPoint>& measuredPeaks)
    {
        size_t closestPeak = 0;
        double closestDistance = std::numeric_limits<double>::max();

        for (size_t peakIdx = 0; peakIdx < measuredPeaks.size(); ++peakIdx)
        {
            const double distance = std::abs(pixel - measuredPeaks[peakIdx].pixel);
            if (distance < 0.001)
            {
                return peakIdx; // no doubt
            }
            if (distance < closestDistance)
            {
                closestDistance = distance;
                closestPeak = peakIdx;
            }
        }

        return closestPeak;
    }

    // Try to figure out which known mercury line corresponds to which measured peak, listing possible correspondences
    std::vector<Correspondence> ListMercurySpectrumCorrespondences(
        const std::vector<novac::SpectrumDataPoint>& measuredPeaks,
        const std::vector<double>& initialPixelToWavelength,
        double initialWavelengthTolerance,
        double relativePositionTolerance)
    {
        // Figure out how many points are saturated (if any)
        // sortedPeaks is the emission lines sorted in decreasing intensity
        std::vector<novac::SpectrumDataPoint> sortedPeaks{ begin(measuredPeaks), end(measuredPeaks) };
        std::sort(begin(sortedPeaks), end(sortedPeaks), [](const novac::SpectrumDataPoint& a, const novac::SpectrumDataPoint& b) {return a.intensity > b.intensity; });

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
            emissionLine(312.5668, 0.70),
            emissionLine(404.6563, 0.47),
            emissionLine(296.7280, 0.41),
            emissionLine(366.3279, 0.16),
            emissionLine(302.1498, 0.08),
            emissionLine(289.4,    0.04),
            emissionLine(334.148,  0.04),
            emissionLine(407.783,  0.03)
        };
        /* More known lines, could be added later if needed:
            253.652,  296.7280,
            302.1498, 312.5668, 313.1548, 313.1839, 334.148, 365.0153, 365.4836, 366.3279, 398.3931,
            404.6563, 407.783, 433.92, 434.75, 435.84
            546.074, 576.960, 579.066 };
        */

        std::vector<Correspondence> correspondences;
        for (size_t sortedPeaksIdx = 0; sortedPeaksIdx < sortedPeaks.size(); ++sortedPeaksIdx)
        {
            double initialWavelengthOfPeak = 0.0;
            if (!novac::LinearInterpolation(initialPixelToWavelength, sortedPeaks[sortedPeaksIdx].pixel, initialWavelengthOfPeak))
            {
                continue;
            }
            const double positionInList = sortedPeaksIdx / (double)sortedPeaks.size();
            int nofKnownLinesAssociated = 0;
            double wavelengthTolerance = initialWavelengthTolerance;

            while (nofKnownLinesAssociated == 0 && wavelengthTolerance < 100.0)
            {
                for (size_t lineIdx = 0; lineIdx < knownMercuryLines.size(); ++lineIdx)
                {
                    const double linePositionInList = (lineIdx / (double)knownMercuryLines.size());

                    if (std::abs(initialWavelengthOfPeak - knownMercuryLines[lineIdx].wavelength) < wavelengthTolerance &&
                        std::abs(positionInList - linePositionInList) < relativePositionTolerance)
                    {
                        novac::Correspondence corr;
                        corr.measuredIdx = FindPeakWithPixel(sortedPeaks[sortedPeaksIdx].pixel, measuredPeaks);
                        corr.measuredValue = sortedPeaks[sortedPeaksIdx].pixel;
                        corr.theoreticalIdx = lineIdx;
                        corr.theoreticalValue = knownMercuryLines[lineIdx].wavelength;
                        corr.error = std::abs(positionInList - linePositionInList);

                        correspondences.push_back(corr);

                        ++nofKnownLinesAssociated;
                    }
                }

                if (nofKnownLinesAssociated == 0)
                {
                    wavelengthTolerance += 10.0;
                }
            }
        }

        return correspondences;
    }

    bool ExistsSaturatedPeak(const std::vector<novac::SpectrumDataPoint>& peaks)
    {
        if (peaks.size() == 0)
        {
            return false;
        }
        double maximumIntensity = peaks.front().intensity;
        for (size_t ii = 1; ii < peaks.size(); ++ii)
        {
            maximumIntensity = std::max(maximumIntensity, peaks[ii].intensity);
        }

        for (size_t ii = 0; ii < peaks.size(); ++ii)
        {
            if (peaks[ii].flatTop && peaks[ii].intensity > 0.8 * maximumIntensity)
            {
                return true;
            }
        }

        return false;
    }

    bool MercuryCalibration(
        const CSpectrum& measuredMercurySpectrum,
        int polynomialOrder,
        const std::vector<double>& initialPixelToWavelength,
        SpectrometerCalibrationResult& result,
        MercurySpectrumCalibrationState* state)
    {
        if (measuredMercurySpectrum.m_length < 50)
        {
            if (state != nullptr)
            {
                state->errorMessage = "To short measured spectrum, at least 50 data poins is needed";
            }
            return false;
        }
        else if (polynomialOrder < 1 || polynomialOrder > 3)
        {
            if (state != nullptr)
            {
                state->errorMessage = "Polynomial order must be between 1 and 3 inclusive";
            }
            return false;
        }
        else if (initialPixelToWavelength.size() != static_cast<size_t>(measuredMercurySpectrum.m_length))
        {
            if (state != nullptr)
            {
                state->errorMessage = "Bad initial wavelength range provided, the length of the mapping must match the length of the measured spectrum";
            }
            return false;
        }

        // measuredPeaks is the emission lines of the spectrum, sorted in increasing pixel order.
        // TODO: Since we are now able to remove not fully resolved spectrum peaks (setting the last flag here to 'false')
        //  then we should also be able to improve the fitting below by only using fully-resolved peaks (and also getting info
        //  on which peaks have been removed (e.g. by calling FindEmissionLines twice, once with 'true' and once with 'false')
        std::vector<novac::SpectrumDataPoint> measuredPeaks;
        FindEmissionLines(measuredMercurySpectrum, measuredPeaks, true);

        auto unresolvedPeaks = novac::FilterByType(measuredPeaks, novac::SpectrumDataPointType::UnresolvedPeak);
        measuredPeaks = novac::FilterByType(measuredPeaks, novac::SpectrumDataPointType::Peak);

        if (measuredPeaks.size() < static_cast<size_t>(polynomialOrder) + 1)
        {
            if (state != nullptr)
            {
                state->errorMessage = "Not enough resolved emission lines could be found";
            }
            return false;
        }

        const double relativePositionTolerance = (unresolvedPeaks.size() > 0 || ExistsSaturatedPeak(measuredPeaks)) ? 0.5 : 0.25;

        // Run the calibration using Ransac
        RansacWavelengthCalibrationSettings calibrationSettings;
        calibrationSettings.modelPolynomialOrder = polynomialOrder;
        calibrationSettings.numberOfRansacIterations = 50000;
        calibrationSettings.inlierLimitInWavelength = 0.1;
        calibrationSettings.detectorSize = static_cast<size_t>(measuredMercurySpectrum.m_length);
        calibrationSettings.refine = false;

#ifdef DEBUG
        calibrationSettings.numberOfThreads = 4; // auto
#else
        calibrationSettings.numberOfThreads = 1;
#endif // DEBUG

        RansacWavelengthCalibrationSetup calibrationSetup{ calibrationSettings };

        double initialWavelengthTolerance = 10.0; //< allowed tolerance in the initial calibration
        novac::RansacWavelengthCalibrationResult ransacResult;
        std::vector<Correspondence> correspondences;

        while (true)
        {
            // Try to figure out which known mercury line corresponds to which measured peak, listing possible correspondences
            correspondences = ListMercurySpectrumCorrespondences(measuredPeaks, initialPixelToWavelength, initialWavelengthTolerance, relativePositionTolerance);

            // Do the calibration
            ransacResult = calibrationSetup.DoWavelengthCalibration(correspondences);

            // Check if the calibration succeeded.
            if (ransacResult.highestNumberOfInliers > 0)
            {
                // we're done
                break;
            }

            // Increase the tolerance and attempt again.
            initialWavelengthTolerance *= 3.0;

            if (initialWavelengthTolerance > 100.0)
            {
                // Time to give up..
                if (state != nullptr)
                {
                    state->errorMessage = "Wavelength calibration failed";
                }
                return false;
            }
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
            // Extract the measured peaks which we have identified in the measured spectrum (and the ones we haven't identified)
            // Remap the points from the (possible) wavelength calibration of the measured spectrum to
            //  our freshly produced pixel-to-wavelength calibration.
            state->peaks.clear();
            for (size_t ii = 0; ii < measuredPeaks.size(); ++ii)
            {
                novac::LinearInterpolation(result.pixelToWavelengthMapping, measuredPeaks[ii].pixel, measuredPeaks[ii].wavelength);

                state->peaks.push_back(measuredPeaks[ii]);
            }

            state->rejectedPeaks.clear();
            for (size_t ii = 0; ii < unresolvedPeaks.size(); ++ii)
            {
                novac::LinearInterpolation(result.pixelToWavelengthMapping, unresolvedPeaks[ii].pixel, unresolvedPeaks[ii].wavelength);

                state->rejectedPeaks.push_back(unresolvedPeaks[ii]);
            }
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
        if (settings.initialPixelToWavelengthMapping.size() != static_cast<size_t>(measuredSpectrum.m_length))
        {
            std::stringstream message;
            message << "The initial pixel to wavelength mapping must have as many points as the measured spectrum.";
            message << " Measured spectrum length is: " << measuredSpectrum.m_length;
            message << " initial pixel to wavelength mapping length is: " << settings.initialPixelToWavelengthMapping.size();
            throw std::invalid_argument(message.str());
        }
        else if (settings.initialInstrumentLineShape.GetSize() < 10)
        {
            throw std::invalid_argument("The initial instrument line shape must not be empty.");
        }

        // Magic parameters...
        const double minimumPeakIntensityInMeasuredSpectrum = 0.02; // in the normalized units.
        const double minimumPeakIntensityInFraunhoferReference = 0.01; // in the normalized units.
        novac::RansacWavelengthCalibrationSettings ransacSettings;
        novac::CorrespondenceSelectionSettings correspondenceSelectionSettings;

        // Setup
        novac::RansacWavelengthCalibrationSetup ransacCalibrationSetup{ ransacSettings };

        // Start by removing any remaining baseline from the measured spectrum and normalizing the intensity of it, such that we can compare it to the fraunhofer spectrum.
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
        calibrationState.fraunhoferSpectrum = std::make_unique<CSpectrum>(*calibrationState.originalFraunhoferSpectrum); // create a copy which we can modify

        // Get the ozone setup.
        std::unique_ptr<CrossSectionSpectrumGenerator> ozoneSetup;
        if (settings.crossSectionsForInstrumentLineShapeFitting.size() > 0)
        {
            ozoneSetup = std::make_unique< CrossSectionSpectrumGenerator>(settings.crossSectionsForInstrumentLineShapeFitting.front());
        }

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
            if (this->calibrationState.allCorrespondences.size() <= ransacSettings.modelPolynomialOrder)
            {
                throw WavelengthCalibrationFailureException("Wavelength calibration failed, no possible correspondences were found.");
            }

            // The actual wavelength calibration by ransac
            auto startTime = std::chrono::steady_clock::now();
            auto ransacResult = ransacCalibrationSetup.DoWavelengthCalibration(calibrationState.allCorrespondences);
            auto stopTime = std::chrono::steady_clock::now();

            // Save the result
            result.pixelToWavelengthMappingCoefficients = ransacResult.bestFittingModelCoefficients;
            result.pixelToWavelengthMapping = GetPixelToWavelengthMapping(ransacResult.bestFittingModelCoefficients, settings.initialPixelToWavelengthMapping.size());
            result.pixelToWavelengthMappingError = ransacResult.smallestError;
            result.pixelToWavelengthMappingInliers = ransacResult.highestNumberOfInliers;
            result.pixelToWavelengthMappingPixelRange = ransacResult.largestPixelSpan;
            calibrationState.correspondenceIsInlier = ransacResult.correspondenceIsInlier;

            {
                // Output for debugging
                std::cout << "Calibration by Ransac took " << std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime).count() << " ms" << std::endl;
                std::cout << "Best fitting model includes " << ransacResult.highestNumberOfInliers << " inliers out of the " << ransacResult.numberOfPossibleCorrelations << " possible correspondences" << std::endl;
                std::cout << "Best fitting model: " << std::endl;
                for (int orderIdx = 0; orderIdx <= (int)ransacResult.modelPolynomialOrder; ++orderIdx)
                {
                    std::cout << "  c[" << orderIdx << "]: " << ransacResult.bestFittingModelCoefficients[orderIdx] << std::endl;
                }
            }

            // Update the estimated instrument line shape
            if (this->settings.estimateInstrumentLineShape == InstrumentLineshapeEstimationOption::SuperGaussian)
            {
                EstimateInstrumentLineShapeAsSuperGaussian(result, fraunhoferSetup, ozoneSetup.get());
            }
            else if (this->settings.estimateInstrumentLineShape == InstrumentLineshapeEstimationOption::ApproximateGaussian)
            {
                EstimateInstrumentLineShapeAsApproximateGaussian(result, fraunhoferSetup);
            }

            // Adjust the selection parameter for the maximum error in the wavelength calibration such that the search space decreases for each iteration
            correspondenceSelectionSettings.maximumPixelDistanceForPossibleCorrespondence = std::max(10, correspondenceSelectionSettings.maximumPixelDistanceForPossibleCorrespondence / 2);

            if (iterationIdx < numberOfIterations - 1)
            {
                // Adjust the shape of the Fraunhofer spectrum using what we now know:
                //  We have the sensitivity of the spectrometer (pixel -> intesity) from the calculated envelope
                //  We have the pixel-to-wavelength mapping of the spectrometer from the ransac calibration
                // Re-convolve the Fraunhofer spectrum to get it on the new pixel grid and with the new instrument line shape (if updated)

                std::cout << " -- Adjusting Fraunhofer Spectrum -- " << std::endl;
                {
                    novac::CCrossSectionData& currentEstimateOfInstrumentLineShape = (result.estimatedInstrumentLineShape.GetSize() > 0) ? result.estimatedInstrumentLineShape : settings.initialInstrumentLineShape;
                    calibrationState.fraunhoferSpectrum = fraunhoferSetup.GetFraunhoferSpectrum(result.pixelToWavelengthMapping, currentEstimateOfInstrumentLineShape);
                }

                UpdateFraunhoferSpectrumWithApparentSensitivity(ransacResult);
            }
        }

        // Normalize the output, such that other programs may use the data directly.
        if (result.estimatedInstrumentLineShape.GetSize() > 0)
        {
            ::Normalize(result.estimatedInstrumentLineShape.m_crossSection);
        }

        return result;
    }

    void WavelengthCalibrationSetup::UpdateFraunhoferSpectrumWithApparentSensitivity(novac::RansacWavelengthCalibrationResult& ransacResult)
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
        for (size_t pixelIdx = 0; pixelIdx < static_cast<size_t>(calibrationState.fraunhoferSpectrum->m_length); ++pixelIdx)
        {
            const double apparentSensitivity = apparentSensitivitySpline.GetValue((MathFit::TFitData)calibrationState.fraunhoferSpectrum->m_wavelength[pixelIdx]);
            calibrationState.fraunhoferSpectrum->m_data[pixelIdx] *= apparentSensitivity;
        }
        Normalize(*calibrationState.fraunhoferSpectrum);
    }

    void WavelengthCalibrationSetup::EstimateInstrumentLineShapeAsApproximateGaussian(novac::SpectrometerCalibrationResult& result, novac::FraunhoferSpectrumGeneration& fraunhoferSetup)
    {
        InstrumentLineShapeEstimationFromKeypointDistance ilsEstimator{ result.pixelToWavelengthMapping };
        if (settings.initialInstrumentLineShape.GetSize() > 0)
        {
            ilsEstimator.UpdateInitialLineShape(settings.initialInstrumentLineShape);
        }

        // Create an estimation and update the 'estimatedInstrumentLineShape'
        double lineShapeFwhm = 0.0;
        auto estimationResult = ilsEstimator.EstimateInstrumentLineShape(fraunhoferSetup, *calibrationState.measuredSpectrum, result.estimatedInstrumentLineShape, lineShapeFwhm);

        // Save the result
        result.estimatedInstrumentLineShapeParameters = std::make_unique<novac::GaussianLineShape>(estimationResult.lineShape);
    }

    void WavelengthCalibrationSetup::EstimateInstrumentLineShapeAsSuperGaussian(
        novac::SpectrometerCalibrationResult& result,
        novac::FraunhoferSpectrumGeneration& fraunhoferSetup,
        novac::ICrossSectionSpectrumGenerator* ozoneSetup)
    {
        // The super-gaussian estimation requires that there is an initial estimate of the instrument line shape.
        //  If there isn't any, then create an initial guess using the approximate-gaussian approach.
        if (result.estimatedInstrumentLineShape.GetSize() == 0 && settings.initialInstrumentLineShape.GetSize() == 0)
        {
            EstimateInstrumentLineShapeAsApproximateGaussian(result, fraunhoferSetup);
        }

        try
        {
            novac::CCrossSectionData& currentEstimateOfInstrumentLineShape = (result.estimatedInstrumentLineShape.GetSize() > 0) ? result.estimatedInstrumentLineShape : settings.initialInstrumentLineShape;
            InstrumentLineshapeEstimationFromDoas ilsEstimator{ result.pixelToWavelengthMapping, currentEstimateOfInstrumentLineShape };

            InstrumentLineshapeEstimationFromDoas::LineShapeEstimationSettings estimationSettings;
            estimationSettings.startPixel = (size_t)std::round(novac::GetFractionalIndex(result.pixelToWavelengthMapping, settings.estimateInstrumentLineShapeWavelengthRegion.first));
            estimationSettings.endPixel = (size_t)std::round(novac::GetFractionalIndex(result.pixelToWavelengthMapping, settings.estimateInstrumentLineShapeWavelengthRegion.second));
            if (estimationSettings.startPixel > estimationSettings.endPixel)
            {
                std::swap(estimationSettings.startPixel, estimationSettings.endPixel);
            }

            auto startTime = std::chrono::steady_clock::now();
            auto estimationResult = ilsEstimator.EstimateInstrumentLineShape(*calibrationState.measuredSpectrum, estimationSettings, fraunhoferSetup, ozoneSetup);
            auto stopTime = std::chrono::steady_clock::now();

            // Output for debugging
            std::cout << "Instrument line shape estimation took " << std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime).count() << " ms" << std::endl;
            std::cout << " Result: (w: " << estimationResult.result.lineShape.w << ", k: " << estimationResult.result.lineShape.k << ") found in " << estimationResult.attempts.size() << " iterations." << std::endl;
            // Save the result.
            result.estimatedInstrumentLineShapeParameters = std::make_unique<novac::SuperGaussianLineShape>(estimationResult.result.lineShape);
            result.estimatedInstrumentLineShape = SampleInstrumentLineShape(estimationResult.result.lineShape);
            result.estimatedInstrumentLineShapeError = estimationResult.result.error;
            result.estimatedInstrumentLineShapeShift = estimationResult.result.shift;
            result.estimatedInstrumentLineShapePixelRange.first = estimationSettings.startPixel;
            result.estimatedInstrumentLineShapePixelRange.second = estimationSettings.endPixel;
        }
        catch (std::exception& e)
        {
            std::cout << "Instrument line shape estimation failed: " << e.what() << std::endl;
        }
    }

}