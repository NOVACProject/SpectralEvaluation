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

// TODO: Move to separate file, with FraunhoferSpectrumGeneration
#include <SpectralEvaluation/Evaluation/EvaluationBase.h>
#include <SpectralEvaluation/Spectra/Scattering.h>

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

// --------------------------- FraunhoferSpectrumGeneration ---------------------------

std::unique_ptr<CSpectrum> FraunhoferSpectrumGeneration::GetFraunhoferSpectrum(
    const std::vector<double>& pixelToWavelengthMapping,
    const CCrossSectionData& measuredInstrumentLineShape)
{
    return GetFraunhoferSpectrum(pixelToWavelengthMapping, measuredInstrumentLineShape, this->crossSectionsToInclude);
}

std::unique_ptr<CSpectrum> FraunhoferSpectrumGeneration::GetFraunhoferSpectrum(
    const std::vector<double>& pixelToWavelengthMapping,
    const CCrossSectionData& measuredInstrumentLineShape,
    std::vector<AbsorbingCrossSection>& localCrossSectionsToInclude)
{
    ReadSolarCrossSection();

    // Create a local copy which we can scale as we want.
    CCrossSectionData localSolarCrossSection{ *this->solarCrossSection };

    for (auto& absorber : localCrossSectionsToInclude)
    {
        // Turn the molecular absorption into an absorbance spectrum and multiply with the high res solar
        if (std::abs(absorber.totalColumn) > std::numeric_limits<double>::epsilon())
        {
            // Get the high res cross section
            if (absorber.crossSectionData == nullptr)
            {
                absorber.crossSectionData = std::make_unique<CCrossSectionData>();
                ReadCrossSectionFile(absorber.path, *absorber.crossSectionData);
            }

            // Create a local copy which we can scale as we want.
            CCrossSectionData crossSectionCopy{ *absorber.crossSectionData };
            Mult(crossSectionCopy.m_crossSection, -absorber.totalColumn);
            Exp(crossSectionCopy.m_crossSection);
            std::vector<double> resampledCrossSection;
            Resample(crossSectionCopy, localSolarCrossSection.m_waveLength, resampledCrossSection);
            Mult(resampledCrossSection, localSolarCrossSection.m_crossSection);
        }
    }

    // Generate a theoretical solar spectrum by convolving the high-res solar atlas with the measured slf
    auto startTime = std::chrono::steady_clock::now();
    std::vector<double> theoreticalFraunhoferSpectrumData;
    ConvolveReference(pixelToWavelengthMapping, measuredInstrumentLineShape, localSolarCrossSection, theoreticalFraunhoferSpectrumData, WavelengthConversion::None, ConvolutionMethod::Fft);
    auto stopTime = std::chrono::steady_clock::now();
    std::cout << "Convolution of Fraunhofer Reference took " << std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime).count() << " ms" << std::endl;

    std::unique_ptr<CSpectrum> theoreticalFraunhoferSpectrum = std::make_unique<CSpectrum>(pixelToWavelengthMapping, theoreticalFraunhoferSpectrumData);
    Normalize(*theoreticalFraunhoferSpectrum); // normalizes the intensity

    return theoreticalFraunhoferSpectrum;
}

std::unique_ptr<CSpectrum> FraunhoferSpectrumGeneration::GetFraunhoferSpectrumMatching(
    const std::vector<double>& pixelToWavelengthMapping,
    const novac::CSpectrum& measuredSpectrum,
    const novac::CCrossSectionData& measuredInstrumentLineShape)
{
    if (this->crossSectionsToInclude.size() == 0)
    {
        throw std::invalid_argument("Cannot create a Fraunhofer spectrum with matching absorption if no absorption cross sections have been provided.");
    }
    if (pixelToWavelengthMapping.size() == 0 || pixelToWavelengthMapping.size() != static_cast<size_t>(measuredSpectrum.m_length))
    {
        throw std::invalid_argument("Cannot create a Fraunhofer spectrum without a pixel-to-wavelength mapping.");
    }

    ReadSolarCrossSection();

    // Create a local copy which we can scale as we want.
    CCrossSectionData localSolarCrossSection{ *this->solarCrossSection };

    // 1. Create a Fraunhofer spectrum without any absorption cross sections
    CCrossSectionData initialFraunhoferSpectrum;
    initialFraunhoferSpectrum.m_waveLength = pixelToWavelengthMapping;
    ConvolveReference(initialFraunhoferSpectrum.m_waveLength, measuredInstrumentLineShape, localSolarCrossSection, initialFraunhoferSpectrum.m_crossSection, WavelengthConversion::None, ConvolutionMethod::Fft);

    // DEBUG
    novac::SaveCrossSectionFile("D:/NOVAC/SpectrometerCalibration/Output/Convolution_Meas.xs", measuredSpectrum);
    novac::SaveCrossSectionFile("D:/NOVAC/SpectrometerCalibration/Output/Convolution_Solar.xs", initialFraunhoferSpectrum.m_crossSection);

    // 2. Setup a DOAS fit to match the newly created fraunhofer spectrum to the measured
    {
        novac::CFitWindow fitWindow;
        fitWindow.fitLow = static_cast<int>(std::round(initialFraunhoferSpectrum.FindWavelength(310.0)));  // TODO: Determine!
        fitWindow.fitHigh = static_cast<int>(std::round(initialFraunhoferSpectrum.FindWavelength(330.0))); // TODO: Determine!
        fitWindow.fitType = novac::FIT_POLY;
        fitWindow.polyOrder = 3;
        fitWindow.specLength = measuredSpectrum.m_length;
        fitWindow.ringCalculation = novac::RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING;
        for (size_t ii = 0; ii < this->crossSectionsToInclude.size(); ++ii)
        {
            if (this->crossSectionsToInclude[ii].crossSectionData == nullptr)
            {
                this->crossSectionsToInclude[ii].crossSectionData = std::make_unique<CCrossSectionData>();
                ReadCrossSectionFile(this->crossSectionsToInclude[ii].path, *this->crossSectionsToInclude[ii].crossSectionData);
            }

            fitWindow.ref[ii].m_data = std::make_unique<novac::CCrossSectionData>();
            fitWindow.ref[ii].m_data->m_waveLength = initialFraunhoferSpectrum.m_waveLength;

            // Convolve to create a usable reference.
            std::vector<double> lowResCrossSection;
            ConvolveReference(
                initialFraunhoferSpectrum.m_waveLength,
                measuredInstrumentLineShape,
                *this->crossSectionsToInclude[ii].crossSectionData,
                fitWindow.ref[ii].m_data->m_crossSection,
                WavelengthConversion::None,
                ConvolutionMethod::Fft);

            // DEBUG
            {
                std::stringstream filename;
                filename << "D:/NOVAC/SpectrometerCalibration/Output/Convolution_Ref" << ii << ".xs";
                novac::SaveCrossSectionFile(filename.str(), *(fitWindow.ref[ii].m_data));
            }

            if (ii == 0)
            {
                fitWindow.ref[ii].m_shiftOption = novac::SHIFT_FREE;
            }
            else
            {
                fitWindow.ref[ii].m_shiftOption = novac::SHIFT_LINK;
                fitWindow.ref[ii].m_shiftValue = 0; // i.e. link to #0
            }
            fitWindow.ref[ii].m_squeezeOption = novac::SHIFT_FIX;
            fitWindow.ref[ii].m_squeezeValue = 1.0; // fix the squeeze to 1.0 to avoid a too complex fit.
        }

        // Also calculate a ring-spectrum and include
        std::unique_ptr<CCrossSectionData> ring;
        {
            novac::CSpectrum modifiableSpectrum{ measuredSpectrum };
            modifiableSpectrum.m_wavelength = pixelToWavelengthMapping;
            Doasis::Scattering scatteringCalc;
            auto ringSpectrum = scatteringCalc.CalcRingSpectrum(modifiableSpectrum);

            // Convert the CSpectrum to a CCrossSectionData
            size_t idx = this->crossSectionsToInclude.size();

            ring = std::make_unique<novac::CCrossSectionData>();
            ring->m_crossSection = std::vector<double>(ringSpectrum.m_data, ringSpectrum.m_data + ringSpectrum.m_length);
            ring->m_waveLength = initialFraunhoferSpectrum.m_waveLength;

            fitWindow.ref[idx].m_data = std::make_unique<novac::CCrossSectionData>(*ring);
            fitWindow.ref[idx].m_shiftOption = novac::SHIFT_LINK;
            fitWindow.ref[idx].m_shiftValue = 0; // i.e. link to #0
            fitWindow.ref[idx].m_squeezeOption = novac::SHIFT_FIX;
            fitWindow.ref[idx].m_squeezeValue = 1.0; // fix the squeeze to 1.0 to avoid a too complex fit.
        }

        // And finally include the sky spectrum into the fit
        {
            size_t idx = this->crossSectionsToInclude.size() + 1;
            CBasicMath math;

            fitWindow.ref[idx].m_data = std::make_unique<novac::CCrossSectionData>(initialFraunhoferSpectrum);
            math.Log(fitWindow.ref[idx].m_data->m_crossSection.data(), static_cast<int>(fitWindow.ref[idx].m_data->m_crossSection.size()));
            fitWindow.ref[idx].m_columnOption = novac::SHIFT_FIX;
            fitWindow.ref[idx].m_columnValue = -1;
            fitWindow.ref[idx].m_shiftOption = novac::SHIFT_LINK;
            fitWindow.ref[idx].m_shiftValue = 0; // i.e. link to #0
            fitWindow.ref[idx].m_squeezeOption = novac::SHIFT_FIX;
            fitWindow.ref[idx].m_squeezeValue = 1.0; // fix the squeeze to 1.0 to avoid a too complex fit.
        }


        // DEBUG
        {
            novac::SaveCrossSectionFile("D:/NOVAC/SpectrometerCalibration/Output/Convolution_Ring.xs", *(ring));
        }

        fitWindow.nRef = 2 + static_cast<int>(this->crossSectionsToInclude.size());

        novac::CEvaluationBase eval;
        eval.SetFitWindow(fitWindow);
        // eval.SetSkySpectrum(initialFraunhoferSpectrum);
        if (eval.Evaluate(measuredSpectrum))
        {
            throw std::exception();
        }

        // we now have the desired column value for each reference
        std::vector<AbsorbingCrossSection> absorbingCrossSection;

        for (size_t ii = 0; ii < this->crossSectionsToInclude.size(); ++ii)
        {
            AbsorbingCrossSection xs;
            xs.path = this->crossSectionsToInclude[ii].path;
            xs.totalColumn = eval.m_result.m_referenceResult[ii].m_column;
            absorbingCrossSection.push_back(std::move(xs));
        }
        // also include the Ring spectrum
        {
            AbsorbingCrossSection xs;
            xs.path = "";
            xs.crossSectionData = std::move(ring);
            xs.totalColumn = eval.m_result.m_referenceResult[this->crossSectionsToInclude.size()].m_column;
            absorbingCrossSection.push_back(std::move(xs));
        }

        std::cout << "Determined total column of absorbers in measured spectrum using a DOAS fit. The results are:" << std::endl;
        std::cout << " Chi2:   " << eval.m_result.m_chiSquare << std::endl;
        std::cout << " Lo chn: " << fitWindow.fitLow << std::endl;
        std::cout << " Hi chn: " << fitWindow.fitHigh << std::endl;
        for (size_t ii = 0; ii < fitWindow.nRef; ++ii)
        {
            std::cout << " Ref" << ii << std::endl;
            std::cout << "   Column   " << eval.m_result.m_referenceResult[ii].m_column << " +- " << eval.m_result.m_referenceResult[ii].m_columnError << std::endl;
            std::cout << "   Shift    " << eval.m_result.m_referenceResult[ii].m_shift << " +- " << eval.m_result.m_referenceResult[ii].m_shiftError << std::endl;
            std::cout << "   Squeeze  " << eval.m_result.m_referenceResult[ii].m_squeeze << " +- " << eval.m_result.m_referenceResult[ii].m_squeezeError << std::endl;
        }

        return GetFraunhoferSpectrum(pixelToWavelengthMapping, measuredInstrumentLineShape, absorbingCrossSection);
    }
}

void FraunhoferSpectrumGeneration::ReadSolarCrossSection()
{
    if (this->solarCrossSection == nullptr)
    {
        this->solarCrossSection = std::make_unique<CCrossSectionData>();
        ReadCrossSectionFile(this->solarAtlasFile, *solarCrossSection);
    }
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


WavelengthCalibrationSetup::WavelengthCalibrationSetup(const WavelengthCalibrationSettings& calibrationSettings)
    : settings(calibrationSettings)
{
}

SpectrometerCalibrationResult WavelengthCalibrationSetup::DoWavelengthCalibration(const CSpectrum& measuredSpectrum, const CCrossSectionData& measuredInstrumentLineShape)
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
    // calibrationState.originalFraunhoferSpectrum = fraunhoferSetup.GetFraunhoferSpectrum(settings.initialPixelToWavelengthMapping, measuredInstrumentLineShape);
    calibrationState.originalFraunhoferSpectrum = fraunhoferSetup.GetFraunhoferSpectrumMatching(settings.initialPixelToWavelengthMapping, *calibrationState.measuredSpectrum, measuredInstrumentLineShape);
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
            // Adjust the selection parameter for the maximum error in the wavelength calibration
            correspondenceSelectionSettings.maximumPixelDistanceForPossibleCorrespondence = std::max(10, correspondenceSelectionSettings.maximumPixelDistanceForPossibleCorrespondence / 2);

            // Re-convolve the Fraunhofer spectrum to get it on the new pixel grid
            // calibrationState.fraunhoferSpectrum = fraunhoferSetup.GetFraunhoferSpectrum(result.pixelToWavelengthMapping, measuredInstrumentLineShape);
            calibrationState.fraunhoferSpectrum = fraunhoferSetup.GetFraunhoferSpectrumMatching(result.pixelToWavelengthMapping, *calibrationState.measuredSpectrum, measuredInstrumentLineShape);

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