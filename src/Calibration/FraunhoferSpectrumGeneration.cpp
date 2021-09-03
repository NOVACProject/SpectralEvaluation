#include <SpectralEvaluation/Calibration/FraunhoferSpectrumGeneration.h>
#include <SpectralEvaluation/Calibration/ReferenceSpectrumConvolution.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/VectorUtils.h>

#include <chrono>
#include <iostream>

// Used for the experimentation below
#ifdef USE_DOAS_FIT
#include <SpectralEvaluation/Evaluation/EvaluationBase.h>
#include <SpectralEvaluation/Spectra/Scattering.h>
#endif

namespace novac
{

std::unique_ptr<CSpectrum> FraunhoferSpectrumGeneration::GetFraunhoferSpectrum(
    const std::vector<double>& pixelToWavelengthMapping,
    const CCrossSectionData& measuredInstrumentLineShape)
{
    return GetFraunhoferSpectrum(pixelToWavelengthMapping, measuredInstrumentLineShape, this->crossSectionsToInclude, 0.0, true);
}

std::unique_ptr<CSpectrum> FraunhoferSpectrumGeneration::GetFraunhoferSpectrum(
    const std::vector<double>& pixelToWavelengthMapping,
    const CCrossSectionData& measuredInstrumentLineShape,
    bool normalize)
{
    return GetFraunhoferSpectrum(pixelToWavelengthMapping, measuredInstrumentLineShape, this->crossSectionsToInclude, 0.0, normalize);
}

std::unique_ptr<CSpectrum> FraunhoferSpectrumGeneration::GetDifferentialFraunhoferSpectrum(
    const std::vector<double>& pixelToWavelengthMapping,
    const CCrossSectionData& measuredInstrumentLineShape,
    double fwhmOfInstrumentLineShape)
{
    return GetDifferentialFraunhoferSpectrum(pixelToWavelengthMapping, measuredInstrumentLineShape, this->crossSectionsToInclude, fwhmOfInstrumentLineShape);
}

std::unique_ptr<CSpectrum> FraunhoferSpectrumGeneration::GetFraunhoferSpectrum(
    const std::vector<double>& pixelToWavelengthMapping,
    const CCrossSectionData& measuredInstrumentLineShape,
    std::vector<AbsorbingCrossSection>& localCrossSectionsToInclude,
    double fwhmOfInstrumentLineShape,
    bool normalize)
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
                absorber.crossSectionData->ReadCrossSectionFile(absorber.path);
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

    const bool normalizeInstrumentLineShape = true;

    // Generate a theoretical solar spectrum by convolving the high-res solar atlas with the measured slf
    auto startTime = std::chrono::steady_clock::now();
    std::vector<double> theoreticalFraunhoferSpectrumData;
    ConvolveReference(
        pixelToWavelengthMapping,
        measuredInstrumentLineShape,
        localSolarCrossSection,
        theoreticalFraunhoferSpectrumData,
        WavelengthConversion::None,
        ConvolutionMethod::Fft,
        fwhmOfInstrumentLineShape,
        normalizeInstrumentLineShape);

    auto stopTime = std::chrono::steady_clock::now();
    std::cout << "Convolution of Fraunhofer Reference took " << std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime).count() << " ms" << std::endl;

    std::unique_ptr<CSpectrum> theoreticalFraunhoferSpectrum = std::make_unique<CSpectrum>(pixelToWavelengthMapping, theoreticalFraunhoferSpectrumData);
    if (normalize)
    {
        Normalize(*theoreticalFraunhoferSpectrum); // normalizes the intensity
    }

    return theoreticalFraunhoferSpectrum;
}

std::unique_ptr<CSpectrum> FraunhoferSpectrumGeneration::GetDifferentialFraunhoferSpectrum(
    const std::vector<double>& pixelToWavelengthMapping,
    const CCrossSectionData& measuredInstrumentLineShape,
    std::vector<AbsorbingCrossSection>& localCrossSectionsToInclude,
    double fwhmOfInstrumentLineShape)
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
                absorber.crossSectionData->ReadCrossSectionFile(absorber.path);
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

    const bool normalizeInstrumentLineShape = false;

    // Generate a theoretical solar spectrum by convolving the high-res solar atlas with the measured slf
    auto startTime = std::chrono::steady_clock::now();
    std::vector<double> theoreticalFraunhoferSpectrumData;
    ConvolveReference(
        pixelToWavelengthMapping,
        measuredInstrumentLineShape,
        localSolarCrossSection,
        theoreticalFraunhoferSpectrumData,
        WavelengthConversion::None,
        ConvolutionMethod::Direct,
        fwhmOfInstrumentLineShape,
        normalizeInstrumentLineShape);

    auto stopTime = std::chrono::steady_clock::now();
    std::cout << "Convolution of differential Fraunhofer Reference took " << std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime).count() << " ms" << std::endl;

    std::unique_ptr<CSpectrum> theoreticalFraunhoferSpectrum = std::make_unique<CSpectrum>(pixelToWavelengthMapping, theoreticalFraunhoferSpectrumData);

    // Make a normalization of the intensity, without changing the mean.
    // const double range = theoreticalFraunhoferSpectrum->MaxValue() - theoreticalFraunhoferSpectrum->MinValue();
    // for (long ii = 0; ii < theoreticalFraunhoferSpectrum->m_length; ++ii)
    // {
    //     theoreticalFraunhoferSpectrum->m_data[ii] = theoreticalFraunhoferSpectrum->m_data[ii] / range;
    // }

    return theoreticalFraunhoferSpectrum;
}

#ifdef USE_DOAS_FIT
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

#endif

void FraunhoferSpectrumGeneration::ReadSolarCrossSection()
{
    if (this->solarCrossSection == nullptr)
    {
        if (this->solarAtlasFile.size() == 0)
        {
            throw std::invalid_argument("Missing input: a solar atlas file needs to be passed to Fraunhofer spectrum generation.");
        }

        this->solarCrossSection = std::make_unique<CCrossSectionData>();
        this->solarCrossSection->ReadCrossSectionFile(this->solarAtlasFile);

        if (this->solarCrossSection->m_crossSection.size() == 0)
        {
            throw std::invalid_argument("Invalid solar atlas file passed to Fraunhofer spectrum generation, file could not be read.");
        }
    }
}

}