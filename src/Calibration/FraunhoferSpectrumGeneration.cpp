#include <SpectralEvaluation/Calibration/FraunhoferSpectrumGeneration.h>
#include <SpectralEvaluation/Calibration/ReferenceSpectrumConvolution.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Spectra/WavelengthRange.h>
#include <SpectralEvaluation/VectorUtils.h>
#include <chrono>
#include <iostream>
#include <algorithm>

#undef min
#undef max

namespace novac
{

    WavelengthRange FraunhoferSpectrumGeneration::GetFraunhoferRange(const std::vector<double>& pixelToWavelengthMapping)
    {
        ReadSolarCrossSection();

        const WavelengthRange resultingRange(
            std::max(pixelToWavelengthMapping.front(), this->solarCrossSection->m_waveLength.front()),
            std::min(pixelToWavelengthMapping.back(), this->solarCrossSection->m_waveLength.back()));

        return resultingRange;
    }

    std::unique_ptr<CSpectrum> FraunhoferSpectrumGeneration::GetFraunhoferSpectrum(
        const std::vector<double>& pixelToWavelengthMapping,
        const CCrossSectionData& measuredInstrumentLineShape)
    {
        return GetFraunhoferSpectrum(pixelToWavelengthMapping, measuredInstrumentLineShape, this->crossSectionsToInclude, 0.0, true);
    }

    std::unique_ptr<CSpectrum> FraunhoferSpectrumGeneration::GetFraunhoferSpectrum(
        const std::vector<double>& pixelToWavelengthMapping,
        const CCrossSectionData& measuredInstrumentLineShape,
        double fwhmOfInstrumentLineShape,
        bool normalize)
    {
        return GetFraunhoferSpectrum(pixelToWavelengthMapping, measuredInstrumentLineShape, this->crossSectionsToInclude, fwhmOfInstrumentLineShape, normalize);
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
        if (debugOutput)
        {
            std::cout << "Convolution of Fraunhofer Reference took " << std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime).count() << " ms" << std::endl;
        }

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

        if (debugOutput)
        {
            std::cout << "Convolution of differential Fraunhofer Reference took " << std::chrono::duration_cast<std::chrono::milliseconds>(stopTime - startTime).count() << " ms" << std::endl;
        }

        std::unique_ptr<CSpectrum> theoreticalFraunhoferSpectrum = std::make_unique<CSpectrum>(pixelToWavelengthMapping, theoreticalFraunhoferSpectrumData);

        // Make a normalization of the intensity, without changing the mean.
        // const double range = theoreticalFraunhoferSpectrum->MaxValue() - theoreticalFraunhoferSpectrum->MinValue();
        // for (long ii = 0; ii < theoreticalFraunhoferSpectrum->m_length; ++ii)
        // {
        //     theoreticalFraunhoferSpectrum->m_data[ii] = theoreticalFraunhoferSpectrum->m_data[ii] / range;
        // }

        return theoreticalFraunhoferSpectrum;
    }

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