#include <SpectralEvaluation/Calibration/CrossSectionSpectrumGenerator.h>
#include <SpectralEvaluation/Calibration/ReferenceSpectrumConvolution.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Spectra/WavelengthRange.h>
#include <stdexcept>
#include <algorithm>

using namespace novac;

WavelengthRange CrossSectionSpectrumGenerator::GetSpectrumRange(const std::vector<double>& pixelToWavelengthMapping)
{
    ReadCrossSection();

    const WavelengthRange resultingRange(
        std::max(pixelToWavelengthMapping.front(), m_highResolutionCrossSection->m_waveLength.front()),
        std::min(pixelToWavelengthMapping.back(), m_highResolutionCrossSection->m_waveLength.back()));

    return resultingRange;
}

std::unique_ptr<CSpectrum> CrossSectionSpectrumGenerator::GetCrossSection(
    const std::vector<double>& pixelToWavelengthMapping,
    const novac::CCrossSectionData& measuredInstrumentLineShape)
{
    return GetCrossSection(pixelToWavelengthMapping, measuredInstrumentLineShape, 0.0, true);
}

std::unique_ptr<CSpectrum> CrossSectionSpectrumGenerator::GetCrossSection(
    const std::vector<double>& pixelToWavelengthMapping,
    const novac::CCrossSectionData& measuredInstrumentLineShape,
    double fwhmOfInstrumentLineShape,
    bool normalize)
{
    ReadCrossSection();

    // Generate a theoretical solar spectrum by convolving the high-res solar atlas with the measured slf
    std::vector<double> convolvedReferenceSpectrumData;
    ConvolveReference(
        pixelToWavelengthMapping,
        measuredInstrumentLineShape,
        *m_highResolutionCrossSection,
        convolvedReferenceSpectrumData,
        WavelengthConversion::None,
        ConvolutionMethod::Fft,
        fwhmOfInstrumentLineShape,
        normalize);

    std::unique_ptr<CSpectrum> convolvedReferenceSpectrum = std::make_unique<CSpectrum>(pixelToWavelengthMapping, convolvedReferenceSpectrumData);
    if (normalize)
    {
        Normalize(*convolvedReferenceSpectrum); // normalizes the intensity
    }

    return convolvedReferenceSpectrum;
}

void CrossSectionSpectrumGenerator::ReadCrossSection()
{
    if (m_highResolutionCrossSection == nullptr)
    {
        if (m_crossSectionFile.size() == 0)
        {
            throw std::invalid_argument("Missing input: a solar atlas file needs to be passed to Fraunhofer spectrum generation.");
        }

        m_highResolutionCrossSection = std::make_unique<CCrossSectionData>();
        m_highResolutionCrossSection->ReadCrossSectionFile(m_crossSectionFile);

        if (m_highResolutionCrossSection->m_crossSection.size() == 0)
        {
            throw std::invalid_argument("Invalid solar atlas file passed to Fraunhofer spectrum generation, file could not be read.");
        }
    }
}

