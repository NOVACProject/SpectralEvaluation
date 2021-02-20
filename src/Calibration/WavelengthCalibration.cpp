#include <algorithm>
#include <chrono>
#include <SpectralEvaluation/Calibration/WavelengthCalibration.h>
#include <SpectralEvaluation/Calibration/ReferenceSpectrumConvolution.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Evaluation/WavelengthFit.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
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
    const ::Evaluation::CCrossSectionData& measuredSlf,
    double ozoneTotalColumn)
{
    // Get the high res solar spectrum
    Evaluation::CCrossSectionData solarCrossSection;
    FileIo::ReadCrossSectionFile(this->solarAtlasFile, solarCrossSection);

    // Turn the ozone into an absorbance spectrum and multiply with the high res solar
    if (std::abs(ozoneTotalColumn) > std::numeric_limits<double>::epsilon())
    {
        // Get the high res ozone spectrum
        Evaluation::CCrossSectionData ozoneCrossSection;
        FileIo::ReadCrossSectionFile(this->ozoneCrossSectionFile, ozoneCrossSection);

        Mult(ozoneCrossSection.m_crossSection, -ozoneTotalColumn);
        Exp(ozoneCrossSection.m_crossSection);
        std::vector<double> resampledOzoneCrossSection;
        Evaluation::Resample(ozoneCrossSection, solarCrossSection.m_waveLength, resampledOzoneCrossSection);
        Mult(resampledOzoneCrossSection, solarCrossSection.m_crossSection);
    }

    // Generate a theoretical solar spectrum by convolving the high-res solar atlas with the measured slf
    auto startTime = std::chrono::steady_clock::now();
    std::vector<double> theoreticalFraunhoferSpectrumData;
    ::Evaluation::ConvolveReference(pixelToWavelengthMapping, measuredSlf, solarCrossSection, theoreticalFraunhoferSpectrumData, WavelengthConversion::None, ConvolutionMethod::Fft);
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


bool EstimateWavelengthToPixelMapping(
    const WavelengthCalibrationSetup& /*calibrationSetup*/,
    const SpectrometerCalibration& /*initialCalibration*/,
    const CSpectrum& /*measuredspectrum*/,
    SpectrometerCalibration& /*result*/)
{



    return true;
}

}