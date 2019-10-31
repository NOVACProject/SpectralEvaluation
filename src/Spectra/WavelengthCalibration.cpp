#include <algorithm>
#include <SpectralEvaluation/Spectra/WavelengthCalibration.h>
#include <SpectralEvaluation/Spectra/ReferenceSpectrumConvolution.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/VectorUtils.h>
#include <SpectralEvaluation/File/File.h>

namespace Evaluation
{
    bool EstimateWavelengthToPixelMapping(
        const std::string& measuredSpectrumFile,
        const std::string& initialWavelengthToPixelMappingFile,
        const std::string& solarAtlasFile,
        SpectrometerCalibration& result)
    {
        CSpectrum measuredSpectrum;
        if (!FileIo::ReadSpectrum(measuredSpectrumFile, measuredSpectrum))
        {
            return false;
        }

        CCrossSectionData pixelToWavelengthMapping;
        if (!FileIo::ReadCrossSectionFile(initialWavelengthToPixelMappingFile, pixelToWavelengthMapping, true))
        {
            return false;
        }
        if (pixelToWavelengthMapping.m_waveLength.size() != (size_t)measuredSpectrum.m_length)
        {
            return false; // not same length
        }

        CCrossSectionData solarAtlas;
        if (!FileIo::ReadCrossSectionFile(solarAtlasFile, solarAtlas))
        {
            return false;
        }
        if (solarAtlas.m_waveLength.size() == 0 || solarAtlas.m_crossSection.size() == 0)
        {
            return false;
        }

        return EstimateWavelengthToPixelMapping(measuredSpectrum, pixelToWavelengthMapping, solarAtlas, result);
    }

    bool EstimateWavelengthToPixelMapping(const CSpectrum& /*measuredspectrum*/, const CCrossSectionData& initialWavelengthToPixelMapping, const CCrossSectionData& solarAtlas, SpectrometerCalibration& result)
    {
        CCrossSectionData convolvedReference;
        CCrossSectionData pixelToWavelengthMapping{ initialWavelengthToPixelMapping };

        double gaussianSigma = 0.7; // [nm] good initial guess for Novac Instruments.

        // TODO: we need a way to pass in the wavelength range as parameter here
        double lambdaLow = 310.0;
        double lambdaHigh = 320.0;
        CCrossSectionData localSolarAtlas{solarAtlas, lambdaLow, lambdaHigh};
        CCrossSectionData localPixelToWavelenthMap{pixelToWavelengthMapping, lambdaLow, lambdaHigh};

        // Generate a gaussian with high enough resolution
        CCrossSectionData slf;
        const double slfDeltaLambda = std::min(0.25 * Resolution(localSolarAtlas.m_waveLength), gaussianSigma * 0.125);
        CreateGaussian(gaussianSigma, slfDeltaLambda, slf);

        if (!ConvolveReference(localPixelToWavelenthMap.m_waveLength, slf, localSolarAtlas, result.wavelengthToPixelMapping))
        {
            return false;
        }


        return true;
    }
}
