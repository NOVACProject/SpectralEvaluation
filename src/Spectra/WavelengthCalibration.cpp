#include <algorithm>
#include <SpectralEvaluation/Spectra/WavelengthCalibration.h>
#include <SpectralEvaluation/Spectra/ReferenceSpectrumConvolution.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/VectorUtils.h>
#include <SpectralEvaluation/File/File.h>

namespace Evaluation
{
    bool EstimateWavelengthToPixelMapping(const std::string& measuredSpectrumFile, const std::string& initialWavelengthToPixelMappingFile, const std::string& solarAtlasFile, SpectrometerCalibration& result)
    {
        CCrossSectionData measuredSpectrum;
        if (!FileIo::ReadCrossSectionFile(measuredSpectrumFile, measuredSpectrum))
        {
            return false;
        }

        CCrossSectionData pixelToWavelengthMapping;
        if (!FileIo::ReadCrossSectionFile(initialWavelengthToPixelMappingFile, pixelToWavelengthMapping, true))
        {
            return false;
        }
        if (pixelToWavelengthMapping.m_waveLength.size() != measuredSpectrum.m_crossSection.size())
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

    bool EstimateWavelengthToPixelMapping(const CCrossSectionData& measuredspectrum, const CCrossSectionData& initialWavelengthToPixelMapping, const CCrossSectionData& solarAtlas, SpectrometerCalibration& result)
    {
        CCrossSectionData convolvedReference;
        CCrossSectionData pixelToWavelengthMapping{ initialWavelengthToPixelMapping };

        double gaussianSigma = 0.7; // [nm] good initial guess for Novac Instruments.

        // Generate a gaussian with high enough resolution
        CCrossSectionData slf;
        const double slfDeltaLambda = std::min(0.25 * Resolution(solarAtlas.m_waveLength), gaussianSigma * 0.125);
        CreateGaussian(gaussianSigma, slfDeltaLambda, slf);

        if (!ConvolveReference(pixelToWavelengthMapping.m_waveLength, slf, solarAtlas, result.wavelengthToPixelMapping))
        {
            return false;
        }


        return true;
    }
}
