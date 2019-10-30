#pragma once

#include <vector>
#include <string>

// ---------------------------------------------------------------------------------------------------------------
// ----------- This header contains methods used to perform wavelength calibration of a spectrometer -------------
// ---------------------------------------------------------------------------------------------------------------

namespace Evaluation
{
    class CCrossSectionData;

    struct SpectrometerCalibration
    {
        /** The wavelength for each pixel on the detector */
        std::vector<double> wavelengthToPixelMapping;

        /** The estimated slit function of the spectrometer */
        std::vector<double> slf;
    };

    /** Estimates the wavelength to pixel mapping for a given measured spectrum by fitting  */
    bool EstimateWavelengthToPixelMapping(const std::string& measuredspectrum, const std::string& initialWavelengthToPixelMapping, const std::string& solarAtlas, SpectrometerCalibration& result);

    /** Estimates the wavelength to pixel mapping for a given measured spectrum by fitting  */
    bool EstimateWavelengthToPixelMapping(const CCrossSectionData& measuredspectrum, const CCrossSectionData& initialWavelengthToPixelMapping, const CCrossSectionData& solarAtlas, SpectrometerCalibration& result);

}