#pragma once

#include <vector>
#include <string>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>

// ---------------------------------------------------------------------------------------------------------------
// ----------- This header contains methods used to perform wavelength calibration of a spectrometer -------------
// ---------------------------------------------------------------------------------------------------------------

class CSpectrum;

namespace Evaluation
{
    struct SpectrometerCalibration
    {
        /** The wavelength for each pixel on the detector */
        std::vector<double> wavelengthToPixelMapping;

        /** The estimated slit function of the spectrometer */
        CCrossSectionData slf;
    };

    struct WavelengthCalibrationSetup
    {
        // A high-resolved Kurucz spektrum
        CCrossSectionData solarAtlas;

        // The high-resolution absorption cross sections necessary to get a good fit.
        // Typically just O3 (the ring spectrum will be calculated in the fit routine).
        std::vector<CCrossSectionData> crossSections;
    };

    /** Estimates the wavelength to pixel mapping for a given measured spectrum by fitting
        the solar atlas (convolved with the instrument slit function) towards the measured spectrum.
        This will estimate the wavelength-to-pixel mapping but not alter the slit function.
        @param measuredspectrum A (dark-corrected) measured spectrum.
        @param initialCalibration The initial wavelength-to-pixel mapping and the measured / estimated slit-function. */
    bool EstimateWavelengthToPixelMapping(
        const WavelengthCalibrationSetup& calibrationSetup,
        const SpectrometerCalibration& initialCalibration,
        const CSpectrum& measuredspectrum,
        SpectrometerCalibration& result);

}