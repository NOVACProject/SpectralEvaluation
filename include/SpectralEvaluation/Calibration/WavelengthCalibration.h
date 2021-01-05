#pragma once

#include <vector>
#include <string>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>

// ---------------------------------------------------------------------------------------------------------------
// ----------- This header contains methods used to perform wavelength calibration of a spectrometer -------------
// ---------------------------------------------------------------------------------------------------------------

class CSpectrum;
struct SpectrumDataPoint;

namespace Calibration
{
struct SpectrometerCalibration
{
    /** The wavelength for each pixel on the detector */
    std::vector<double> wavelengthToPixelMapping;

    /** The estimated slit function of the spectrometer */
    ::Evaluation::CCrossSectionData slf;
};

struct WavelengthCalibrationSetup
{
    // A high-resolved Kurucz spektrum
    ::Evaluation::CCrossSectionData solarAtlas;

    // The high-resolution absorption cross sections necessary to get a good fit.
    // Typically just O3 (the ring spectrum will be calculated in the fit routine).
    std::vector<::Evaluation::CCrossSectionData> crossSections;
};

/**
 * @brief Performs a wavelength calibration of a spectrometer using a measured mercury spectrum.
 *  This will identify the peaks in the measured spectrum and fit a polynomial to them such that the 
 *  pixel-to-wavelength mapping for the entire spectrum can be calculated.
 * @param measuredMercurySpectrum The measured spectrum
 * @param polynomialOrder The order of the polynomial to fit.
 * @param minimumWavelength The initial guess for the shortest wavelength present in the spectrum.
 *  This will be used to constrain the identification of the mercury lines.
 * @param maximumWavelength The initial guess for the longest wavelenght present in the spectrum.
 *  This will be used to constrain the identification of the mercury lines.
 * @param foundPeaks Will on successful return be filled with the location and wavelength of the
 *  successfully identified mercury lines.
 * @param pixelToWavelengthPolynomial The final pixel-to-wavelength mapping polynomial.
 *  This will have (polynomialOrder + 1) coefficients and store the 0:th order coefficient first.
 * @return True if the calibration was successful.
*/
bool MercuryCalibration(
    const CSpectrum& measuredMercurySpectrum,
    int polynomialOrder,
    double minimumWavelength,
    double maximumWavelength,
    std::vector<SpectrumDataPoint>& foundPeaks,
    std::vector<double> pixelToWavelengthPolynomial);

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