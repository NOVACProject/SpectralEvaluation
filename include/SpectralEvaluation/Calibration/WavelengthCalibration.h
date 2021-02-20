#pragma once

#include <vector>
#include <memory>
#include <string>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>

// ---------------------------------------------------------------------------------------------------------------
// ----------- This header contains methods used to perform wavelength calibration of a spectrometer -------------
// ---------------------------------------------------------------------------------------------------------------

class CSpectrum;
struct SpectrumDataPoint;

namespace novac
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

/// <summary>
/// Reads the pixel-to-wavelength mapping from a .clb calibration 
///     data file which is expected to contain a single column of data.
/// Notice, if the file contains two columns then this will return the second!
/// </summary>
std::vector<double> GetPixelToWavelengthMappingFromFile(const std::string& clbFile);

/// <summary>
/// This is a helper class for generating a Fraunhofer spectrum from a high resolved
/// solar spectrum, a likewise high resolved ozone spectrum and a given instrument setup.
/// </summary>
class FraunhoferSpectrumGeneration
{
public:
    /// <summary>
    /// Sets up the generation parameters
    /// </summary>
    /// <param name="highResolutionSolarAtlas">The full path to the high resolved solar atlas. This must be in nm air.</param>
    /// <param name="highResolutionOzoneCrossSection">The full path to the high resolved ozone cross section. 
    ///     This must have x-axis unit of nm air and y-axis unit of molecules / cm2</param>
    FraunhoferSpectrumGeneration(const std::string& highResolutionSolarAtlas, const std::string& highResolutionOzoneCrossSection)
        : solarAtlasFile(highResolutionSolarAtlas), ozoneCrossSectionFile(highResolutionOzoneCrossSection)
    {
    }

    /// <summary>
    /// Creates a Fraunhofer reference spectrum using the provided pixel-to-wavelength mapping and measured instrument line shape.
    /// </summary>
    /// <param name="pixelToWavelengthMapping">The wavelength (in nm air) for each pixel on the detector.</param>
    /// <param name="measuredSlf">A measurement of the instrument line shape</param>
    /// <param name="ozoneColumn">The total column of ozone in the resulting Fraunhofer spectrum. In molecules / cm2</param>
    /// <returns></returns>
    std::unique_ptr<CSpectrum> GetFraunhoferSpectrum(
        const std::vector<double>& pixelToWavelengthMapping,
        const ::Evaluation::CCrossSectionData& measuredSlf,
        double ozoneColumn);

private:
    const std::string solarAtlasFile;
    const std::string ozoneCrossSectionFile;
};



/// <summary>
/// Generates a Fraunhofer reference spectrum from the provided data.
/// </summary>
/// <param name="pixelToWavelengthMapping"></param>
/// <param name="initialSlfFile"></param>
/// <returns></returns>

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
