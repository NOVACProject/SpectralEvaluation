#pragma once

#include <utility>
#include <memory>
#include <vector>
#include <string>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>

// ---------------------------------------------------------------------------------------------------------------
// ----------- This header contains methods used to perform wavelength calibration of a spectrometer -------------
// ---------------------------------------------------------------------------------------------------------------

class CSpectrum;
struct SpectrumDataPoint;

namespace novac
{

/// <summary>
/// This is the result of running the wavelength calibration routine.
/// </summary>
struct SpectrometerCalibrationResult
{
    /// <summary>
    /// The final estimate for the pixel to wavelength mapping.
    /// </summary>
    std::vector<double> pixelToWavelengthMapping;

    /// <summary>
    /// The coefficients of the pixel-to-wavelength mapping polynomial
    /// </summary>
    std::vector<double> pixelToWavelengthMappingCoefficients;

};

struct WavelengthCalibrationSettings
{
    /// <summary>
    /// The intial estimate for the pixel to wavelength mapping.
    /// </summary>
    std::vector<double> initialPixelToWavelengthMapping;

    /// <summary>
    /// The path to the high resolved solar atlas to be used.
    ///  Assumed that this has the wavelength unit of nm air.
    /// </summary>
    std::string highResSolarAtlas;

    /// <summary>
    /// A set of high resolved cross sections to include into the generation
    /// of the fraunhofer spectrum. Each pair makes up a path to the cross section
    /// file on disk and a total column of the cross section to use.
    /// Only cross sections with a total column != 0 will be included.
    /// </summary>
    std::vector<std::pair<std::string, double>> crossSections;
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
    FraunhoferSpectrumGeneration(const std::string& highResolutionSolarAtlas)
        : solarAtlasFile(highResolutionSolarAtlas)
    {
    }

    /// <summary>
    /// Creates a Fraunhofer reference spectrum using the provided pixel-to-wavelength mapping and measured instrument line shape.
    /// </summary>
    /// <param name="pixelToWavelengthMapping">The wavelength (in nm air) for each pixel on the detector.</param>
    /// <param name="measuredInstrumentLineShape">A measurement of the instrument line shape</param>
    /// <param name="highResolutionOzoneCrossSection">The full path to the high resolved ozone cross section. 
    ///     This must have x-axis unit of nm air and y-axis unit of molecules / cm2</param>
    /// <param name="ozoneColumn">The total column of ozone in the resulting Fraunhofer spectrum. In molecules / cm2</param>
    /// <returns>The high resolution solar spectrum convolved with the measured slf and resample to the provided grid.</returns>
    std::unique_ptr<CSpectrum> GetFraunhoferSpectrum(
        const std::vector<double>& pixelToWavelengthMapping,
        const ::Evaluation::CCrossSectionData& measuredInstrumentLineShape,
        const std::string& highResolutionOzoneCrossSection,
        double ozoneColumn);

    /// <summary>
    /// Creates a Fraunhofer reference spectrum using the provided pixel-to-wavelength mapping and measured instrument line shape.
    /// </summary>
    /// <param name="pixelToWavelengthMapping">The wavelength (in nm air) for each pixel on the detector.</param>
    /// <param name="measuredInstrumentLineShape">A measurement of the instrument line shape</param>
    /// <param name="highResolutionCrossSections">The full path to a set of high resolved molecular cross section together with the total column for them.
    ///     These must have x-axis unit of nm air and y-axis unit of molecules / cm2</param>
    /// <returns>The high resolution solar spectrum convolved with the measured slf and resample to the provided grid.</returns>
    std::unique_ptr<CSpectrum> GetFraunhoferSpectrum(
        const std::vector<double>& pixelToWavelengthMapping,
        const ::Evaluation::CCrossSectionData& measuredInstrumentLineShape,
        const std::vector<std::pair<std::string, double>>& highResolutionCrossSections);

private:
    const std::string solarAtlasFile;
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



/// <summary>
/// WavelengthCalibrationSetup is the setup of a calibration run
///     and contains all necessary elements to perform the calibration.
/// </summary>
class WavelengthCalibrationSetup
{
public:
    WavelengthCalibrationSetup(const WavelengthCalibrationSettings& calibrationSettings);

    /// <summary>
    /// This performs the actual calibration of a measured spectrum against a 
    ///   high resolution fraunhofer spectrum assuming that the provided measured instrument line shape
    ///   is the correct line shape for the instrument.
    /// This will modify the measured spectrum, by normalizing its intensity to the range [0, 1]
    /// </summary>
    SpectrometerCalibrationResult DoWavelengthCalibration(CSpectrum& measuredSpectrum, const Evaluation::CCrossSectionData& measuredInstrumentLineShape);

    // TODO: Create a way to get out some more internal information regarding the calibration, for debugging

private:

    WavelengthCalibrationSettings settings;

    static std::vector<double> GetPixelToWavelengthMapping(const std::vector<double>& polynomialCoefficients, size_t detectorSize);

};

}
