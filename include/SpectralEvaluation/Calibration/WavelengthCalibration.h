#pragma once

#include <utility>
#include <memory>
#include <vector>
#include <string>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Spectra/SpectrumUtils.h>

// ---------------------------------------------------------------------------------------------------------------
// ----------- This header contains methods used to perform wavelength calibration of a spectrometer -------------
// ---------------------------------------------------------------------------------------------------------------

namespace novac
{

class CSpectrum;
struct SpectrumDataPoint;
struct Correspondence;

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
/// Notice that this class will read in the high-resolved solar spectrum when needed (calling GetFraunhoferSpectrum)
/// and will keep it in memory to save loading time. If memory is a consern, then make sure that this object gets destructed when no longer needed.
/// </summary>
class FraunhoferSpectrumGeneration
{
public:
    /// <summary>
    /// Sets up the generation parameters
    /// </summary>
    /// <param name="highResolutionSolarAtlas">The full path to the high resolved solar atlas. This must be in nm air.</param>
    /// <param name="highResolutionCrossSections">The full path to a set of high resolved molecular cross section together with the total column for them.
    ///     These must have x-axis unit of nm air and y-axis unit of molecules / cm2</param>
    FraunhoferSpectrumGeneration(const std::string& highResolutionSolarAtlas, const std::vector<std::pair<std::string, double>>& highResolutionCrossSections)
        : solarAtlasFile(highResolutionSolarAtlas), crossSectionsToInclude(highResolutionCrossSections.size())
    {
        for (size_t ii = 0; ii < highResolutionCrossSections.size(); ++ii)
        {
            crossSectionsToInclude[ii].path = highResolutionCrossSections[ii].first;
            crossSectionsToInclude[ii].totalColumn = highResolutionCrossSections[ii].second;
        }
    }

    /// <summary>
    /// Creates a Fraunhofer reference spectrum using the provided pixel-to-wavelength mapping and measured instrument line shape.
    /// </summary>
    /// <param name="pixelToWavelengthMapping">The wavelength (in nm air) for each pixel on the detector.</param>
    /// <param name="measuredInstrumentLineShape">A measurement of the instrument line shape</param>
    /// <returns>The high resolution solar spectrum convolved with the measured slf and resample to the provided grid.</returns>
    std::unique_ptr<CSpectrum> GetFraunhoferSpectrum(
        const std::vector<double>& pixelToWavelengthMapping,
        const novac::CCrossSectionData& measuredInstrumentLineShape);

    /// <summary>
    /// Creates a Fraunhofer reference spectrum using the provided pixel-to-wavelength mapping and measured instrument line shape.
    ///     A DOAS fit is applied in order to determine the total columns of each high resolution cross section in this setup.
    ///     I.e., the total columns passed in to the constructor are ignored.
    /// This will throw an std::invalid_argument if no high resolution cross section was provided at creation (indicating an invalid setup)
    /// </summary>
    /// <param name="measuredInstrumentLineShape">A measurement of the instrument line shape</param>
    /// <param name="measuredSpectrum">A measured spectrum. </param>
    /// <returns>The high resolution solar spectrum convolved with the measured slf and resample to the provided grid.</returns>
    std::unique_ptr<CSpectrum> GetFraunhoferSpectrumMatching(
        const std::vector<double>& pixelToWavelengthMapping,
        const novac::CSpectrum& measuredSpectrum,
        const novac::CCrossSectionData& measuredInstrumentLineShape);

private:
    /// <summary>
    /// The path and filename of the solar atlas file to use.
    /// </summary>
    const std::string solarAtlasFile;

    struct AbsorbingCrossSection
    {
        AbsorbingCrossSection() = default;

        AbsorbingCrossSection(const std::pair<std::string, double>& value)
            : path(value.first), totalColumn(value.second)
        {
        }

        std::string path;
        double totalColumn = 0.0;
        std::unique_ptr<novac::CCrossSectionData> crossSectionData;
    };

    /// <summary>
    /// The path and total column of the high resolved absorption cross section files to include.
    /// </summary>
    std::vector<AbsorbingCrossSection> crossSectionsToInclude;

    /// <summary>
    /// The read in high resolution solar cross section, saved in order to reduce file-io time.
    /// </summary>
    std::unique_ptr<novac::CCrossSectionData> solarCrossSection;

    /// <summary>
    /// This creates a Fraunhofer spectrum using the provided Absorbing cross sections instead of using the member
    /// </summary>
    /// <param name="pixelToWavelengthMapping"></param>
    /// <param name="measuredInstrumentLineShape"></param>
    /// <param name="crossSectionsToInclude"></param>
    /// <returns></returns>
    std::unique_ptr<CSpectrum> GetFraunhoferSpectrum(
        const std::vector<double>& pixelToWavelengthMapping,
        const novac::CCrossSectionData& measuredInstrumentLineShape,
        std::vector<AbsorbingCrossSection>& crossSectionsToInclude);

    void ReadSolarCrossSection();
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
    /// </summary>
    SpectrometerCalibrationResult DoWavelengthCalibration(const CSpectrum& measuredSpectrum, const novac::CCrossSectionData& measuredInstrumentLineShape);

    /// <summary>
    /// Simple structure used to save the internal state of the wavelength calibration. For inspection and debugging
    /// </summary>
    struct SpectrumeterCalibrationState
    {
        std::unique_ptr<CSpectrum> measuredSpectrum;
        std::unique_ptr<CSpectrum> fraunhoferSpectrum;
        std::unique_ptr<CSpectrum> originalFraunhoferSpectrum;
        std::vector<novac::SpectrumDataPoint> measuredKeypoints;
        std::vector<novac::SpectrumDataPoint> fraunhoferKeypoints;
        std::vector<double> measuredSpectrumEnvelopePixels;
        std::vector<double> measuredSpectrumEnvelopeIntensities;
        std::vector<novac::Correspondence> allCorrespondences;
        std::vector<bool> correspondenceIsInlier;
    };

    /// <summary>
    /// Helper function which retrieves the state and result of the last call to 'DoWavelengthCalibration', for debugging.
    /// </summary>
    const WavelengthCalibrationSetup::SpectrumeterCalibrationState& GetLastCalibrationSetup() { return this->calibrationState; }

private:

    WavelengthCalibrationSettings settings;

    SpectrumeterCalibrationState calibrationState;

    static std::vector<double> GetPixelToWavelengthMapping(const std::vector<double>& polynomialCoefficients, size_t detectorSize);

};

}
