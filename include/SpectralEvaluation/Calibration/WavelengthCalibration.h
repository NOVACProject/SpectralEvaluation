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

enum class InstrumentLineshapeEstimationOption
{
    None = 0,
    Gaussian,
};

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

    /// <summary>
    /// The estimation of the instrument line shape.
    /// This is only set if the instrument line shape is estimated in the 
    /// </summary>
    novac::CCrossSectionData estimatedInstrumentLineShape;
};

struct WavelengthCalibrationSettings
{
    /// <summary>
    /// The intial estimate for the pixel to wavelength mapping.
    /// </summary>
    std::vector<double> initialPixelToWavelengthMapping;

    /// <summary>
    /// The initial estimate for the instrument line shape (measured or estimated).
    /// </summary>
    novac::CCrossSectionData initialInstrumentLineShape;

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

    /// <summary>
    /// The option for how, and if, the instrument line shape should also be estimated 
    /// during the wavelength calibration procedure.
    /// </summary>
    InstrumentLineshapeEstimationOption estimateInstrumentLineShape = InstrumentLineshapeEstimationOption::None;
};

/// <summary>
/// Reads the pixel-to-wavelength mapping from a .clb calibration 
///     data file which is expected to contain a single column of data.
/// Notice, if the file contains two columns then this will return the second!
/// </summary>
std::vector<double> GetPixelToWavelengthMappingFromFile(const std::string& clbFile);

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
 * TODO: Figure out if this is actually used anywhere. If not then this could probably be removed.
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
    ///   high resolution fraunhofer spectrum assuming that the provided instrument line shape
    ///   is the correct line shape for the instrument.
    /// </summary>
    SpectrometerCalibrationResult DoWavelengthCalibration(const CSpectrum& measuredSpectrum);

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
