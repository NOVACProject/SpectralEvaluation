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
    /// This is only set if the instrument line shape is set to be estimated 
    //  in the process, otherwise this is empty
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

class WavelengthCalibrationFailureException : public std::exception
{
private:
    const char* const m_msg = "";

public:
    WavelengthCalibrationFailureException(const char* msg) :
        m_msg(msg)
    {}

    const char* what() const noexcept override final { return m_msg; }
};

/// <summary>
/// Reads the pixel-to-wavelength mapping from a .clb calibration 
///     data file which is expected to contain a single column of data.
/// Notice, if the file contains two columns then this will return the second!
/// </summary>
std::vector<double> GetPixelToWavelengthMappingFromFile(const std::string& clbFile);

/// <summary>
/// This is a helper structure used to extract the internal state of the pixel-to-wavelength calibration
/// of a spectrometer from a measured mercury spectrum.
/// </summary>
struct MercurySpectrumCalibrationState
{
    /// <summary>
    /// This lists all the peaks found in the mercury spectrum (defined in pixels).
    /// </summary>
    std::vector<SpectrumDataPoint> peaks;

    /// <summary>
    /// This is a list of peaks found in the spectrum but for which we didn't find a
    /// suitable emission line in the theoretical spectrum, or the line seems to be not fully resolved.
    /// </summary>
    std::vector<SpectrumDataPoint> rejectedPeaks;

    /// <summary>
    /// If the Mercury calibration failed, then this will be filled in with the reason why
    /// </summary>
    std::string errorMessage;
};

/**
 * @brief Performs a wavelength calibration of a spectrometer using a measured mercury spectrum.
 *  This will identify the peaks in the measured spectrum and fit a polynomial to them such that the
 *  pixel-to-wavelength mapping for the entire spectrum can be calculated.
 * @param measuredMercurySpectrum The measured spectrum
 * @param polynomialOrder The order of the polynomial to fit (in the range 1-3).
 * @param initialPixelToWavelength The initial guess for the pixel to wavelength calibration of the spectrum,
 *  This will be used to constrain the identification of the mercury lines.
 * @param result Will on successful return be filled with the resulting pixel-to-wavelength calibration.
 * @param state If not null, then this will be filled with information on how the calibration did perform.
 * @return True if the calibration was successful. */
bool MercuryCalibration(
    const CSpectrum& measuredMercurySpectrum,
    int polynomialOrder,
    const std::vector<double>& initialPixelToWavelength,
    SpectrometerCalibrationResult& result,
    MercurySpectrumCalibrationState* state = nullptr);

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
    /// @throws std::invalid_argument if any of the incoming parameters is invalid.
    /// @throws WavelengthCalibrationFailureException if the calibration fails.
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
