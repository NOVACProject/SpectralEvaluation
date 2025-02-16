#pragma once

#include <memory>
#include <string>
#include <vector>
#include <SpectralEvaluation/Log.h>
#include <SpectralEvaluation/Spectra/WavelengthRange.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Spectra/SpectrometerModel.h>

namespace novac
{
    class InstrumentCalibration;
}

/**
 * @brief  Abstract base class for the
    MobileDoasWavelengthCalibrationController and
    NovacProgramWavelengthCalibrationController
*/
class WavelengthCalibrationController
{
public:
    WavelengthCalibrationController(novac::ILogger& log);
    virtual ~WavelengthCalibrationController();

    enum class InstrumentLineShapeFitOption
    {
        None = 0,
        SuperGaussian = 1
    };

    /** The full path to the high resolved solar spectrum */
    std::string m_solarSpectrumFile;

    /** The full path to a file which contains an initial calibration
         This can be either a full calibration file, as saved from another novac program,
         or just the pixel-to-wavelength mapping file. */
    std::string m_initialCalibrationFile;

    /** The full path to a file which contains an initial measured line shape.
         This may be left out if m_initialCalibrationFile does contain an instrument line shape. */
    std::string m_initialLineShapeFile;

    /** The option for if an instrument line shape should be fitted as well during
         the retrieval of the pixel-to-wavelength calibration. */
    InstrumentLineShapeFitOption m_instrumentLineShapeFitOption;

    /** The wavelength region in which the instrument line shape should be fitted (in nm). */
    novac::WavelengthRange m_instrumentLineShapeFitRegion;

    /** One additional cross section to be included when performing the instrument line shape fitting.
        Typically this is the high resolution cross section for ozone. */
    std::string m_crossSectionsForInstrumentLineShapeFitting;

    /** This is the initial calibration, used as a starting point for the calibration routine.
        This will be set from the contents of m_initialCalibrationFile and m_initialLineShapeFile
        when RunCalibration() is called. */
    std::unique_ptr<novac::InstrumentCalibration> m_initialCalibration;

    /** This is the result of the wavelength calibration.
         Can only be set after calling RunCalibration().
         Notice that this only contains the portions of the calibration which have been determined by RunCalibration(),
         i.e. if no instrument line shape was fitted then this will not contain an instrument line shape description.
         To get a full calibration in that case, call GetFinalCalibration() which will combine m_resultingCalibration with m_initialCalibration
            in order to get a full calibration. */
    std::unique_ptr<novac::InstrumentCalibration> m_resultingCalibration;

    // The model of the spectrometer which generated the spectra.
    // This can be nullptr in which case the model is determined from the measurement.
    std::unique_ptr<novac::SpectrometerModel> m_spectrometerModel;

    /** User friendly description of the fitted parameters for the instrument line shape function. */
    std::vector<std::pair<std::string, std::string>> m_instrumentLineShapeParameterDescriptions;

    /** If the calibration fails, for some reason, then this message should be set to indicate why. */
    std::string m_errorMessage;

    /** An elementary log, will contain debugging information from running the calibration. */
    std::vector<std::string> m_logMessages;

    /** An optional override of the check for the maximum intensity of the spectrometer model.
        Sometimes the intensity is specified by the user directly and not taken from the device (e.g. directory reading mode in MobileDoas)
        and in this case we do not want to use the information in the spectrum about the model for calculating a saturation ratio.
        This value will be used if set to a value > 0.0 . */
    double m_spectrometerMaximumIntensityForSingleReadout = -1.0;

    /** Set to false if the spectrometer (or the software running the spectrometer) averages the spectra
        instead of adding them.
        Typically, MobileDoas averages and NovacProgram adds. */
    bool m_spectraAreAverages = false;

    struct WavelengthCalibrationDebugState
    {
        WavelengthCalibrationDebugState(size_t estimatedSize)
        {
            inlierCorrespondencePixels.reserve(estimatedSize);
            inlierCorrespondenceWavelengths.reserve(estimatedSize);
            outlierCorrespondencePixels.reserve(estimatedSize);
            outlierCorrespondenceWavelengths.reserve(estimatedSize);
        }

        std::vector<double> initialPixelToWavelengthMapping;

        std::vector<double> inlierCorrespondencePixels;
        std::vector<double> inlierCorrespondenceWavelengths;
        std::vector<double> inlierCorrespondenceMeasuredIntensity; //< the intensity of the measured spectrum
        std::vector<double> inlierCorrespondenceFraunhoferIntensity; //< the intensity of the Fraunhofer spectrum

        std::vector<double> outlierCorrespondencePixels;
        std::vector<double> outlierCorrespondenceWavelengths;

        std::vector<double> measuredSpectrum;
        novac::CSpectrumInfo spectrumInfo;

        // All the keypoints from the measured spectrum
        std::vector<double> measuredSpectrumKeypointPixels;
        std::vector<double> measuredSpectrumKeypointIntensities;

        // The inlier keypoints from the measured spectrum
        std::vector<double> measuredSpectrumInlierKeypointPixels;
        std::vector<double> measuredSpectrumInlierKeypointIntensities;

        std::vector<double> fraunhoferSpectrum;

        // All the keypoints from the Fraunhofer spectrum
        std::vector<double> fraunhoferSpectrumKeypointWavelength;
        std::vector<double> fraunhoferSpectrumKeypointIntensities;

        // The inlier keypoints from the Fraunhofer spectrum
        std::vector<double> fraunhoferSpectrumInlierKeypointWavelength;
        std::vector<double> fraunhoferSpectrumInlierKeypointIntensities;
    };

    WavelengthCalibrationDebugState m_calibrationDebug;

    /** Performs the actual wavelength calibration.
        @throws std::invalid_argument with an explanatory error message if the calibration fails. */
    void RunCalibration();

    /** Returns a full instrument calibration by combining the portions determined by RunCalibration() with m_initialCalibration */
    std::unique_ptr<novac::InstrumentCalibration> GetFinalCalibration() const;

    /** Saves the resulting pixel-to-wavelength and instrument-line-shape information as a .std file. */
    void SaveResultAsStd(const std::string& filename);

    /** Saves the resulting pixel-to-wavelength mapping information as a .clb file. */
    void SaveResultAsClb(const std::string& filename);

    /** Saves the resulting instrument line shape information as a .slf file. */
    void SaveResultAsSlf(const std::string& filename);

    /** Empties the previous result (if any). */
    void ClearResult();

protected:

    novac::ILogger& m_log;

    virtual void ReadInput(novac::CSpectrum& measurement) = 0;

    /** Guesses for an instrment line shape from the measured spectrum. Useful if no measured instrument line shape exists */
    void CreateGuessForInstrumentLineShape(const std::string& solarSpectrumFile, const novac::CSpectrum& measuredSpectrum, novac::InstrumentCalibration& calibration);

    /** Checks the provided (dark corrected) spectrum and makes sure that it is good enough for the calibration to succeed.
        @throws an std::invalid_argument exception if the spectrum isn't good enough.  */
    void CheckSpectrumQuality(const novac::CSpectrum& spectrum, const novac::SpectrometerModel& spectrometerModel) const;

    double GetSpectrometerMaxIntensityForSingleReadout(const novac::CSpectrum& spectrum, const novac::SpectrometerModel& spectrometerModel) const;

    // Retrieves a spectrometer model to be used. Either from the 'm_spectrometerModel' if that is set, or by guessing the model from the serial.
    novac::SpectrometerModel GetModelForMeasurement(const std::string& deviceSerial) const;

    void Log(const std::string& message);
    void Log(const std::string& message, double value);
    void Log(const std::string& message, const std::string& messagePart2);

};

/** Subclass of WavelengthCalibrationController which takes the path to a measured spectrum
    and a measured dark spectrum as input */
class MobileDoasWavelengthCalibrationController : public WavelengthCalibrationController
{
public:
    MobileDoasWavelengthCalibrationController(novac::ILogger& log)
        : WavelengthCalibrationController(log)
    {
        // MobileDOAS will always average spectra together
        m_spectraAreAverages = true;
    }

    /** The full path to the spectrum to calibrate */
    std::string m_inputSpectrumFile;

    /** The full path to the dark spectrum of the spectrum to calibrate */
    std::string m_darkSpectrumFile;

protected:

    virtual void ReadInput(novac::CSpectrum& measurement) override;
};

/** Subclass of WavelengthCalibrationController which doesn't read the data from disk but
    uses a spectrum in memory as target for the calibration. */
class InMemoryWavelengthCalibrationController : public WavelengthCalibrationController
{
public:
    InMemoryWavelengthCalibrationController(novac::ILogger& log)
        : WavelengthCalibrationController(log) { }

    /** The spectrum to calibrate. Notice that this may be modified during the calibration. */
    novac::CSpectrum m_measuredSpectrum;

protected:

    virtual void ReadInput(novac::CSpectrum& measurement) override;
};

