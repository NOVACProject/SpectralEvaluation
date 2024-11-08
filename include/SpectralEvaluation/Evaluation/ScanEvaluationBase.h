#pragma once

#include <string>
#include <SpectralEvaluation/Log.h>

namespace Configuration
{
struct CDarkSettings;
struct CSkySettings;
}

namespace novac
{
struct SpectrometerModel;
class CScanFileHandler;
class CSpectrum;
class CEvaluationBase;
class CFitWindow;

/** Reads a spectrum from .std or .txt files. @return true if successful */
bool ReadSpectrumFromFile(const std::string& fullFilename, CSpectrum& spec);

/** ScanEvaluationBase is the base class for the ScanEvaluation-classes found in
    NovacPPP and NovacProgram. This collects the common elements between the two program */
class ScanEvaluationBase
{
public:
    ScanEvaluationBase(novac::ILogger& log);
    virtual ~ScanEvaluationBase();

    /** Setting the option for wheather the spectra are averaged or not.
        averaged=true means that the spectra are averaged instead of summed.
        The default Novac option is that the spectra are summed (averaged=false). */
    void SetOption_AveragedSpectra(bool averaged) { this->m_averagedSpectra = averaged; }

    /** Attempts to determine the optimium shift and squeeze value to use when evaluating the provided scan.
            This assumes that the provided fitWindow has a defined fitWindow.fraunhoferRef member with an already read in cross section.
            This will set the shift/squeeze of the fraunhofer-ref to 'free' and link all other references to it.
        @param fitWindow Defines the references and fit-range to use.
        @param scan The scan to evaluate.
        @return a new evaluator (with a new fit-window set) to use to evaluate this scan with the shift/squeeze fixed to a good value.
        @return nullptr if the determination failed. */
    CEvaluationBase* FindOptimumShiftAndSqueezeFromFraunhoferReference(
        novac::LogContext context,
        const CFitWindow& fitWindow,
        const novac::SpectrometerModel& spectrometerModel,
        const Configuration::CDarkSettings& darkSettings,
        CScanFileHandler& scan);

    /** @return the index of the spectrum with the 'most suitable intensity for fitting', i.e. this is the
        spectrum with the highest intensity which isn't (close to being) saturated. */
    static int GetIndexOfSpectrumWithBestIntensity(
        novac::LogContext context,
        const CFitWindow& fitWindow,
        const novac::SpectrometerModel& spectrometerModel,
        CScanFileHandler& scan);

protected:

    novac::ILogger& m_log;

    /** This is the lower edge of the fit region used in the last evaluation performed (in pixels).
        Used for reference and further processing. */
    int m_fitLow = 320;

    /** This is the upper edge of the fit region used in the last evaluation performed (in pixels).
        Used for reference and further processing. */
    int m_fitHigh = 460;

    /** True if the spectra are averaged, not summed as they normally are. Default value is false. */
    bool m_averagedSpectra = false;

    /** The last error message, formatted as a string, set by any function which can return true/false to indicate success or failure. */
    std::string m_lastErrorMessage = "";

    /** Retrieves the dark spectrum to use for the provided spectrum given the settings from the user.
        @return true on success.
        This will set m_lastErrorMessage if failed to get the spectrum. */
    bool GetDark(CScanFileHandler& scan, const CSpectrum& spec, CSpectrum& dark, const Configuration::CDarkSettings* darkSettings);

    /** This returns the sky spectrum that is to be used in the fitting.
        Which spectrum to be used is taken from the given settings.
        @return true on success.
        This will set m_lastErrorMessage if failed to get the spectrum. */
    bool GetSky(CScanFileHandler& scan, const Configuration::CSkySettings& settings, CSpectrum& sky);

    /** Creates a spectrum which is the average of all non-saturated and not-dark spectra in the given scan.
        @return true on successful spectrum creation. */
    bool GetSkySpectrumFromAverageOfGoodSpectra(CScanFileHandler& scan, CSpectrum& sky) const;

    /** Reads a 'sky' spectrum from file. This can read the 'sky' spectrum in one .pak
        file or one spectrum from a .std file.
        @return true on successful spectrum read.
        Sets m_lastErrorMessage if the reading failed */
    bool GetSkySpectrumFromFile(const std::string& filename, CSpectrum& sky);
};
}

