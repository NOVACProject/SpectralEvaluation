#pragma once

#include <string>

class CSpectrum;
namespace FileHandler
{
    class CScanFileHandler;
}
namespace Configuration
{
    struct CDarkSettings;
}

namespace Evaluation
{
    /** ScanEvaluationBase is the base class for the ScanEvaluation-classes found in 
        NovacPPP and NovacProgram. This collects the common elements between the two program */
    class ScanEvaluationBase
    {
    public:
        ScanEvaluationBase();
        virtual ~ScanEvaluationBase();

        /** Setting the option for wheather the spectra are averaged or not.
            averaged=true means that the spectra are averaged instead of summed.
            The default Novac option is that the spectra are summed (averaged=false). */
        void SetOption_AveragedSpectra(bool averaged) { this->m_averagedSpectra = averaged; }

    protected:
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

        /** Models a dark spectrum for the provided measured spectrum from an offset and a dark-current spectrum.
            @param spec The spectrum which needs dark-correcting. This does not modify that spectrum...
            @param dark Will on successful return be filled with the dark spectrum.
            @return true if successful or false if the spectrum could not be modelled. This sets m_lastErrorMessage. */
        bool ModelDarkSpectrum(FileHandler::CScanFileHandler &scan, const Configuration::CDarkSettings& darkSettings, const CSpectrum &spec, CSpectrum &dark);

        /** Retrieves an offset-spectrum (spectrum with very short exposure time) from the given scan using the options.
            If the options says to always use the measured then the spectrum is taken from the scan. 
            If the option says to use an user-supplied offset spectrum then the offset-spectrum is read from disk. */
        bool GetOffsetSpectrum(FileHandler::CScanFileHandler &scan, const Configuration::CDarkSettings& darkSettings, CSpectrum &offsetSpectrum);

        /** Retrieves a dark-current-spectrum (spectrum with very long exposure time) from the given scan using the options.
            @param needsOffsetCorrection Will be set to true if the provided spectrum does need to be offset-corrected.
            If the options says to always use the measured then the spectrum is taken from the scan.
            If the option says to use an user-supplied offset spectrum then the offset-spectrum is read from disk. */
        bool GetDarkCurrentSpectrum(FileHandler::CScanFileHandler &scan, const Configuration::CDarkSettings& darkSettings, CSpectrum &darkCurrent, bool& needsOffsetCorrection);

        /** Reads a spectrum from .std or .txt files. @return true if successful */
        bool ReadSpectrumFromFile(const std::string& fullFilename, CSpectrum& spec);

        /** Retrieves the dark spectrum to use for the provided spectrum given the settings from the user.
            @return true on success. */
        bool GetDark(FileHandler::CScanFileHandler& scan, const CSpectrum &spec, CSpectrum &dark, const Configuration::CDarkSettings* darkSettings);

        /** Creates a spectrum which is the average of all non-saturated and not-dark spectra in the given scan.
            @return true on successful spectrum creation. */
        bool GetSkySpectrumFromAverageOfGoodSpectra(FileHandler::CScanFileHandler& scan, CSpectrum &sky) const;

        /** Reads a 'sky' spectrum from file. This can read the 'sky' spectrum in one .pak
            file or one spectrum from a .std file.
            @return true on successful spectrum read.
            Sets m_lastErrorMessage if the reading failed */
        bool GetSkySpectrumFromFile(const std::string& filename, CSpectrum& sky);
    };
}