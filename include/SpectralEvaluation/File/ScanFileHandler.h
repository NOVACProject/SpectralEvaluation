#pragma once

#include <memory>
#include <vector>
#include <SpectralEvaluation/DateTime.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Spectra/IScanSpectrumSource.h>

namespace novac
{

/** ScanFileHandler is a class to read in information from the scan-files
    (all the spectra from one scan are supposed to be packed together in one file in Manne's 'pak'-format.
    Each instance of 'CScanFileHandler' is capable of reading data from one .pak-file.
    Each instance of this class should be initialized by first calling 'CheckScanFile' which will read in
    the data on the file and then create*/
class CScanFileHandler : public IScanSpectrumSource
{
public:
    CScanFileHandler(ILogger& log);

    // ----------------------------------------------------------------------
    // ---------------------- PUBLIC DATA -----------------------------------
    // ----------------------------------------------------------------------

    /** The serial number of the spectrometer which collected the spectra in this scan.
        Set by CheckScanFile */
    std::string m_device = "";

    /** The channel of the spectrometer which was used for collecting this scan.
         (if a SD2000 with multiple channels is used, one spectrometer should be
        configured for each channel).
        Set by CheckScanFile */
    unsigned char m_channel = 0;

    /** The time (UTC) when the measurement started.
        Set by CheckScanFile */
    CDateTime m_startTime;

    /** The time (UTC) when the measurement was finished */
    CDateTime m_stopTime;

    /** If any error occurs in the reading of the file, this int is set to
        any of the errors defined int 'SpectrumIO.h. */
    int m_lastError = 0;

    /** If any error occurs in the reading of the file, this is set to a
    *   human readable error message describing the problem. */
    std::string m_lastErrorMessage = "";

    // ----------------------------------------------------------------------
    // --------------------- PUBLIC METHODS ---------------------------------
    // ----------------------------------------------------------------------

    /** Checks the scan saved in the given filename
        If any file-error occurs the parameter 'm_lastError' will be set.
        @param fileName - the name of the file in which the spectra of the scan are saved.
        @return true on success.
        @return false if any error occurs */
    bool CheckScanFile(novac::LogContext context, const std::string& fileName);

    /** Gets the next spectrum in the scan.
        If any file-error occurs the parameter 'm_lastError' will be set.
        @param spec - will on successful return be filled with the newly read spectrum.
        @return the number of spectra read (1 if successful, otherwise 0).*/
    int GetNextSpectrum(novac::LogContext context, CSpectrum& spec);

    virtual int GetNextMeasuredSpectrum(novac::LogContext context, CSpectrum& spec) override
    {
        return 1 - GetNextSpectrum(context, spec);
    }

    /** Returns the desired spectrum in the scan.
        If any file-error occurs the parameter 'm_lastError' will be set.
        @param spec - will on successful return be filled with the newly read spectrum.
        @param specNo - The zero-based index into the scan-file.
        @return the number of spectra read (1 if successful, otherwise 0) */
    int GetSpectrum(novac::LogContext context, CSpectrum& spec, long specNo);

    virtual int GetSpectrum(novac::LogContext context, int spectrumNumber, CSpectrum& spec) override
    {
        return 1 - GetSpectrum(context, spec, (long)spectrumNumber);
    }

    /** Gets the dark spectrum of the scan
        @return 0 if there is a dark spectrum, else non-zero */
    virtual int GetDark(CSpectrum& spec) const override;

    /** Gets the sky spectrum of the scan
        @return 0 if there is a sky spectrum, else non-zero */
    virtual int GetSky(CSpectrum& spec) const override;

    /** Gets the offset spectrum of the scan - if any
        @return 0 if there is an offset spectrum, else non-zero */
    virtual int GetOffset(CSpectrum& spec) const override;

    /** Gets the dark-current spectrum of the scan - if any
        @return 0 if there is a dark-current spectrum, else non-zero */
    virtual int GetDarkCurrent(CSpectrum& spec) const override;

    /** Returns the interlace steps for the spectra in this scan-file.
             @return the interlace steps for the spectra in this scan.
             @return -1 if the function 'CheckScanFile' has not been called */
    int	GetInterlaceSteps() const;

    /** Returns the length of the spectra in this scan-file.
             @return the spectrum-length for the spectra in this scan.
             @return -1 if the function 'CheckScanFile' has not been called */
    int	GetSpectrumLength() const;

    /** Returns the start-channel for the spectra in this scan-file.
            This is the pixel on the detector for which corresponds to the first
                datapoint in the spectra (normally 0).
             @return the start-channel for the spectra in this scan.
             @return -1 if the function 'CheckScanFile' has not been called */
    int GetStartChannel() const;

    /** Retrieves GPS-information from the spectrum files */
    const CGPSData GetGPS() const;

    /** Retrieves compass-information from the spectrum files */
    double GetCompass() const;

    /** Retrieves the name of the file that this object is working on */
    virtual std::string GetFileName() const override { return m_fileName; }

    /** Resets the m_specReadSoFarNum to start reading from the first spectrum again */
    void ResetCounter();

    /** Retrieves the total number of spectra in the .pak-file (including sky and dark) */
    virtual int GetSpectrumNumInFile() const override;

    /** Retrieves the time (UTC) when the first spectrum was collected */
    virtual CDateTime GetScanStartTime() const override { return m_startTime; }

    /** Retrieves the time (UTC) when the last spectrum was collected */
    virtual CDateTime GetScanStopTime() const override { return m_stopTime; }

    /** Retrieves the serial number of the device which collected this scan. */
    virtual std::string GetDeviceSerial() const override { return m_device; }

private:
    // ----------------------------------------------------------------------
    // ---------------------- PRIVATE DATA ----------------------------------
    // ----------------------------------------------------------------------

    ILogger& m_log;

    /** The dark spectrum. this may be null if there's no dark spectrum in the scan. */
    std::unique_ptr<CSpectrum> m_dark;

    /** The sky spectrum */
    std::unique_ptr<CSpectrum> m_sky;

    /** The offset spectrum - if any */
    std::unique_ptr<CSpectrum> m_offset;

    /** The dark-current spectrum - if any */
    std::unique_ptr<CSpectrum> m_darkCurrent;

    /** Remember how many spectra we have read from the scan */
    unsigned int m_specReadSoFarNum = 0U;

    /** True if the function 'CheckScanFile' has been called, and the scan-file handler
            has been initialized */
    bool m_initialized = false;

    /** The filename of the spectrum file */
    std::string m_fileName = "";

    /** The total number of spectra in the current .pak-file */
    std::uint32_t m_specNum = 0;

    /** An array containing the spectra in the current spectrum file.
        These are read in when 'CheckScanFile' is called and retrieved
        by GetSpectrum(...)
        The buffer is introduced to save some read/writes from hard-disk */
    std::vector<CSpectrum> m_spectrumBuffer;

    /** The number of spectra read in to the m_spectrumBuffer
        This might not be the same as 'm_specNum' */
    unsigned int m_spectrumBufferNum = 0;

    /** Updates the m_startTime and m_stopTime to include the timestamp of the provided spectrum */
    void UpdateStartAndStopTimeOfScan(novac::CSpectrum& spec);
};
}