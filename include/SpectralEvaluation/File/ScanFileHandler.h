#pragma once

#include <vector>
#include <SpectralEvaluation/DateTime.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>

namespace novac
{

    class CGPSData;

    /** <b>CScanFileHandler</b> is a class to read in information from the scan-files
        (all the spectra from one scan are supposed to be packed together in one file in Manne's 'pak'-format.
         Each instance of 'CScanFileHandler' is capable of reading data from one .pak-file. */
    class CScanFileHandler
    {
    public:
        CScanFileHandler();

        // ----------------------------------------------------------------------
        // ---------------------- PUBLIC DATA -----------------------------------
        // ----------------------------------------------------------------------

        /** The serial number of the spectrometer which has collected the spectra
            in this scan */
        std::string m_device = "";

        /** The channel of the spectrometer which was used for collecting this scan.
            (if a SD2000 with multiple channels is used, one spectrometer should be
            configured for each channel). */
        unsigned char m_channel = 0;

        /** The time (UMT) when the measurement started */
        CDateTime m_startTime;

        /** The time (UMT) when the measurement was finished */
        CDateTime m_stopTime;

        /** If any error occurs in the reading of the file, this int is set to
            any of the errors defined int 'SpectrumIO.h. */
        int m_lastError;

        // ----------------------------------------------------------------------
        // --------------------- PUBLIC METHODS ---------------------------------
        // ----------------------------------------------------------------------

        /** Checks the scan saved in the given filename
            If any file-error occurs the parameter 'm_lastError' will be set.
            @param fileName - the name of the file in which the spectra of the scan are saved.
            @return true on success.
            @return false if any error occurs */
        bool CheckScanFile(const std::string& fileName);

        /** Gets the next spectrum in the scan.
            If any file-error occurs the parameter 'm_lastError' will be set.
            @param spec - will on successful return be filled with the newly read spectrum.
            @return the number of spectra read (0 if failure and 1 on success).*/
        int GetNextSpectrum(CSpectrum& spec);

        /** Returns the desired spectrum in the scan.
            If any file-error occurs the parameter 'm_lastError' will be set.
            @param spec - will on successful return be filled with the newly read spectrum.
            @param specNo - The zero-based index into the scan-file.
            @return the number of spectra read (1 if successful, otherwise 0) */
        int GetSpectrum(CSpectrum& spec, long specNo);

        /** Gets the dark spectrum of the scan
            @return 0 if there is a dark spectrum, else non-zero */
        int GetDark(CSpectrum& spec) const;

        /** Gets the sky spectrum of the scan
            @return 0 if there is a sky spectrum, else non-zero */
        int GetSky(CSpectrum& spec) const;

        /** Gets the offset spectrum of the scan - if any
            @return 0 if there is an offset spectrum, else non-zero */
        int GetOffset(CSpectrum& spec) const;

        /** Gets the dark-current spectrum of the scan - if any
            @return 0 if there is a dark-current spectrum, else non-zero */
        int GetDarkCurrent(CSpectrum& spec) const;

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
        const CGPSData& GetGPS() const;

        /** Retrieves compass-information from the spectrum files */
        double GetCompass() const;

        /** Retrieves the name of the file that this object is working on */
        const std::string& GetFileName() const { return m_fileName; }

        /** Resets the m_specReadSoFarNum to start reading from the first spectrum again */
        void  ResetCounter();

        /** Retrieves the total number of spectra in the .pak-file (including sky and dark) */
        int GetSpectrumNumInFile() const;

        /** Returns a spectrum which is the average of the provided indices.
            @return the number of spectra co-added (may be less than indices.size()
            if some spectrum/spectra could not be read). */
        int AddSpectra(const std::vector<size_t>& indices, CSpectrum& result);

    private:
        // ----------------------------------------------------------------------
        // ---------------------- PRIVATE DATA ----------------------------------
        // ----------------------------------------------------------------------

        /** The dark spectrum */
        CSpectrum m_dark;
        bool m_fHasDark;

        /** The sky spectrum */
        CSpectrum m_sky;
        bool m_fHasSky;

        /** The offset spectrum - if any */
        CSpectrum m_offset;
        bool m_fHasOffset;

        /** The dark-current spectrum - if any */
        CSpectrum m_darkCurrent;
        bool m_fHasDarkCurrent;

        /** Remember how many spectra we have read from the scan */
        unsigned int m_specReadSoFarNum;

        /** True if the function 'CheckScanFile' has been called, and the scan-file handler
                has been initialized */
        bool m_initialized;

        /** The filename of the spectrum file */
        std::string m_fileName;

        /** The total number of spectra in the current .pak-file */
        std::uint32_t m_specNum = 0;

        /** An array containing the spectra in the current spectrum file.
            These are read in when 'CheckScanFile' is called and retrieved
            by GetSpectrum(...)
            The buffer is introduced to save some read/writes from hard-disk */
        std::vector<CSpectrum> m_spectrumBuffer;

        /** The number of spectra read in to the m_spectrumBuffer
            This might not be the same as 'm_specNum' */
        unsigned int m_spectrumBufferNum;
    };
}