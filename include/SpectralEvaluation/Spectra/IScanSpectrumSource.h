#pragma once

#include <vector>
#include <SpectralEvaluation/Spectra/Spectrum.h>

namespace novac
{
    /** IScanSpectrumSource is an interface by which it is possible to retrieve spectra, possibly from file, possibly from memory.
        All the spectra contained in the IScanSpectrumSource are from one single scan and hence collected using one single instrument. */
    class IScanSpectrumSource
    {
    public:
        /** Returns the desired spectrum in the scan. Notice the first spectra may be 'sky' or 'dark'.
            @param specNumber The zero-based index into the scan-file (including sky and dark).
            @param spec will on successful return be filled with the requested spectrum in the scan.
            @return zero if successful. */
        virtual int GetSpectrum(int specNumber, CSpectrum& spec) = 0;

        /** Retrieves the total number of spectra in the file (including sky and dark) */
        virtual int GetSpectrumNumInFile() const = 0;

        /** Retrieves the time when the first spectrum was collected */
        virtual CDateTime GetScanStartTime() const = 0;

        /** Retrieves the time when the last spectrum was collected */
        virtual CDateTime GetScanStopTime() const = 0;

        /** Retrieves the serial number of the device which collected this scan. */
        virtual std::string GetDeviceSerial() const = 0;

        /** Resets the counter associated with 'GetNextSpectrum' */
        virtual void ResetCounter() = 0;

        /** Retrieves the next measured spectrum in the scan, ignoring special spectra such as sky and dark.
            @return zero on success. */
        virtual int GetNextMeasuredSpectrum(CSpectrum& spec) = 0;

        /** Retrieves the measured sky spectrum in the scan and copies it to 'result'.
            @return zero if there is a sky spectum. */
        virtual int GetSky(CSpectrum& spec) const = 0;

        /** Retrieves the measured dark spectrum in the scan and copies it to 'result'.
            @return zero if there is a dark spectum. */
        virtual int GetDark(CSpectrum& result) const = 0;

        /** Gets the measured offset spectrum of the scan and copies it to 'result'.
            @return zero if there is an offset spectrum. */
        virtual int GetOffset(CSpectrum& spec) const = 0;

        /** Gets the measured dark-current spectrum of the scan and copies it to 'result'.
            @return zero if there is a dark-current spectrum. */
        virtual int GetDarkCurrent(CSpectrum& spec) const = 0;

        /** @return the name of the file which contains this scan, if any. */
        virtual std::string GetFileName() const = 0;
    };
}