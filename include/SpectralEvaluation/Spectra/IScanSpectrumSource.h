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
        /** Returns the desired spectrum in the scan.
            @param specNumber - The zero-based index into the scan-file.
            @param spec - will on successful return be filled with the requested spectrum in the scan.
            @return zero if successful. */
        virtual int GetSpectrum(int specNumber, CSpectrum& spec) = 0;

        /** Retrieves the measured sky spectrum in the scan and copies it to 'result'.
            @return zero if there is a sky spectum. */
        virtual int GetSky(CSpectrum& spec) = 0;

        /** Retrieves the measured dark spectrum in the scan and copies it to 'result'.
            @return zero if there is a dark spectum. */
        virtual int GetDark(CSpectrum& result) = 0;

        /** Gets the measured offset spectrum of the scan and copies it to 'result'.
            @return zero if there is an offset spectrum. */
        virtual int GetOffset(CSpectrum& spec) = 0;

        /** Gets the measured dark-current spectrum of the scan and copies it to 'result'.
            @return zero if there is a dark-current spectrum. */
        virtual int GetDarkCurrent(CSpectrum& spec) = 0;

        /** @return the name of the file which contains this scan, if any. */
        virtual std::string GetFileName() const = 0;
    };
}