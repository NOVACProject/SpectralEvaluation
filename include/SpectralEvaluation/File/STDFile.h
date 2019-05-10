#pragma once

#include <string>

class CSpectrum;

namespace SpectrumIO
{
    /** <b>CSTDFile</b> is a simple class for reading/writing spectra
            from/to .std-files */
    class CSTDFile
    {
    public:
        CSTDFile(void);
        ~CSTDFile(void);

        /** Reads a spectrum from a STD-file. @return true if the reading succeeds. */
        static bool ReadSpectrum(CSpectrum &spec, const std::string &fileName);

        /** Writes a spectrum to a STD-file. @return true if the writing succeeds. */
        static bool WriteSpectrum(const CSpectrum &spec, const std::string &fileName, int extendedFormat = 0);

        /** Writes a spectrum to a STD-file. @return true if the writing succeeds. */
        static bool WriteSpectrum(const CSpectrum *spec, const std::string &fileName, int extendedFormat = 0);
    };
}