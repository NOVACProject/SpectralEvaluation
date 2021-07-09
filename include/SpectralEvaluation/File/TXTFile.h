#pragma once

#include <string>

namespace novac
{
class CSpectrum;

/** <b>CTXTFile</b> is a simple class for reading/writing spectra
        from/to .txt/.xs - files */
class CTXTFile
{
public:

    /** Reads a spectrum from a TXT-file.
        @return true if the reading succeeds. */
    static bool ReadSpectrum(CSpectrum& spec, const std::string& fileName);

    /** Writes a spectrum to a TXT-file.
        If the spectrum contains a wavelength description then the file will be generated with two tab separated columns
            the first with the wavelength data and the second with the intensity data.
        @return true if the writing succeeds. */
    static bool WriteSpectrum(const CSpectrum& spec, const std::string& fileName);

    /** Writes a spectrum to a TXT-file.
        If the spectrum contains a wavelength description then the file will be generated with two tab separated columns
            the first with the wavelength data and the second with the intensity data.
        @return true if the writing succeeds. */
    static bool WriteSpectrum(const CSpectrum* spec, const std::string& fileName);

};
}