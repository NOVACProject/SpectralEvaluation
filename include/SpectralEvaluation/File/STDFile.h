#pragma once

#include <string>
#include <vector>

namespace novac
{

class CSpectrum;

/** <b>CSTDFile</b> is a simple class for reading/writing spectra
            from/to .std-files */
class CSTDFile
{
public:
    /** Additional information on the spectrum which may be saved to file */
    struct ExtendedFormatInformation
    {
        // The location of the marker, in pixels
        double Marker = 0.0;

        // The minimum pixel to display
        int MinChannel = 0;

        // The maximum pixel to display. If negative then the spectrum-length will be used
        int MaxChannel = -1;

        // The minimum pixel in the math region
        int MathLow = 0;

        // The maximum pixel in the math region. If negative then the spectrum-length will be used
        int MathHigh = -1;

        // If set, then the provided calibration will be written to the file as well.
        // This should be stored with the 0th order coefficient first, followed by the 1st order coefficient, etc
        std::vector<double> calibrationPolynomial;

        // Listing of any additional (non-default) properties to write to the file
        std::vector<std::pair<std::string, std::string>> additionalProperties;
    };

    /** Reads a spectrum from a STD-file.
        @return true if the reading succeeds, else returns false. */
    static bool ReadSpectrum(CSpectrum& spec, const std::string& fileName);

    /** Writes a spectrum to a STD-file.
        @return true if the writing succeeds, else returns false. */
    static bool WriteSpectrum(const CSpectrum& spec, const std::string& fileName, int extendedFormat = 0);

    /** Writes a spectrum to a STD-file.
        @return true if the writing succeeds, else returns false. */
    static bool WriteSpectrum(const CSpectrum* spec, const std::string& fileName, int extendedFormat = 0);

    /** Writes a spectrum to a STD-file in extended format with the additional supplied information.
        @return true if the writing succeeds, else returns false. */
    static bool WriteSpectrum(const CSpectrum& spec, const std::string& fileName, const ExtendedFormatInformation& extendedInformation);
};
}