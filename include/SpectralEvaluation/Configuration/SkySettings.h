#pragma once

#include <string>

namespace Configuration
{
    // Defining the options for how to get the 'sky' spectrum
    enum class SKY_OPTION
    {
        MEASURED_IN_SCAN,
        AVERAGE_OF_GOOD_SPECTRA_IN_SCAN,
        SPECTRUM_INDEX_IN_SCAN,
        USER_SUPPLIED
    };

    struct CSkySettings
    {
    public:
        CSkySettings() = default;
        ~CSkySettings() = default;

        /** Copying */
        CSkySettings& operator=(const CSkySettings&) = default;
        CSkySettings(const CSkySettings&) = default;

        /** This is the main option for how to get the sky-spectrum */
        SKY_OPTION skyOption = SKY_OPTION::MEASURED_IN_SCAN;

        /** If skyOption Equals SPECTRUM_INDEX_IN_SCAN, then this points out
            the (zero-based) index in the scan to the spectrum to use as sky. */
        int indexInScan = 0;

        /** If skyOption equals USER_SUPPLIED, then this points out the file
            to get the sky-spectrum from. */
        std::string skySpectrumFile = "";
    };
}