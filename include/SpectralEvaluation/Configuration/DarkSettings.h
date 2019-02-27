#pragma once

#include <string>

namespace Configuration
{
    // defining how to get the dark-spectrum
    enum class DARK_SPEC_OPTION { MEASURE, MODEL_SOMETIMES, MODEL_ALWAYS, DARK_USER_SUPPLIED };

    // defining how to use the dark-current and offset spectra
    enum class DARK_MODEL_OPTION { MEASURED, USER_SUPPLIED };

    struct CDarkSettings
    {
    public:
        CDarkSettings() = default;
        ~CDarkSettings() = default;

        /** Copying */
        CDarkSettings& operator=(const CDarkSettings&) = default;
        CDarkSettings(const CDarkSettings&) = default;

        /** Resets all values to default */
        void Clear()
        {
            m_darkSpecOption = DARK_SPEC_OPTION::MEASURE;
            m_offsetOption = DARK_MODEL_OPTION::MEASURED;
            m_darkCurrentOption = DARK_MODEL_OPTION::MEASURED;
            m_offsetSpec = "";
            m_darkCurrentSpec = "";
        }

        /** The options for the how to get the dark.
            Can be: MEASURED - use measured (DEFAULT)
                    1 - model of no measured is available
                    2 - always model
                    3 - the dark-spectrum is given by the user, do not model */
        DARK_SPEC_OPTION m_darkSpecOption = DARK_SPEC_OPTION::MEASURE;

        /** The offset-spectrum, only useful if 'm_darkSpecOption' is not 0.
                When this should be used is determined by 'm_offsetOption'.
                If 'm_darkSpecOption' is '3' then this is the dark-spectrum to use */
        std::string m_offsetSpec = "";

        /** The option for how to use the offset-spectrum.
            Can be:	MEASURED - always use measured
                    1 - use the user supplied */
        DARK_MODEL_OPTION m_offsetOption = DARK_MODEL_OPTION::MEASURED;

        /** The dark-current spectrum, only useful if 'm_darkSpecOption' is not 0.
                When this should be used is determined by 'm_darkCurrentOption'. */
        std::string m_darkCurrentSpec = "";

        /** The option for how to use the dark-current spectrum.
                Can be:	MEASURED - always use measured
                        1 - use the user supplied */
        DARK_MODEL_OPTION m_darkCurrentOption = DARK_MODEL_OPTION::MEASURED;
    };

}