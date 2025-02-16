#pragma once

#include <string>

namespace Configuration
{

// defining how to get the dark-spectrum
enum class DARK_SPEC_OPTION { MEASURED_IN_SCAN = 0, MODEL_WHEN_NOT_MEASURED_IN_SCAN = 1, MODEL_ALWAYS = 2, USER_SUPPLIED = 3 };

// defining how to use the dark-current and offset spectra
enum class DARK_MODEL_OPTION { MEASURED_IN_SCAN = 0, USER_SUPPLIED = 1 };

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
        m_darkSpecOption = DARK_SPEC_OPTION::MEASURED_IN_SCAN;
        m_offsetOption = DARK_MODEL_OPTION::MEASURED_IN_SCAN;
        m_darkCurrentOption = DARK_MODEL_OPTION::MEASURED_IN_SCAN;
        m_offsetSpec = "";
        m_darkCurrentSpec = "";
    }

    /** The options for the how to get the dark.
        Can be: MEASURED_IN_SCAN - use measured spectrum from the scan which is being evaluated (DEFAULT)
                MODEL_WHEN_NOT_MEASURED_IN_SCAN - model of no measured is available.
                MODEL_ALWAYS - always model, from the m_offsetSpec and m_darCurrentspec.
                USER_SUPPLIED - the dark-spectrum is given by the user, do not model (read from m_offsetSpec) */
    DARK_SPEC_OPTION m_darkSpecOption = DARK_SPEC_OPTION::MEASURED_IN_SCAN;

    /** The offset-spectrum.
            If m_darkSpecOption is MEASURED_IN_SCAN then this is ignored.
            If m_darkSpecOption is MODEL_ALWAYS and m_offsetOption is USER_SUPPLIED then this is the offset spectrum to model from.
            If m_darkSpecOption is USER_SUPPLIED then this is the dark-spectrum to use. */
    std::string m_offsetSpec = "";

    /** The option for how to use the offset-spectrum.
        Can be:	MEASURED_IN_SCAN - always use measured (DEFAULT)
                USER_SUPPLIED - use the user supplied */
    DARK_MODEL_OPTION m_offsetOption = DARK_MODEL_OPTION::MEASURED_IN_SCAN;

    /** The dark-current spectrum, only useful if 'm_darkSpecOption' is not MEASURED_IN_SCAN.
            When this should be used is determined by 'm_darkCurrentOption'. */
    std::string m_darkCurrentSpec = "";

    /** The option for how to use the dark-current spectrum.
            Can be:	MEASURED_IN_SCAN - always use measured (DEFAULT)
                    USER_SUPPLIED - use the user supplied */
    DARK_MODEL_OPTION m_darkCurrentOption = DARK_MODEL_OPTION::MEASURED_IN_SCAN;
};

}