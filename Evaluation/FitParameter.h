#pragma once

namespace Evaluation
{
    /** constants for selection of fit-parameters*/
    enum FIT_PARAMETER{ 
        COLUMN, 
        COLUMN_ERROR, 
        SHIFT, 
        SHIFT_ERROR, 
        SQUEEZE, 
        SQUEEZE_ERROR,
        DELTA};

    // The options for how to ignore spectra (spectra which are too dark or saturated)
    enum IgnoreType { 
        IGNORE_DARK,
        IGNORE_LIMIT,
        IGNORE_NOTHING };

    typedef struct IgnoreOption {
        IgnoreType m_type;
        double     m_intensity;
        long       m_channel;
    } IgnoreOption;

}