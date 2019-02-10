#pragma once

    /** Not initialized values are set to 'NOT_A_NUMBER' to indicate that a value is missing.
    This can e.g. be the case if only the wind-direction (and not the wind speed)
    is known by an element in the database (which can be the case if a the
    wind direction was calculated by combining two scans). */
#define NOT_A_NUMBER -9999.0