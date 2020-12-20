#pragma once

class CDateTime;

// --------------------------------------------------------------------
// -------------------- SUN - FUNCTIONS -------------------------------
// --------------------------------------------------------------------

/** Retrieves the solar zenith angle (SZA) and the solar azimuth angle (SAZ)
        for the site specified by (lat, lon) and for the time given in gmtTime.
        Note that the returned angles are in degrees and that the specified
        time _must_ be GMT-time. */
void GetSunPosition(const CDateTime& gmtTime, double lat, double lon, double& SZA, double& SAZ);

