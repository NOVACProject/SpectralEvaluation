#pragma once

#include <SpectralEvaluation/Units.h>
#include <SpectralEvaluation/GPSData.h>

namespace novac
{

class CDateTime;

// --------------------------------------------------------------------
// -------------------- SUN - FUNCTIONS -------------------------------
// --------------------------------------------------------------------

struct SolarPosition
{
    // degrees from zenith
    angle_degrees_t zenithAngle;

    // azumith angle
    angle_degrees_t azimuth;
};

/** Retrieves the solar zenith angle (SZA) and the solar azimuth angle (SAZ)
        for the site specified by (lat, lon) and for the time given in gmtTime.
        Note that the returned angles are in degrees and that the specified
        time _must_ be GMT-time.
    The returned angles are in degrees. */
SolarPosition GetSunPosition(const CDateTime& gmtTime, CGPSData position);

}
