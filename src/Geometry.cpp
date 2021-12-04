#include <SpectralEvaluation/DateTime.h>
#include <SpectralEvaluation/Geometry.h>
#include <SpectralEvaluation/Definitions.h>
#include <cmath>

namespace novac
{
// The following functions are taken from the NovacProgram.
//  Be sure to document further what is being done here.

/// <summary>
/// Calculates and returns the Right Ascension, declination and Equation of Time
/// See also https://en.wikipedia.org/wiki/Position_of_the_Sun
/// </summary>
/// <param name="D"></param>
/// <param name="RA"></param>
/// <param name="dec"></param>
/// <param name="EQT"></param>
void EquatorialCoordinates(double D, double& RA, double& dec, double& EQT)
{
    const double g_deg = fmod(357.529 + 0.98560028 * D, 360.0);
    const double g_rad = g_deg * DEGREETORAD;
    const double q_deg = fmod(280.459 + 0.98564736 * D, 360.0);

    const double L_deg = q_deg + 1.915 * std::sin(g_rad) + 0.02 * std::sin(2 * g_rad);
    const double L_rad = L_deg * DEGREETORAD;

    // The distance between the sun and the earth (in Astronomical Units)
    // const double R = 1.00014 - 0.01671 * cos(g_rad) - 0.00014 * cos(2 * g_rad);

    // The obliquity of the earth's orbit:
    const double obliq_deg = 23.439 - 0.00000036 * D;
    const double obliq_rad = obliq_deg * DEGREETORAD;

    // The right ascension (RA)
    double RA_rad = std::atan(std::cos(obliq_rad) * std::sin(L_rad) / std::cos(L_rad));
    if (RA_rad < 0)
        RA_rad = TWO_PI + RA_rad;

    if (fabs(RA_rad - L_rad) > 1.570796)
        RA_rad = M_PI + RA_rad;

    const double dec_rad = std::asin(std::sin(obliq_rad) * std::sin(L_rad));
    RA = fmod(RA_rad * RADTODEGREE, 360.0); // The right ascension

    // The declination
    dec = dec_rad * RADTODEGREE;

    // The Equation of Time
    EQT = q_deg / 15.0 - RA / 15.0;
}

void HorizontalCoordinates(double lat, double H, double dec, double& elev, double& azim)
{
    double H_rad = H * DEGREETORAD;
    double lat_rad = lat * DEGREETORAD;
    double dec_rad = dec * DEGREETORAD;

    // The elevation angle
    double elev_rad = asin(cos(H_rad) * cos(dec_rad) * cos(lat_rad) + sin(dec_rad) * sin(lat_rad));

    // The cosine of the azimuth - angle
    double cazim_rad = (cos(H_rad) * cos(dec_rad) * sin(lat_rad) - sin(dec_rad) * cos(lat_rad)) / cos(elev_rad);

    // The sine of the azimuth - angle
    double sazim_rad = (sin(H_rad) * cos(dec_rad)) / cos(elev_rad);

    double azim_rad = 0.0;
    if (cazim_rad > 0 && sazim_rad > 0)
        azim_rad = asin(sazim_rad);						// azim is in the range 0 - 90 degrees
    else if (cazim_rad < 0 && sazim_rad > 0)
        azim_rad = M_PI - asin(sazim_rad);			// azim is in the range 90 - 180 degrees
    else if (cazim_rad < 0 && sazim_rad < 0)
        azim_rad = M_PI - asin(sazim_rad);		// azim is in the range 180 - 270 degrees
    else if (cazim_rad > 0 && sazim_rad < 0)
        azim_rad = TWO_PI + asin(sazim_rad);		// azim is in the range 270 - 360 degrees

    elev = elev_rad * RADTODEGREE;
    azim = azim_rad * RADTODEGREE;
}

/** Returns the hour angle given the longitude and equation of time. */
double GetHourAngle(double hr, double lon, double EqT)
{
    return 15.0 * (hr + lon / 15 + EqT - 12);
}

void GetSunPosition(const CDateTime& gmtTime, double lat, double lon, double& SZA, double& SAZ)
{
    SZA = 0.0;
    SAZ = 0.0;

    // Get the julian day
    double D = JulianDay(gmtTime) - 2451545.0;

    // Get the Equatorial coordinates...
    double	RA; //	the right ascension (deg)
    double	dec; // the declination	(deg)
    double	EqT;	// the equation of time (hours)
    EquatorialCoordinates(D, RA, dec, EqT);

    // Get the hour angle
    double fractionalHour = (double)gmtTime.hour + gmtTime.minute / 60.0 + gmtTime.second / 3600.0;
    double H = GetHourAngle(fractionalHour, lon, EqT);

    // Get the horizontal coordinates
    double	elev, sAzim; // The elevation and azimuth (towards south);
    HorizontalCoordinates(lat, H, dec, elev, sAzim);

    // Convert the elevation into sza
    SZA = 90.0 - elev;

    // Convert the azimuth to a value counted from the north and 
    SAZ = fmod(180.0 + sAzim, 360.0);
}
}
