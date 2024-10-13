#include <SpectralEvaluation/GPSData.h>
#include <SpectralEvaluation/Definitions.h>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <iostream>

namespace novac
{
CGPSData::CGPSData()
    : m_altitude(0.0), m_latitude(0.0), m_longitude(0.0)
{
}

CGPSData::CGPSData(const CGPSData& other) :
    m_altitude(other.m_altitude),
    m_latitude(other.m_latitude),
    m_longitude(other.m_longitude)
{
}

CGPSData::CGPSData(double lat, double lon, double alt) :
    m_altitude(alt),
    m_latitude(lat),
    m_longitude(lon)
{
}

CGPSData& CGPSData::operator =(const CGPSData& other)
{
    this->m_altitude = other.m_altitude;
    this->m_latitude = other.m_latitude;
    this->m_longitude = other.m_longitude;
    return *this;
}

bool CGPSData::operator==(const CGPSData& gps2) const
{
    if (std::abs(this->m_longitude - gps2.m_longitude) > 1e-4)
    {
        return false;
    }
    else if (std::abs(this->m_latitude - gps2.m_latitude) > 1e-4)
    {
        return false;
    }
    else if (std::abs(this->m_altitude - gps2.m_altitude) > 1e-4)
    {
        return false;
    }
    else
    {
        return true;
    }
}

std::ostream& operator<<(std::ostream& os, const CGPSData& gps)
{
    char buffer[128];
    sprintf(buffer, "(lat: %.4lf; lon: %.4lf; alt: %.1lf)", gps.m_latitude, gps.m_longitude, gps.m_altitude);

    os << std::string(buffer);

    return os;
}

/* The GPS reports latitude and longitude in the format ddmm.mmmm
  , this function converts this to the format dd.dddd */
double CGPSData::DoubleToAngle(double rawData)
{
    const double remainder = fmod(rawData, 100.0);
    const int degree = (int)(rawData / 100);
    const double fDegree = degree + remainder / 60.0;

    return fDegree;
}


CNamedLocation::CNamedLocation(double lat, double lon, double alt, const std::string& name)
    : CGPSData(lat, lon, alt)
{
    this->m_name = name;
}

CNamedLocation::CNamedLocation(const CNamedLocation& other) :
    CGPSData(other.m_latitude, other.m_longitude, other.m_altitude),
    m_name(other.m_name)
{
}

CNamedLocation& CNamedLocation::operator =(const CNamedLocation& other)
{
    this->m_altitude = other.m_altitude;
    this->m_latitude = other.m_latitude;
    this->m_longitude = other.m_longitude;
    this->m_name = other.m_name;
    return *this;
}

// R_Earth is the radius of the Earth, in meters.
static const double R_Earth = 6367000;

double GpsMath::Distance(double lat1, double lon1, double lat2, double lon2)
{
    lat1 = lat1 * DEGREETORAD;
    lat2 = lat2 * DEGREETORAD;
    lon1 = lon1 * DEGREETORAD;
    lon2 = lon2 * DEGREETORAD;

    const double dLon = lon2 - lon1;
    const double dLat = lat2 - lat1;

    if ((dLon == 0) && (dLat == 0))
    {
        return 0;
    }

    const double a = std::pow((std::sin(dLat / 2.0)), 2.0) + std::cos(lat1) * std::cos(lat2) * std::pow((std::sin(dLon / 2.0)), 2.0);
    const double c = 2 * std::asin(std::min(1.0, std::sqrt(a)));
    const double distance = R_Earth * c;

    return distance;
}

double GpsMath::Bearing(double lat1, double lon1, double lat2, double lon2)
{
    lat1 = lat1 * DEGREETORAD;
    lat2 = lat2 * DEGREETORAD;
    lon1 = lon1 * DEGREETORAD;
    lon2 = lon2 * DEGREETORAD;
    const double dLat = lat1 - lat2;
    const double dLon = lon1 - lon2;

    if ((dLon == 0) && (dLat == 0))
    {
        return 0;
    }

    double tmpAngle = atan2(-sin(dLon) * cos(lat2),
        cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(dLon));

    if (tmpAngle < 0)
    {
        tmpAngle = TWO_PI + tmpAngle;
    }

    tmpAngle = RADTODEGREE * tmpAngle;
    return tmpAngle;
}

void GpsMath::CalculateDestination(double lat1, double lon1, double dist, double az, double& lat2, double& lon2) {

    const double dR = dist / R_Earth;

    // convert to radians
    lat1 = lat1 * DEGREETORAD;
    lon1 = lon1 * DEGREETORAD;
    az = az * DEGREETORAD;

    // calculate the second point
    lat2 = asin(sin(lat1) * cos(dR) + cos(lat1) * sin(dR) * cos(az));

    lon2 = lon1 + atan2(sin(az) * sin(dR) * cos(lat1), cos(dR) - sin(lat1) * sin(lat2));

    // convert back to degrees
    lat2 = lat2 * RADTODEGREE;
    lon2 = lon2 * RADTODEGREE;
}

} // namespace novac