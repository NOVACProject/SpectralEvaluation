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

double CGPSData::DoubleToAngle(double rawData)
{
    // conversion of values on the format 'ddmm.mmmmm'
    const int degree = static_cast<int>(rawData / 100.0); // extract the degrees part.

    const double minutes = fmod(rawData, 100.0); // extract the minutes part

    // combine the two
    return degree + minutes / 60.0;
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

double GpsMath::Distance(CGPSData position1, CGPSData position2)
{
    const double lat1 = position1.m_latitude * DEGREETORAD;
    const double lat2 = position2.m_latitude * DEGREETORAD;
    const double lon1 = position1.m_longitude * DEGREETORAD;
    const double lon2 = position2.m_longitude * DEGREETORAD;

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

    double tmpAngle = atan2(-std::sin(dLon) * std::cos(lat2),
        std::cos(lat1) * std::sin(lat2) - sin(lat1) * std::cos(lat2) * std::cos(dLon));

    if (tmpAngle < 0)
    {
        tmpAngle = TWO_PI + tmpAngle;
    }

    tmpAngle = RADTODEGREE * tmpAngle;
    return tmpAngle;
}

double GpsMath::Bearing(CGPSData position1, CGPSData position2)
{
    const double lat1 = position1.m_latitude * DEGREETORAD;
    const double lat2 = position2.m_latitude * DEGREETORAD;
    const double lon1 = position1.m_longitude * DEGREETORAD;
    const double lon2 = position2.m_longitude * DEGREETORAD;
    const double dLat = lat1 - lat2;
    const double dLon = lon1 - lon2;

    if ((dLon == 0) && (dLat == 0))
    {
        return 0;
    }

    double tmpAngle = atan2(-std::sin(dLon) * std::cos(lat2),
        std::cos(lat1) * std::sin(lat2) - std::sin(lat1) * std::cos(lat2) * std::cos(dLon));

    if (tmpAngle < 0)
    {
        tmpAngle = TWO_PI + tmpAngle;
    }

    tmpAngle = RADTODEGREE * tmpAngle;
    return tmpAngle;
}

void GpsMath::CalculateDestination(double lat1, double lon1, double dist, double az, double& lat2, double& lon2)
{

    const double dR = dist / R_Earth;

    // convert to radians
    lat1 = lat1 * DEGREETORAD;
    lon1 = lon1 * DEGREETORAD;
    az = az * DEGREETORAD;

    // calculate the second point
    lat2 = asin(std::sin(lat1) * std::cos(dR) + std::cos(lat1) * std::sin(dR) * std::cos(az));

    lon2 = lon1 + atan2(std::sin(az) * std::sin(dR) * std::cos(lat1), std::cos(dR) - std::sin(lat1) * std::sin(lat2));

    // convert back to degrees
    lat2 = lat2 * RADTODEGREE;
    lon2 = lon2 * RADTODEGREE;
}

CGPSData GpsMath::CalculateDestination(CGPSData origin, double dist, double az)
{
    const double dR = dist / R_Earth;

    // convert to radians
    const double lat1 = origin.m_latitude * DEGREETORAD;
    const double lon1 = origin.m_longitude * DEGREETORAD;
    az = az * DEGREETORAD;

    // calculate the second point
    const double lat2 = asin(std::sin(lat1) * std::cos(dR) + std::cos(lat1) * std::sin(dR) * std::cos(az));

    const double lon2 = lon1 + atan2(std::sin(az) * std::sin(dR) * std::cos(lat1), std::cos(dR) - std::sin(lat1) * std::sin(lat2));

    // convert back to degrees
    return CGPSData(lat2 * RADTODEGREE, lon2 * RADTODEGREE, 0.0);
}

} // namespace novac