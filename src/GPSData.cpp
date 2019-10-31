#include <SpectralEvaluation/GPSData.h>
#include <math.h>
#include <cstdlib>

CGPSData::CGPSData()
    : m_altitude(0.0), m_latitude(0.0), m_longitude(0.0)
{
}

CGPSData::CGPSData(const CGPSData &gps2)
{
    this->m_altitude = gps2.m_altitude;
    this->m_latitude = gps2.m_latitude;
    this->m_longitude = gps2.m_longitude;
}

CGPSData::CGPSData(double lat, double lon, double alt)
{
    this->m_altitude = alt;
    this->m_latitude = lat;
    this->m_longitude = lon;
}

CGPSData &CGPSData::operator =(const CGPSData &gps2)
{
    this->m_altitude = gps2.m_altitude;
    this->m_latitude = gps2.m_latitude;
    this->m_longitude = gps2.m_longitude;
    return *this;
}

bool CGPSData::operator==(const CGPSData &gps2) const
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

/* The GPS reports latitude and longitude in the format ddmm.mmmm
  , this function converts this to the format dd.dddd */
double CGPSData::DoubleToAngle(double rawData)
{
    double remainder    = fmod(rawData, 100.0);
    int degree          = (int)(rawData / 100);
    double fDegree      = degree + remainder / 60.0;

    return fDegree;
}