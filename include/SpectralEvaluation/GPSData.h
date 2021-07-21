#pragma once

#include <string>

namespace novac
{

/** CGPSData represents a location on the Earth */
class CGPSData
{
public:
    CGPSData();

    /** Make a new CGPSData object with the supplied latitude, longitude and altitude */
    CGPSData(double lat, double lon, double alt);

    CGPSData(const CGPSData& gps2);
    CGPSData& operator=(const CGPSData& gps2);

    bool operator==(const CGPSData& gps2) const;

    ~CGPSData() = default;

    /** The altitude from the GPS data, in meters above sea-level */
    double  m_altitude = 0.0;

    /** The latitude, in degrees (with decimals) with positive values on the northern hemisphere */
    double  m_latitude = 0.0;

    /** The longitude, in degrees (with decimals) with positive values on the easter hemisphere */
    double  m_longitude = 0.0;

    /* The GPS reports latitude and longitude in the format ddmm.mmmm
            this function converts this to the format dd.dddd */
    static double DoubleToAngle(double rawData);
};

/** CNamedLocation corresponds to a location with a given name.
    Commonly used to pin-point the peak of a volcano (with the name of the volcano). */
class CNamedLocation : public CGPSData
{
public:
    CNamedLocation() : CGPSData() { }

    CNamedLocation(double lat, double lon, double alt, const std::string& name);

    /** The name given to this location */
    std::string m_name;

    CNamedLocation(const CNamedLocation& gps2);
    CNamedLocation& operator=(const CNamedLocation& gps2);

};
}