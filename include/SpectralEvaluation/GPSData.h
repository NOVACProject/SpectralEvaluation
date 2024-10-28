#pragma once

#include <string>
#include <iosfwd>

namespace novac
{

/** CGPSData represents a location on the Earth */
class CGPSData
{
public:
    CGPSData();
    CGPSData(double lat, double lon, double alt);
    CGPSData(const CGPSData& gps2);
    CGPSData& operator=(const CGPSData& gps2);

    bool operator==(const CGPSData& gps2) const;

    friend std::ostream& operator<<(std::ostream& os, const CGPSData& gps);

    ~CGPSData() = default;

    /** The altitude from the GPS data, in meters above sea-level */
    double m_altitude = 0.0;

    /** The latitude, in degrees (with decimals) with positive values on the northern hemisphere */
    double m_latitude = 0.0;

    /** The longitude, in degrees (with decimals) with positive values on the easter hemisphere */
    double m_longitude = 0.0;

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

/* GpsMath contains helper methods for performing calculations on Gps Data */
class GpsMath
{
public:

    /* This returns the distance in METERS between the two points defined
        by (lat1,lon1) and (lat2, lon2). ALL ANGLES MUST BE IN DEGREES. */
    static double Distance(double lat1, double lon1, double lat2, double lon2);

    /* This returns the distance in METERS between the two points.
        ALL ANGLES MUST BE IN DEGREES. */
    static double Distance(CGPSData position1, CGPSData position2);

    /* This function returns the initial bearing (<b>degrees</b>) when travelling from
      the point defined by (lat1, lon1) to the point (lat2, lon2). <b>All angles must be in degrees</b> */
    static double Bearing(double lat1, double lon1, double lat2, double lon2);

    /* This function returns the initial bearing (<b>degrees</b>) when travelling from
      the point defined by (lat1, lon1) to the point (lat2, lon2). <b>All angles must be in degrees</b> */
    static double Bearing(CGPSData position1, CGPSData position2);

    /** This function calculates the latitude and longitude for point
            which is the distance 'dist' m and bearing 'az' degrees from
            the point defied by 'lat1' and 'lon1' */
    static void CalculateDestination(double lat1, double lon1, double dist, double az, double& lat2, double& lon2);

    /** This function calculates the latitude and longitude for point
            which is the distance 'dist' m and bearing 'az' degrees from the origin point */
    static CGPSData CalculateDestination(CGPSData origin, double dist, double az);
};


}