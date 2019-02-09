#pragma once

class CGPSData
{
public:
    CGPSData();

    /** Make a new CGPSData object with the supplied latitude, longitude and altitude */
    CGPSData(double lat, double lon, double alt);

    CGPSData(const CGPSData &gps2);
    CGPSData &operator=(const CGPSData &gps2);

    ~CGPSData();

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
