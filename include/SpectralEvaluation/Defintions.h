#pragma once

/** Not initialized values are set to 'NOT_A_NUMBER' to indicate that a value is missing.
This can e.g. be the case if only the wind-direction (and not the wind speed)
is known by an element in the database (which can be the case if a the
wind direction was calculated by combining two scans). */
#define NOT_A_NUMBER -9999.0

// converts degrees to radians
#define DEGREETORAD 0.017453 

// converts radians to degrees
#define RADTODEGREE 57.295791

// a quite familiar constant
#define TWO_PI 6.28318
#define HALF_PI 1.5708
#ifndef M_PI
    #define M_PI 3.141592
#endif
