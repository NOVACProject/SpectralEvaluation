#pragma once

#include <SpectralEvaluation/DateTime.h>

namespace novac
{

/** The class ScanFluxResult collects necessary data for storing the result of a
*   flux calculation from a scanner. All non-initialized variables are set to -999

    This is an attempt to collext the data used by NovacProgram and NovacPPP for
    representing the calculated flux from a scan. This class is extended with more
    information in both NovacProgram and in NovacPPP. */

class ScanFluxResult
{
public:
    ScanFluxResult();

    /** Clears the results */
    virtual void Clear();

    /** Copying */
    ScanFluxResult& operator=(const ScanFluxResult& other) = default;
    ScanFluxResult(const ScanFluxResult& other) = default;

    /** The calculated flux, in kg/s */
    double m_flux;

    /** True if the flux-value is a good measurement */
    bool m_fluxOk;

    /** The wind-direction used to calculate the flux */
    double m_windDirection;

    /** The wind-speed used to calculate the flux */
    double m_windSpeed;

    /** The plume-height used to calculate the flux */
    double m_plumeHeight;

    /** The cone-angle of the scanner that collected this scan */
    double m_coneAngle;

    /** The tilt of the scanner that collected this scan */
    double m_tilt;

    /** The compass-direction of the scanner that collected this scan */
    double m_compass;

    /** The date and time (UTC) when the measurement was started */
    novac::CDateTime m_startTime;

    /** The volcano that this measurement was made at. Set to -1 if unknown */
    int   m_volcano;
};

}
