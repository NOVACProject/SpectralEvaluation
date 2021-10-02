#include <SpectralEvaluation/Flux/ScanFluxResult.h>

using namespace novac;

ScanFluxResult::ScanFluxResult()
{
    m_flux = 0.0;
    m_fluxOk = true;
    m_windDirection = -999.0;
    m_windSpeed = -999.0;
    m_plumeHeight = -999.0;

    m_coneAngle = -999.0;
    m_tilt = -999.0;
    m_compass = -999.0;
    m_volcano = -1;
}

void ScanFluxResult::Clear()
{
    m_flux = 0.0;
    m_fluxOk = true;
    m_windDirection = -999.0;
    m_windSpeed = -999.0;
    m_plumeHeight = -999.0;

    m_coneAngle = -999.0;
    m_tilt = -999.0;
    m_compass = -999.0;
    m_volcano = -1;
}

