#pragma once

#include "../Spectra/DateTime.h"
#include "WindField.h"
#include "PlumeHeight.h"

enum FLUX_QUALITY_FLAG {
    FLUX_QUALITY_GREEN,
    FLUX_QUALITY_YELLOW,
    FLUX_QUALITY_RED
};

/** The class <b>CFluxResult</b> is a generic class for storing the results
        from a flux-calculation of a scan. The class holds the values of all the
        parameters used in the calculation (except for the measurment itself) and
        the result of the measurement. All non-initialized variables are set to -999
        */

namespace Flux
{
    class CFluxResult
    {
    public:
        CFluxResult();
        ~CFluxResult();

        /** Clears the results */
        void Clear();

        /** Assignment operator */
        CFluxResult &operator=(const CFluxResult &fl2);

        /** The calculated flux, in kg/s */
        double	m_flux;

        /** The quality of this flux measurement */
        FLUX_QUALITY_FLAG m_fluxQualityFlag;

        /** The error in flux due to the uncertainty in
        wind speed and wind direction */
        double m_fluxError_Wind;

        /** The error in flux due to the uncertainty in plume height */
        double m_fluxError_PlumeHeight;

        /** The wind field that was used to calculate this flux */
        CWindField m_windField;

        /** The plume height that was used to calculate this flux */
        CPlumeHeight m_plumeHeight;

        /** The number of good spectra in this measurement. 
            This is alsothe number of column-values that were used to calculate the flux */
        int m_numGoodSpectra;

        /** The cone-angle of the scanner that collected this scan */
        double m_coneAngle;

        /** The tilt of the scanner that collected this scan */
        double m_tilt;

        /** The compass-direction of the scanner that collected this scan */
        double m_compass;

        /** The date and time (UTC) when the measurement was started */
        CDateTime m_startTime;

        /** The date and time (UTC) when the measurement was finished */
        CDateTime m_stopTime;

        /** The volcano that this measurement was made at. Set to -1 if unknown */
        int m_volcano;

        /** The instrument that collected this scan */
        std::string  m_instrument;

        /** The type of the instrument that collected this scan */
        // INSTRUMENT_TYPE m_instrumentType;

        /** The calculated offset of the scan */
        double  m_scanOffset;

        /** The calculated plume completeness */
        double  m_completeness;

        /** The calculated plume centre position */
        double  m_plumeCentre[2];
    };
}