#include <SpectralEvaluation/Flux/Flux.h>
#include <SpectralEvaluation/Defintions.h>

#include <cmath>

double CalculateFluxFlatScanner(const double *scanAngle, const double *column, double offset, int nDataPoints, double windSpeed, double windDirection, double relativePlumeHeight, double compass)
{
    double avgVCD, VCD1, VCD2, TAN1, TAN2, distance;
    double flux = 0;
    double partialFlux;

    // the wind factor
    double windFactor = std::abs(std::cos(DEGREETORAD*(windDirection - compass)));

    // now calculate the flux
    for (int i = 0; i < nDataPoints - 1; ++i) {
        if (std::abs(std::abs(scanAngle[i]) - 90.0) < 0.5)
            continue; // the distance-calculation has a singularity at +-90 degrees so just skip those points!
        if (std::abs(std::abs(scanAngle[i + 1]) - 90.0) < 0.5)
            continue; // the distance-calculation has a singularity at +-90 degrees so just skip those points!

                      // The vertical columns
        VCD1 = (column[i] - offset)	 * std::cos(DEGREETORAD * scanAngle[i]);
        VCD2 = (column[i + 1] - offset) * std::cos(DEGREETORAD * scanAngle[i + 1]);

        // calculating the horisontal distance
        TAN1 = std::tan(DEGREETORAD * scanAngle[i]);
        TAN2 = std::tan(DEGREETORAD * scanAngle[i + 1]);
        distance = relativePlumeHeight * std::abs(TAN2 - TAN1);

        // The average vertical column
        avgVCD = (VCD1 + VCD2) * 0.5;

        // The flux...
        partialFlux = distance * avgVCD * windSpeed * windFactor;
        flux += partialFlux;
    }

    return std::abs(flux);
}

double CalculateFluxConicalScanner(const double *scanAngle, const double *column, double offset, int nDataPoints, double windSpeed, double windDirection, double relativePlumeHeight, double compass, double coneAngle, double tilt)
{
    double flux = 0;
    double partialFlux, columnAmplification;

    // convert the angles to radians
    tilt *= DEGREETORAD;
    coneAngle *= DEGREETORAD;

    // local-data buffer to store the intermediate calculations
    double	*alpha = new double[nDataPoints];
    double	*scd = new double[nDataPoints];
    double	*columnCorrection = new double[nDataPoints];
    double	*x = new double[nDataPoints];
    double	*y = new double[nDataPoints];

    // Temporary variables, to do less computations
    double	tan_coneAngle = std::tan(coneAngle);
    double	sin_tilt = std::sin(tilt);
    double	cos_tilt = std::cos(tilt);

    // First prepare the buffers before we calculate anything
    for (int i = 0; i < nDataPoints - 1; ++i) {
        // The slant columns
        scd[i] = column[i] - offset;

        // The scan-angles, in radians
        alpha[i] = scanAngle[i] * DEGREETORAD;

        // std::cosine and sine of the scan-angle
        double cos_alpha = std::cos(alpha[i]);
        double sin_alpha = std::sin(alpha[i]);

        // Calculate the AMF in order to get vertical columns
        double x_term = pow(cos_tilt / tan_coneAngle - cos_alpha*sin_tilt, 2);
        double y_term = pow(sin_alpha, 2);
        double divisor = pow(cos_alpha*cos_tilt + sin_tilt / tan_coneAngle, 2);
        columnAmplification = sqrt((x_term + y_term) / divisor + 1);
        columnCorrection[i] = 1 / columnAmplification;

        // Calculate the projections of the intersection points in the ground-plane
        double commonDenominator = cos_alpha*cos_tilt + sin_tilt / tan_coneAngle;
        x[i] = (cos_tilt / tan_coneAngle - cos_alpha*sin_tilt) / commonDenominator;
        y[i] = (sin_alpha) / commonDenominator;
    }

    // Now make the actual flux-calculation
    for (int i = 0; i < nDataPoints - 2; ++i) {
        if (std::abs(std::abs(alpha[i]) - HALF_PI) < 1e-2 || std::abs(std::abs(alpha[i + 1]) - HALF_PI) < 1e-2)
            continue;// This algorithm does not work very well for scanangles around +-90 degrees

                     // The average vertical column
        double avgVCD = (scd[i] * columnCorrection[i] + scd[i + 1] * columnCorrection[i + 1]) * 0.5;

        // The horizontal distance
        double S = relativePlumeHeight * sqrt(pow(x[i + 1] - x[i], 2) + pow(y[i + 1] - y[i], 2));

        // The local compass-direction [radians] due to the curvature of the cone
        double coneCompass = atan2(y[i + 1] - y[i], x[i + 1] - x[i]);

        // The wind-factor 
        double windFactor = std::abs(std::sin(DEGREETORAD * (windDirection - compass) - coneCompass));

        // The partial flux
        partialFlux = avgVCD * S * windSpeed * windFactor;

        // The total flux
        flux += partialFlux;
    }

    delete[] alpha;
    delete[] scd;
    delete[] columnCorrection;
    delete[] x;
    delete[] y;

    return std::abs(flux);
}

double CalculateFluxHeidelbergScanner(const double *scanAngle1, const double *scanAngle2, const double *column, double offset, int nDataPoints, double windSpeed, double windDirection, double relativePlumeHeight, double /*compass*/)
{
    double flux = 0;
    double partialFlux;

    // local-data buffer to store the intermediate calculations
    // double	*alpha							= new double[nDataPoints];

    double	*elev = new double[nDataPoints];
    double	*azim = new double[nDataPoints];

    double	*scd = new double[nDataPoints];
    double	*columnCorrection = new double[nDataPoints];
    double	*x = new double[nDataPoints];
    double	*y = new double[nDataPoints];

    // First prepare the buffers before we calculate anything
    for (int i = 0; i < nDataPoints - 1; ++i) {
        // The slant columns
        scd[i] = column[i + 1] - offset;

        // The scan-angles, in radians
        elev[i] = scanAngle1[i] * DEGREETORAD;
        azim[i] = scanAngle2[i] * DEGREETORAD;

        double tan_elev = std::tan(elev[i]);
        double cos_elev = std::cos(elev[i]);
        double sin_azim = std::sin(azim[i]);
        double cos_azim = std::cos(azim[i]);

        // Calculate the AMF in order to get vertical columns
        columnCorrection[i] = cos_elev;


        // Calculate the projections of the intersection points in the ground-plane
        double x_term = tan_elev * cos_azim;
        double y_term = tan_elev * sin_azim;

        x[i] = x_term;
        y[i] = y_term;
    }

    // Now make the actual flux-calculation
    /*TODO: flux calculations for Heidelberg instrument differ from Gothenborg instrument because local and global coordinate system do not differ!!!! Define another loop for Heidelberg?*/
    for (int i = 0; i < nDataPoints - 2; ++i) {
        if (std::abs(std::abs(elev[i]) - HALF_PI) < 1e-2 || std::abs(std::abs(elev[i + 1]) - HALF_PI) < 1e-2)
            continue;// This algorithm does not work very well for scanangles around +-90 degrees

                     // The average vertical column
        double avgVCD = (scd[i] * columnCorrection[i] + scd[i + 1] * columnCorrection[i + 1]) * 0.5;

        // The horizontal distance
        double S = relativePlumeHeight * sqrt(pow(x[i + 1] - x[i], 2) + pow(y[i + 1] - y[i], 2));

        // The local compass-direction [radians] due to the curvature of the cone
        double DirectionCompass = atan2(y[i + 1] - y[i], x[i + 1] - x[i]);

        // The wind-factor HD: compass=azim[i]
        double windFactor = std::abs(std::sin(DEGREETORAD * windDirection - DirectionCompass));

        // The partial flux
        partialFlux = avgVCD * S * windSpeed * windFactor;

        // The total flux
        flux += partialFlux;

    }
    //HD constants are to be deleted as well
    delete[] scd;
    delete[] columnCorrection;
    delete[] x;
    delete[] y;
    delete[] elev;
    delete[] azim;

    return std::abs(flux);
}

