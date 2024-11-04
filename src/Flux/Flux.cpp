#include <SpectralEvaluation/Flux/Flux.h>
#include <SpectralEvaluation/Definitions.h>

#include <assert.h>
#include <limits>
#include <cmath>
#include <vector>
#include <iostream>

double CalculateFluxFlatScanner(const double* scanAngle, const double* column, double offset, size_t nDataPoints, double windSpeed, double windDirection, double relativePlumeHeight, double compass, double gasFactor)
{
    if (nDataPoints <= 2)
    {
        return 0.0;
    }

    double totalFlux = 0;

    // the wind factor
    const double windFactor = std::abs(std::cos(DEGREETORAD * (windDirection - compass)));

    // now calculate the flux
    for (size_t i = 0; i < nDataPoints - 1; ++i)
    {
        if (std::abs(std::abs(scanAngle[i]) - 90.0) < 0.5)
            continue; // the distance-calculation has a singularity at +-90 degrees so just skip those points!
        if (std::abs(std::abs(scanAngle[i + 1]) - 90.0) < 0.5)
            continue; // the distance-calculation has a singularity at +-90 degrees so just skip those points!

        // The vertical columns
        const double VCD1 = (column[i] - offset) * std::cos(DEGREETORAD * scanAngle[i]);
        const double VCD2 = (column[i + 1] - offset) * std::cos(DEGREETORAD * scanAngle[i + 1]);

        // calculating the horisontal distance
        const double TAN1 = std::tan(DEGREETORAD * scanAngle[i]);
        const double TAN2 = std::tan(DEGREETORAD * scanAngle[i + 1]);
        const double distance = relativePlumeHeight * std::abs(TAN2 - TAN1);

        // The average vertical column
        const double avgVCD = (1E-6) * gasFactor * (VCD1 + VCD2) * 0.5;

        // The flux...
        const double partialFlux = distance * avgVCD * windSpeed * windFactor;

        totalFlux += partialFlux;
    }

    return std::abs(totalFlux);
}

double CalculateFluxConicalScanner(const double* scanAngle, const double* column, double offset, size_t nDataPoints, double windSpeed, double windDirection, double relativePlumeHeight, double compass, double gasFactor, double coneAngle, double tilt)
{
    if (std::isnan(windSpeed) || std::isnan(windDirection) || std::isnan(relativePlumeHeight) || std::isnan(offset) || std::isnan(compass) || std::isnan(tilt) || std::isnan(gasFactor))
    {
        std::cout << "Cannot calculate flux, Nan in incoming parameters." << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (coneAngle < 10.0)
    {
        std::cout << "Cannot calculate flux with an invalid cone angle (less than 10 degrees)" << std::endl;
        return std::numeric_limits<double>::quiet_NaN();
    }

    double flux = 0;

    // convert the angles to radians
    tilt *= DEGREETORAD;
    coneAngle *= DEGREETORAD;

    // local-data buffer to store the intermediate calculations
    std::vector<double> alpha(nDataPoints, 0.0);
    std::vector<double> scd(nDataPoints, 0.0);
    std::vector<double> columnCorrection(nDataPoints, 0.0);
    std::vector<double> x(nDataPoints, 0.0);
    std::vector<double> y(nDataPoints, 0.0);

    // Temporary variables, to do less computations
    const double tan_coneAngle = std::tan(coneAngle);
    const double sin_tilt = std::sin(tilt);
    const double cos_tilt = std::cos(tilt);

    assert(!std::isnan(coneAngle));
    assert(!std::isnan(sin_tilt));
    assert(!std::isnan(cos_tilt));

    // First prepare the buffers before we calculate anything
    for (size_t i = 0; i < nDataPoints - 1; ++i)
    {
        // The slant columns
        scd[i] = column[i] - offset;

        // The scan-angles, in radians
        alpha[i] = scanAngle[i] * DEGREETORAD;

        // cosine and sine of the scan-angle
        const double cos_alpha = std::cos(alpha[i]);
        const double sin_alpha = std::sin(alpha[i]);
        assert(!std::isnan(cos_alpha));
        assert(!std::isnan(sin_alpha));

        // Calculate the AMF in order to get vertical columns
        const double x_term = std::pow(cos_tilt / tan_coneAngle - cos_alpha * sin_tilt, 2);
        const double y_term = std::pow(sin_alpha, 2);
        const double divisor = std::pow(cos_alpha * cos_tilt + sin_tilt / tan_coneAngle, 2);
        const double columnAmplification = std::sqrt((x_term + y_term) / divisor + 1);
        columnCorrection[i] = 1 / columnAmplification;

        assert(!std::isnan(x_term));
        assert(!std::isnan(y_term));
        assert(!std::isnan(divisor));
        assert(!std::isnan(columnAmplification));
        assert(std::abs(columnAmplification) > 0.01);

        // Calculate the projections of the intersection points in the ground-plane
        const double commonDenominator = cos_alpha * cos_tilt + sin_tilt / tan_coneAngle;
        x[i] = (cos_tilt / tan_coneAngle - cos_alpha * sin_tilt) / commonDenominator;
        y[i] = (sin_alpha) / commonDenominator;
    }

    // Now make the actual flux-calculation
    for (size_t i = 0; i < nDataPoints - 2; ++i)
    {
        if (std::abs(std::abs(alpha[i]) - HALF_PI) < 1e-2 || std::abs(std::abs(alpha[i + 1]) - HALF_PI) < 1e-2)
        {
            continue;// This algorithm does not work very well for scanangles around +-90 degrees
        }

        // The average vertical column
        double avgVCD = (1e-6) * gasFactor * (scd[i] * columnCorrection[i] + scd[i + 1] * columnCorrection[i + 1]) * 0.5;

        // The horizontal distance
        double S = relativePlumeHeight * sqrt(pow(x[i + 1] - x[i], 2) + pow(y[i + 1] - y[i], 2));

        // The local compass-direction [radians] due to the curvature of the cone
        double coneCompass = atan2(y[i + 1] - y[i], x[i + 1] - x[i]);

        // The wind-factor 
        double windFactor = std::abs(std::sin(DEGREETORAD * (windDirection - compass) - coneCompass));

        // The partial flux
        double partialFlux = avgVCD * S * windSpeed * windFactor;

        // The total flux
        flux += partialFlux;
    }

    return std::abs(flux);
}

double CalculateFluxHeidelbergScanner(const double* scanAngle1, const double* scanAngle2, const double* column, double offset, size_t nDataPoints, double windSpeed, double windDirection, double relativePlumeHeight, double /*compass*/)
{
    double flux = 0;
    double partialFlux;

    // local-data buffer to store the intermediate calculations
    // double	*alpha  = new double[nDataPoints];

    double* elev = new double[nDataPoints];
    double* azim = new double[nDataPoints];

    double* scd = new double[nDataPoints];
    double* columnCorrection = new double[nDataPoints];
    double* x = new double[nDataPoints];
    double* y = new double[nDataPoints];

    // First prepare the buffers before we calculate anything
    for (size_t i = 0; i < nDataPoints - 1; ++i)
    {
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
    for (size_t i = 0; i < nDataPoints - 2; ++i)
    {
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

