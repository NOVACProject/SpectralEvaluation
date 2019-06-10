#pragma once

// --------------------------------------------------------------------
// -------------------- CALCULATING FLUX ------------------------------
// --------------------------------------------------------------------

/** Calculates the flux for the supplied data using the old algorithm */
double CalculateFluxFlatScanner(const double *scanAngle, const double *column, double offset, int nDataPoints, double windSpeed, double windDirection, double relativePlumeHeight, double compass, double gasFactor);

/** Calculates the flux for the supplied data using the new algorithm */
double CalculateFluxConicalScanner(const double *scanAngle, const double *column, double offset, int nDataPoints, double windSpeed, double windDirection, double relativePlumeHeight, double compass, double gasFactor, double coneAngle, double tilt);

/** Calculates the flux for the Heidelberg-instrument for the supplied data using the general algorithm */
double CalculateFluxHeidelbergScanner(const double *scanAngle1, const double *scanAngle2, const double *column, double offset, int nDataPoints, double windSpeed, double windDirection, double relativePlumeHeight, double compass);

