#pragma once

// --------------------------------------------------------------------
// -------------------- CALCULATING FLUX ------------------------------
// --------------------------------------------------------------------

/** Calculates the flux for the supplied data using the old algorithm.
    Notice that the angles here must be in degrees, with zero being zenith.
    The column and offset must have the same unit (ppmm or molec/cm2).
    The supplied gasfactor is used to convert from column unit (ppmm or molec/cm2) to weight.
    Returns the calculated flux, the unit of this depends on the unit of the columns and the gas factor. */
double CalculateFluxFlatScanner(const double* scanAngle, const double* column, double offset, size_t nDataPoints, double windSpeed, double windDirection, double relativePlumeHeight, double compass, double gasFactor);

/** Calculates the flux for the supplied data using the new algorithm
    Notice that the angles here must be in degrees, with zero being zenith.
    The column and offset must have the same unit (ppmm or molec/cm2).
    The supplied gasfactor is used to convert from column unit (ppmm or molec/cm2) to weight.
    Returns the calculated flux, the unit of this depends on the unit of the columns and the gas factor. */
double CalculateFluxConicalScanner(const double* scanAngle, const double* column, double offset, size_t nDataPoints, double windSpeed, double windDirection, double relativePlumeHeight, double compass, double gasFactor, double coneAngle, double tilt);

/** Calculates the flux for the Heidelberg-instrument for the supplied data using the general algorithm
    Notice that the angles here must be in degrees, with zero being zenith.
    The column and offset must have the same unit (ppmm or molec/cm2).
    The supplied gasfactor is used to convert from column unit (ppmm or molec/cm2) to weight.
    Returns the calculated flux, the unit of this depends on the unit of the columns and the gas factor. */
double CalculateFluxHeidelbergScanner(const double* scanAngle1, const double* scanAngle2, const double* column, double offset, size_t nDataPoints, double windSpeed, double windDirection, double relativePlumeHeight, double compass);

