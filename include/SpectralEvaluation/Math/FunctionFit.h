#pragma once

#include <vector>

// ---------------------------------------------------------------------
// -------------------- Elementary function fitting --------------------
// ---------------------------------------------------------------------

namespace MathFit
{
class CGaussFunction;
}

namespace novac
{
enum class FUNCTION_FIT_RETURN_CODE
{
    SUCCESS = 0,
    FIT_FAILURE,
    EMPTY_INPUT,
    MISSING_WAVELENGTH_CALIBRATION,
    INVALID_INPUT_DATA
};

/// <summary>
/// Fits a polynomial of the provided order to the given X- and Y-axis data, with the following constraints:
/// 1: The polynomial order must be at least zero.
/// 2: The two vectors xData and yData must have equal length.
/// 3: There must not be any repetions in xData and yData
/// 4: The length of xData (and yData) must be at least polynomialOrder + 1.
/// </summary>
/// <param name="polynomialOrder">The order of the polynomial to fit. Must be >= zero.</param>
/// <param name="xData">The x-axis data, i.e. the input data to the polynomial.</param>
/// <param name="yData">The y-axis data, i.e. the output data of the polynomial.</param>
/// <param name="polynomialCoefficients">Will on successful return be filled with the coefficients of the polynomial, with the 0th order coefficient first.</param>
/// <returns></returns>
FUNCTION_FIT_RETURN_CODE FitPolynomial(int polynomialOrder, std::vector<double>& xData, std::vector<double>& yData, std::vector<double>& polynomialCoefficients);

/// <summary>
/// Creates a basic initial estimate of a gaussian function which could fit to the provided data.
/// This does not do the actual fit but uses heuristics to get an initial estimation, making
/// the Gaussian fit itself more stable.
/// </summary>
void CreateInitialEstimate(const std::vector<double>& x, const std::vector<double>& y, MathFit::CGaussFunction& result);

/// <summary>
/// Fits a Gaussian function to the provided data.
/// </summary>
FUNCTION_FIT_RETURN_CODE FitGaussian(std::vector<double>& x, std::vector<double>& y, MathFit::CGaussFunction& gaussian);

/// <summary>
/// Fits a Gaussian function to the provided data.
/// The x-axis data will be index.
/// </summary>
FUNCTION_FIT_RETURN_CODE FitGaussian(std::vector<double>& y, MathFit::CGaussFunction& gaussian);

}