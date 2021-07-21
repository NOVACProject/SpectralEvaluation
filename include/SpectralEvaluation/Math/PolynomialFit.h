#pragma once

#include <vector>

// ---------------------------------------------------------------------
// -------------------- Elementary function fitting --------------------
// ---------------------------------------------------------------------

namespace MathFit
{
class CPolynomialFunction;
}

namespace novac
{

/// <summary>
/// This is a helper class used to optimize the fitting of a polynomial to data.
/// </summary>
class PolynomialFit
{
private:
    MathFit::CPolynomialFunction* functionToFit;

    const int polynomialOrder;

public:
    PolynomialFit(int order);
    virtual ~PolynomialFit();

    PolynomialFit(const PolynomialFit& other) = delete;
    PolynomialFit& operator=(const PolynomialFit& other) = delete;

    /// <summary>
    /// Fits a polynomial to the provided data. 
    /// Notice that xData and yData must have the same length and not contain any repeated data points.
    /// </summary>
    /// <param name="xData">The input data to the desired polynomial.</param>
    /// <param name="yData">The output data from the desired polynomial.</param>
    /// <param name="polynomialCoefficients">Will on successful return be filled with the coefficients of the polynomial, with the zero-th order coefficient first.s</param>
    /// <returns>True if the fit succeeds.</returns>
    bool FitPolynomial(std::vector<double>& xData, std::vector<double>& yData, std::vector<double>& polynomialCoefficients);

    // Specialization of FitPolynomial which fits a cubic polynomial to four data points.
    bool FitCubicPolynomial(std::vector<double>& xData, std::vector<double>& yData, std::vector<double>& polynomialCoefficients);

};

}