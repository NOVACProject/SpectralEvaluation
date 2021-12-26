#pragma once

#include <vector>
#include <complex>

// ---------------------------------------------------------------------
// -------------------- Elementary function fitting --------------------
// ---------------------------------------------------------------------

namespace MathFit
{
    class CPolynomialFunction;
}

namespace novac
{
    /** This is a helper class used to optimize the fitting of a polynomial to data. */
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

        /** Fits a polynomial to the provided data.
        *   Notice that xData and yData must have the same length and not contain any repeated data points.
        *   @param xData The input data to the desired polynomial.
        *   @param yData The output data from the desired polynomial.
        *   @param polynomialCoefficients Will on successful return be filled with the coefficients of the polynomial, with the zero-th order coefficient first.
        *   @return True if the fit succeeds. */
        bool FitPolynomial(std::vector<double>& xData, std::vector<double>& yData, std::vector<double>& polynomialCoefficients);

        /** Specialization of FitPolynomial which fits a cubic polynomial to exactly four data points. */
        bool FitCubicPolynomial(std::vector<double>& xData, std::vector<double>& yData, std::vector<double>& polynomialCoefficients);

    };

    // -------------- Free functions for polynomial math --------------

    /** Calculates the value of the provided polynomial at the given point. */
    double PolynomialValueAt(const std::vector<double>& coefficients, double x);

    /** Finds the roots (zero crossings) of the provided polynomial. The result is returned in the second parameter.
    *   The number of roots equals the order of the polynomial minus one.
    *   This calculation only supports polynomials of order 0, 1, 2 or 3.
    *   @return true if the polynomial is supported and the root can be calculated, otherwise false. */
    bool FindRoots(const std::vector<double>& polynomial, std::vector<std::complex<double>>& roots);

}