#include "catch.hpp"
#include <SpectralEvaluation/Math/PolynomialFit.h>

using namespace novac;

TEST_CASE("PolynomialFit - PolynomialValueAt", "[Math][PolynomialFit]")
{
    SECTION("Constant")
    {
        std::vector<double> polynomial{ 3.4 };

        REQUIRE(3.4 == PolynomialValueAt(polynomial, 1.0));
        REQUIRE(3.4 == PolynomialValueAt(polynomial, -1.0));
        REQUIRE(3.4 == PolynomialValueAt(polynomial, 123.2));
        REQUIRE(3.4 == PolynomialValueAt(polynomial, -1922.212));
    }

    SECTION("first order polynomial")
    {
        std::vector<double> polynomial{ 3.4, 1.2 }; // i.e. y = 3.4 + 1.2x

        REQUIRE(3.4 == PolynomialValueAt(polynomial, 0.0));
        REQUIRE(4.6 == PolynomialValueAt(polynomial, 1.0));
        REQUIRE(2.2 == PolynomialValueAt(polynomial, -1.0));
        REQUIRE(123.4 == PolynomialValueAt(polynomial, 100));
        REQUIRE(-116.6 == PolynomialValueAt(polynomial, -100));
    }

    SECTION("second order polynomial")
    {
        std::vector<double> polynomial{ 1.0, 2.0, 3.0 }; // i.e. y = 1 + 2x + 3x2

        REQUIRE(1.0 == PolynomialValueAt(polynomial, 0.0));
        REQUIRE(6.0 == PolynomialValueAt(polynomial, 1.0));
        REQUIRE(2.0 == PolynomialValueAt(polynomial, -1.0));
        REQUIRE(1 + 200 + 30000 == PolynomialValueAt(polynomial, 100));
        REQUIRE(1 - 200 + 30000 == PolynomialValueAt(polynomial, -100));
    }

    SECTION("third order polynomial")
    {
        std::vector<double> polynomial{ 1.0, 2.0, 3.0, 4.0 }; // i.e. y = 1 + 2x + 3x2 + 4x3

        REQUIRE(1.0 == PolynomialValueAt(polynomial, 0.0));
        REQUIRE(10.0 == PolynomialValueAt(polynomial, 1.0));
        REQUIRE(-2.0 == PolynomialValueAt(polynomial, -1.0));
        REQUIRE(1 + 20 + 300 + 4000 == PolynomialValueAt(polynomial, 10));
        REQUIRE(1 - 20 + 300 - 4000 == PolynomialValueAt(polynomial, -10));
    }
}


TEST_CASE("PolynomialFit - FindRoots", "[Math][PolynomialFit]")
{
    SECTION("Constant")
    {
        std::vector<double> polynomial{ 3.4 };

        std::vector<std::complex<double>> foundRoots;
        bool result = FindRoots(polynomial, foundRoots);

        REQUIRE(true == result);
        REQUIRE(foundRoots.size() == 0);
    }

    SECTION("First order polynomial")
    {
        std::vector<double> polynomial{ 3.4, 1.2 }; // i.e. y = 3.4 + 1.2x

        std::vector<std::complex<double>> foundRoots;
        bool result = FindRoots(polynomial, foundRoots);

        REQUIRE(true == result);
        REQUIRE(foundRoots.size() == 1);
        REQUIRE(foundRoots[0].real() == -polynomial[0] / polynomial[1]);
        REQUIRE(foundRoots[0].imag() == 0.0);
    }

    SECTION("Second order polynomial with imaginary roots")
    {
        std::vector<double> polynomial{ 1.0, 2.0, 4.0 }; // i.e. y = 1 + 2x + 3x2

        std::vector<std::complex<double>> foundRoots;
        bool result = FindRoots(polynomial, foundRoots);

        REQUIRE(true == result);
        REQUIRE(foundRoots.size() == 2);
        REQUIRE(foundRoots[0].real() == -0.25);
        REQUIRE(foundRoots[0].imag() == Approx(0.433013));
        REQUIRE(foundRoots[1].real() == -0.25);
        REQUIRE(foundRoots[1].imag() == -Approx(0.433013));
    }

    SECTION("Second order polynomial with one root")
    {
        std::vector<double> polynomial{ 0.0, 0.0, 3.0 }; // i.e. y = 3x2

        std::vector<std::complex<double>> foundRoots;
        bool result = FindRoots(polynomial, foundRoots);

        REQUIRE(true == result);
        REQUIRE(foundRoots.size() == 2);
        REQUIRE(foundRoots[0].real() == 0.0);
        REQUIRE(foundRoots[0].imag() == 0.0);
        REQUIRE(foundRoots[1].real() == 0.0);
        REQUIRE(foundRoots[1].imag() == 0.0);
    }

    SECTION("Second order polynomial with two roots")
    {
        std::vector<double> polynomial{ 2.0, 0.0, -2.0 }; // i.e. y = 2 - 2x2

        std::vector<std::complex<double>> foundRoots;
        bool result = FindRoots(polynomial, foundRoots);

        REQUIRE(true == result);
        REQUIRE(foundRoots.size() == 2);
        REQUIRE(foundRoots[0].real() == 1.0);
        REQUIRE(foundRoots[0].imag() == 0.0);
        REQUIRE(foundRoots[1].real() == -1.0);
        REQUIRE(foundRoots[1].imag() == 0.0);
    }
}

TEST_CASE("PolynomialFit - First order polynomial", "[Math][PolynomialFit]")
{
    const int order = 1;

    SECTION("No data points returns false. ")
    {
        std::vector<double> xData;
        std::vector<double> yData;
        PolynomialFit sut{ order };

        std::vector<double> resultingPolynomial;
        bool fitSucceeded = sut.FitPolynomial(xData, yData, resultingPolynomial);

        REQUIRE(false == fitSucceeded);
    }

    SECTION("Uneven number of data points in x and y, returns false. ")
    {
        std::vector<double> xData{ 1.0, 5.0, 19.0 };
        std::vector<double> yData{ 1.0, 3.0 };
        PolynomialFit sut{ order };

        std::vector<double> resultingPolynomial;
        bool fitSucceeded = sut.FitPolynomial(xData, yData, resultingPolynomial);

        REQUIRE(false == fitSucceeded);
    }

    SECTION("Horiontal line fits perfectly")
    {
        std::vector<double> actualPolynomial{ 1.0, 0.0 };
        PolynomialFit sut{ order };
        std::vector<double> xData{ 1.0, 7.0 };
        std::vector<double> yData{ 1.0, 1.0 };

        std::vector<double> resultingPolynomial;
        bool fitSucceeded = sut.FitPolynomial(xData, yData, resultingPolynomial);

        REQUIRE(true == fitSucceeded);
        REQUIRE(std::abs(resultingPolynomial[0] - actualPolynomial[0]) < std::numeric_limits<float>::epsilon());
        REQUIRE(std::abs(resultingPolynomial[1] - actualPolynomial[1]) < std::numeric_limits<float>::epsilon());
    }

    SECTION("Positive slope fits perfectly")
    {
        std::vector<double> actualPolynomial{ 1.3, 3.1 };
        PolynomialFit sut{ order };
        std::vector<double> xData{ 1.0, 7.0 };
        std::vector<double> yData;
        for (double x : xData)
        {
            yData.push_back(PolynomialValueAt(actualPolynomial, x));
        }

        std::vector<double> resultingPolynomial;
        bool fitSucceeded = sut.FitPolynomial(xData, yData, resultingPolynomial);

        REQUIRE(true == fitSucceeded);
        REQUIRE(std::abs(resultingPolynomial[0] - actualPolynomial[0]) < std::numeric_limits<float>::epsilon());
        REQUIRE(std::abs(resultingPolynomial[1] - actualPolynomial[1]) < std::numeric_limits<float>::epsilon());
    }

    SECTION("Negative slope fits perfectly")
    {
        std::vector<double> actualPolynomial{ 1.3, -3.1 };
        PolynomialFit sut{ order };
        std::vector<double> xData{ 1.0, 7.0 };
        std::vector<double> yData;
        for (double x : xData)
        {
            yData.push_back(PolynomialValueAt(actualPolynomial, x));
        }

        std::vector<double> resultingPolynomial;
        bool fitSucceeded = sut.FitPolynomial(xData, yData, resultingPolynomial);

        REQUIRE(true == fitSucceeded);
        REQUIRE(std::abs(resultingPolynomial[0] - actualPolynomial[0]) < std::numeric_limits<float>::epsilon());
        REQUIRE(std::abs(resultingPolynomial[1] - actualPolynomial[1]) < std::numeric_limits<float>::epsilon());
    }
}

