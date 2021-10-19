#include "catch.hpp"
#include <SpectralEvaluation/Fit/Vector.h>

TEST_CASE("CVector creation from vector data - Basic Operations", "[Fit][CVector]")
{
    const int stepSize = 1;
    const bool takeOwnershipOfData = false;

    SECTION("GetSize returns correct length")
    { 
        std::vector<double> initialValues {9, 8, 7, 6};
        MathFit::CVector sut(initialValues.data(), (int)initialValues.size(), stepSize, takeOwnershipOfData);

        REQUIRE(4 == sut.GetSize());
    }

    SECTION("GetAt returns expected element value")
    {
        std::vector<double> initialValues {9, 8, 7, 6};
        MathFit::CVector sut(initialValues.data(), (int)initialValues.size(), stepSize, takeOwnershipOfData);

        REQUIRE(initialValues[0] == sut.GetAt(0));
        REQUIRE(initialValues[1] == sut.GetAt(1));
        REQUIRE(initialValues[2] == sut.GetAt(2));
        REQUIRE(initialValues[3] == sut.GetAt(3));   
    }

    SECTION("operator[] returns expected element value")
    {
        std::vector<double> initialValues {9, 8, 7, 6};
        MathFit::CVector sut(initialValues.data(), (int)initialValues.size(), stepSize, takeOwnershipOfData);

        REQUIRE(initialValues[0] == sut[0]);
        REQUIRE(initialValues[1] == sut[1]);
        REQUIRE(initialValues[2] == sut[2]);
        REQUIRE(initialValues[3] == sut[3]);
    }

    SECTION("GetSafePtr - returns pointer to original vector")
    {
        std::vector<double> initialValues {9, 8, 7, 6};
        MathFit::CVector sut(initialValues.data(), (int)initialValues.size(), stepSize, takeOwnershipOfData);

        const double* result = sut.GetSafePtr();

        REQUIRE(initialValues.data() == result);
    }

    SECTION("SetAt - updates sut and original vector")
    {
        std::vector<double> initialValues {9, 8, 7, 6};
        MathFit::CVector sut(initialValues.data(), (int)initialValues.size(), stepSize, takeOwnershipOfData);

        sut.SetAt(2, 3.0);

        REQUIRE(3.0 == sut.GetAt(2));
        REQUIRE(3.0 == initialValues[2]);
    }

    SECTION("Zero - fills vector with all zeroes")
    {
        std::vector<double> initialValues {9, 8, 7, 6};
        MathFit::CVector sut(initialValues.data(), (int)initialValues.size(), stepSize, takeOwnershipOfData);

        sut.Zero();

        REQUIRE(0.0 == sut.GetAt(0));
        REQUIRE(0.0 == sut.GetAt(1));
        REQUIRE(0.0 == sut.GetAt(2));
        REQUIRE(0.0 == sut.GetAt(3));
        REQUIRE(sut.IsZero());
    }

    SECTION("Min - returns minimum value")
    {
        std::vector<double> initialValues {9, 8, 7, 6};
        MathFit::CVector sut(initialValues.data(), (int)initialValues.size(), stepSize, takeOwnershipOfData);

        const double result = sut.Min();

        REQUIRE(6.0 == result);
    }

    SECTION("Max - returns maximum value")
    {
        std::vector<double> initialValues {9, 8, 7, 6};
        MathFit::CVector sut(initialValues.data(), (int)initialValues.size(), stepSize, takeOwnershipOfData);

        const double result = sut.Max();

        REQUIRE(9.0 == result);
    }

    SECTION("Max with offset - returns maximum value in selected range")
    {
        std::vector<double> initialValues {9, 8, 7, 6};
        MathFit::CVector sut(initialValues.data(), (int)initialValues.size(), stepSize, takeOwnershipOfData);

        const double result = sut.Max(2);

        REQUIRE(7.0 == result);
    }

}


TEST_CASE("CVector creation from vector data - Polynomial Operations", "[Fit][CVector]")
{
    const int stepSize = 1;
    const bool takeOwnershipOfData = false;

    SECTION("CalcPoly at 1.0 - returns polynomial value at this point")
    {
        const double x = 1.0;
        std::vector<double> initialValues {9, 8, 7, 6}; // interpreted as a polynomial with zeroth order coefficient first
        MathFit::CVector sut(initialValues.data(), (int)initialValues.size(), stepSize, takeOwnershipOfData);
        const double expectedValue = initialValues[0] + x * initialValues[1] + x * x * initialValues[2] + x * x * x * initialValues[3];

        const double result = sut.CalcPoly(x);

        REQUIRE(expectedValue == result);
    }

    SECTION("CalcPoly at 2.0 - returns polynomial value at this point")
    {
        const double x = 2.0;
        std::vector<double> initialValues {9, 8, 7, 6}; // interpreted as a polynomial with zeroth order coefficient first
        MathFit::CVector sut(initialValues.data(), (int)initialValues.size(), stepSize, takeOwnershipOfData);
        const double expectedValue = initialValues[0] + x * initialValues[1] + x * x * initialValues[2] + x * x * x * initialValues[3];

        const double result = sut.CalcPoly(x);

        REQUIRE(expectedValue == result);
    }

    SECTION("CalcPoly at -2.0 - returns polynomial value at this point")
    {
        const double x = -2.0;
        std::vector<double> initialValues {9, 8, 7, 6}; // interpreted as a polynomial with zeroth order coefficient first
        MathFit::CVector sut(initialValues.data(), (int)initialValues.size(), stepSize, takeOwnershipOfData);
        const double expectedValue = initialValues[0] + x * initialValues[1] + x * x * initialValues[2] + x * x * x * initialValues[3];

        const double result = sut.CalcPoly(x);

        REQUIRE(expectedValue == result);
    }


    SECTION("CalcPolySlope at 1.0 - returns derivative of polynomial value at this point")
    {
        const double x = 1.0;
        std::vector<double> initialValues {9, 8, 7, 6}; // interpreted as a polynomial with zeroth order coefficient first
        MathFit::CVector sut(initialValues.data(), (int)initialValues.size(), stepSize, takeOwnershipOfData);
        // calculate the derivative of (9 + 8x + 7x2 + 6x3) at the point x=1
        const double expectedValue = initialValues[1] + 2 * x * initialValues[2] + 3 * x * x * initialValues[3];

        const double result = sut.CalcPolySlope(x);

        REQUIRE(expectedValue == result);
    }

    SECTION("CalcPolySlope at 2.0 - returns polynomial value at this point")
    {
        const double x = 2.0;
        std::vector<double> initialValues {9, 8, 7, 6}; // interpreted as a polynomial with zeroth order coefficient first
        MathFit::CVector sut(initialValues.data(), (int)initialValues.size(), stepSize, takeOwnershipOfData);
        // calculate the derivative of (9 + 8x + 7x2 + 6x3) at the point x=2
        const double expectedValue = initialValues[1] + 2 * x * initialValues[2] + 3 * x * x * initialValues[3];

        const double result = sut.CalcPolySlope(x);

        REQUIRE(expectedValue == result);
    }

    SECTION("CalcPolySlope at -2.0 - returns polynomial value at this point")
    {
        const double x = -2.0;
        std::vector<double> initialValues {9, 8, 7, 6}; // interpreted as a polynomial with zeroth order coefficient first
        MathFit::CVector sut(initialValues.data(), (int)initialValues.size(), stepSize, takeOwnershipOfData);
        // calculate the derivative of (9 + 8x + 7x2 + 6x3) at the point x=-2
        const double expectedValue = initialValues[1] + 2 * x * initialValues[2] + 3 * x * x * initialValues[3];

        const double result = sut.CalcPolySlope(x);

        REQUIRE(expectedValue == result);
    }

}