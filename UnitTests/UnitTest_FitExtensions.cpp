#include "catch.hpp"
#include <SpectralEvaluation/FitExtensions/SuperGaussFunction.h>
#include <SpectralEvaluation/Fit/GaussFunction.h>

// -------- SuperGaussFunction --------
TEST_CASE("SuperGaussFunction : Power equals two, idential to Gaussian Function", "[FitExtensions]")
{
    const double gaussianSigma = 2.3;
    
    MathFit::CGaussFunction gauss;
    gauss.SetCenter(5.5);
    gauss.SetSigma(gaussianSigma);

    MathFit::CSuperGaussFunction sut;
    sut.SetCenter(5.5);
    sut.SetW(gaussianSigma * std::sqrt(2.0)); // there's a difference in the interpretation of the parameters.
    sut.SetK(2.0);

    SECTION("GetValue returns same")
    {
        REQUIRE(std::abs(gauss.GetValue(0.0) - sut.GetValue(0.0)) < 1e-4);
        REQUIRE(std::abs(gauss.GetValue(5.0) - sut.GetValue(5.0)) < 1e-4);
        REQUIRE(std::abs(gauss.GetValue(14.0) - sut.GetValue(14.0)) < 1e-4);
    }

    SECTION("GetValues returns same")
    {
        const double xMin = -15.0;
        const double xDelta = 1.5;
        double curX = xMin;
        std::vector<double> xValues(15, 0);
        std::generate_n(begin(xValues), 15, [&] { curX += xDelta; return curX; });
        MathFit::CVector xVector(xValues.data(), (int)xValues.size(), 1, false);
        MathFit::CVector gaussValues(xVector.GetSize());
        gauss.GetValues(xVector, gaussValues);

        MathFit::CVector sutValues(xVector.GetSize());
        sut.GetValues(xVector, sutValues);

        for (int ii = 0; ii < sutValues.GetSize(); ++ii)
        {
            REQUIRE(std::abs(gaussValues.GetAt(ii) - sutValues.GetAt(ii)) < 1e-4);
        }
    }

    SECTION("GetSlope returns same")
    {
        double gaussValue = gauss.GetSlope(0.0);
        double sutValue = sut.GetSlope(0.0);
        REQUIRE(std::abs(gaussValue - sutValue) < 1e-4);

        REQUIRE(std::abs(gauss.GetSlope(5.0) - sut.GetSlope(5.0)) < 1e-4);
        REQUIRE(std::abs(gauss.GetSlope(14.0) - sut.GetSlope(14.0)) < 1e-4);
    }

    SECTION("GetSlopes returns same")
    {
        const double xMin = -15.0;
        const double xDelta = 1.5;
        double curX = xMin;
        std::vector<double> xValues(15, 0);
        std::generate_n(begin(xValues), 15, [&] { curX += xDelta; return curX; });
        MathFit::CVector xVector(xValues.data(), (int)xValues.size(), 1, false);
        MathFit::CVector gaussValues(xVector.GetSize());
        gauss.GetSlopes(xVector, gaussValues);

        MathFit::CVector sutValues(xVector.GetSize());
        sut.GetSlopes(xVector, sutValues);

        for (int ii = 0; ii < sutValues.GetSize(); ++ii)
        {
            REQUIRE(std::abs(gaussValues.GetAt(ii) - sutValues.GetAt(ii)) < 1e-4);
        }
    }
}