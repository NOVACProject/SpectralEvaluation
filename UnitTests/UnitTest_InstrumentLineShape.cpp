#include "catch.hpp"
#include <string.h>
#include <SpectralEvaluation/Calibration/InstrumentLineShape.h>
#include <SpectralEvaluation/Calibration/InstrumentLineShapeEstimation.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/VectorUtils.h>

// ---------- Unit tests of the free functions found in InstrumentLineShape.h ----------

using namespace novac;

TEST_CASE("SampleInstrumentLineShape Gaussian", "[InstrumentLineShape]")
{
    SECTION("Input has Fwhm of 0.5 nm - result is as expected")
    {
        const double setSigma = 0.5 / 2.35482004;
        const GaussianLineShape gaussian{ setSigma };

        auto result = SampleInstrumentLineShape(gaussian);

        REQUIRE(result.GetSize() > 0);
        REQUIRE(0.5 == Approx(GetFwhm(result)).margin(0.005));
        REQUIRE(1.0 == Approx(Sum(result.m_crossSection)));
    }

    SECTION("Input has negative sigma - same result as for positive sigma")
    {
        const double setSigma = -0.5 / 2.35482004;
        const GaussianLineShape gaussian{ setSigma };

        auto result = SampleInstrumentLineShape(gaussian);

        REQUIRE(result.GetSize() > 0);
        REQUIRE(0.5 == Approx(GetFwhm(result)).margin(0.005));
        REQUIRE(1.0 == Approx(Sum(result.m_crossSection)));
    }
}

TEST_CASE("SampleInstrumentLineShape Super Gaussian", "[InstrumentLineShape]")
{
    SECTION("Input has Fwhm of 0.5nm and power equals two - same result as Gaussian")
    {
        const double fwhm = 0.5;
        const double setWidthParam = 0.5 * fwhm / std::pow(0.69314718056, 1.0 / 2.0);
        const GaussianLineShape gaussian{ fwhm / 2.35482004 };
        auto comparison = SampleInstrumentLineShape(gaussian);
        const SuperGaussianLineShape superGaussian{ setWidthParam, 2.0};

        auto result = SampleInstrumentLineShape(superGaussian);

        REQUIRE(result.GetSize() == comparison.GetSize());
        REQUIRE(0.5 == Approx(GetFwhm(result)).margin(0.005));
        REQUIRE(1.0 == Approx(Sum(result.m_crossSection)));
    }

    SECTION("Input has negative sigma - same result as for positive sigma")
    {
        const double fwhm = 0.5;
        const double setWidthParam = -0.5 * fwhm / std::pow(0.69314718056, 1.0 / 2.0);
        const GaussianLineShape gaussian{ fwhm / 2.35482004 };
        auto comparison = SampleInstrumentLineShape(gaussian);
        const SuperGaussianLineShape superGaussian{ setWidthParam, 2.0 };

        auto result = SampleInstrumentLineShape(superGaussian);

        REQUIRE(result.GetSize() == comparison.GetSize());
        REQUIRE(0.5 == Approx(GetFwhm(result)).margin(0.005));
        REQUIRE(1.0 == Approx(Sum(result.m_crossSection)));
    }
}