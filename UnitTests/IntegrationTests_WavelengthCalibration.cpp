#include <SpectralEvaluation/Calibration/WavelengthCalibration.h>
#include "catch.hpp"
#include "TestData.h"

namespace novac
{
    TEST_CASE("GetPixelToWavelengthMappingFromFile Clb file with only one column.",
    "[GetPixelToWavelengthMappingFromFile][WavelengthCalibration][IntegrationTest]")
    {
        const std::string filename = TestData::GetInitialPixelToWavelengthCalibration_D2J2200();

        const auto result = GetPixelToWavelengthMappingFromFile(filename);

        REQUIRE(2048 == result.size());
        REQUIRE(278.4631392 == Approx(result.front()));
        REQUIRE(425.2065349 == Approx(result.back()));
    }

    TEST_CASE("GetPixelToWavelengthMappingFromFile Txt file with two columns.",
    "[GetPixelToWavelengthMappingFromFile][WavelengthCalibration][IntegrationTest]")
    {
        const std::string filename = TestData::GetMeasuredMercurySpectrum_D2J2200();

        const auto result = GetPixelToWavelengthMappingFromFile(filename);

        REQUIRE(2048 == result.size());
        REQUIRE(278.4631391 == Approx(result.front()));
        REQUIRE(425.2065349 == Approx(result.back()));
    }

}

