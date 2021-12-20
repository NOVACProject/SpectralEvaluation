#include <SpectralEvaluation/Calibration/WavelengthCalibration.h>
#include "catch.hpp"

namespace novac
{
    static std::string GetTestDataDirectory()
    {
#ifdef _MSC_VER
        return std::string("../TestData/");
#else
        return std::string("TestData/");
#endif // _MSC_VER 
    }

    TEST_CASE("GetPixelToWavelengthMappingFromFile Clb file with only one column.",
    "[GetPixelToWavelengthMappingFromFile][WavelengthCalibration][IntegrationTest]")
    {
        const std::string filename = GetTestDataDirectory() + "D2J2200/D2J2200_Master.clb";

        const auto result = GetPixelToWavelengthMappingFromFile(filename);

        REQUIRE(2048 == result.size());
        REQUIRE(278.4631392 == Approx(result.front()));
        REQUIRE(425.2065349 == Approx(result.back()));
    }

    TEST_CASE("GetPixelToWavelengthMappingFromFile Txt file with two columns.",
    "[GetPixelToWavelengthMappingFromFile][WavelengthCalibration][IntegrationTest]")
    {
        const std::string filename = GetTestDataDirectory() + "D2J2200/Hg_D2J2200_all.Master.Sample.txt";

        const auto result = GetPixelToWavelengthMappingFromFile(filename);

        REQUIRE(2048 == result.size());
        REQUIRE(278.4631391 == Approx(result.front()));
        REQUIRE(425.2065349 == Approx(result.back()));
    }

}

