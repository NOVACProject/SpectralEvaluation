#include "catch.hpp"
#include <string.h>
#include <SpectralEvaluation/Spectra/SpectrometerModel.h>

std::vector<double> CreateGaussian(double center, double sigma, const std::vector<double>& x); // located elsewhere

namespace novac
{
    TEST_CASE("GuessModelFromSerial correct result.", "[CSpectrometerDatabase][SpectrometerModel][GuessModelFromSerial]")
    {
        SECTION("D2J1840 is S2000")
        {
            const auto result = CSpectrometerDatabase::GuessModelFromSerial("D2J1840");

            REQUIRE("S2000" == result.modelName);
        }

        SECTION("I2J4969 is S2000")
        {
            const auto result = CSpectrometerDatabase::GuessModelFromSerial("I2J4969");

            REQUIRE("S2000" == result.modelName);
        }

        SECTION("I2P0093 is S2000")
        {
            const auto result = CSpectrometerDatabase::GuessModelFromSerial("I2P0093");

            REQUIRE("S2000" == result.modelName);
        }

        SECTION("1904156M1 is AVASPEC")
        {
            const auto result = CSpectrometerDatabase::GuessModelFromSerial("1904156M1");

            REQUIRE("AVASPEC" == result.modelName);
        }

        SECTION("FLMS14634 is FLAME")
        {
            const auto result = CSpectrometerDatabase::GuessModelFromSerial("FLMS14634");

            REQUIRE("FLAME" == result.modelName);
        }

        // SECTION("H05705 is FLAME")
        // {
        //     const auto result = CSpectrometerDatabase::GuessModelFromSerial("H05705");
        // 
        //     REQUIRE("FLAME" == result.modelName);
        // }

        SECTION("HR2B2203 is HR2000")
        {
            const auto result = CSpectrometerDatabase::GuessModelFromSerial("HR2B2203");

            REQUIRE("HR2000" == result.modelName);
        }

        SECTION("MAYP11440 is MAYAPRO")
        {
            const auto result = CSpectrometerDatabase::GuessModelFromSerial("MAYP11440");

            REQUIRE("MAYAPRO" == result.modelName);
        }

        SECTION("USB+U12835 is USB2000+")
        {
            const auto result = CSpectrometerDatabase::GuessModelFromSerial("USB+U12835");

            REQUIRE("USB2000+" == result.modelName);
        }

        SECTION("USB2+F03188 is USB2000+")
        {
            const auto result = CSpectrometerDatabase::GuessModelFromSerial("USB2+F03188");

            REQUIRE("USB2000+" == result.modelName);
        }

        SECTION("USB2+U12835 is USB2000+")
        {
            const auto result = CSpectrometerDatabase::GuessModelFromSerial("USB2+U12835");

            REQUIRE("USB2000+" == result.modelName);
        }

        SECTION("USB2G789 is USB2000")
        {
            const auto result = CSpectrometerDatabase::GuessModelFromSerial("USB2G789");

            REQUIRE("USB2000" == result.modelName);
        }

        SECTION("NOT SET is UNKNOWN")
        {
            const auto result = CSpectrometerDatabase::GuessModelFromSerial("NOT SET");

            REQUIRE("UNKNOWN" == result.modelName);
        }

        // NOT SET
    }
}

