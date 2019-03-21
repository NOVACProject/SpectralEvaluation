#include "catch.hpp"
#include <SpectralEvaluation/Air.h>

TEST_CASE("NmVacuumToNmAir")
{
    SECTION("Red HeNe returns correct value.")
    {
        REQUIRE(Approx(632.819) == NmVacuumToNmAir(632.9915));
    }
}

TEST_CASE("RefractiveIndexOfAir_Edlen")
{
    SECTION("633nm")
    {
        REQUIRE(fabs(1.000271 - RefractiveIndexOfAir_Edlen(633.0)) < 2e-6);
    }
}