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
    // Test cases created witht the help of https://emtoolbox.nist.gov/Wavelength/Edlen.asp
    SECTION("633nm")
    {
        REQUIRE(fabs(1.000271374 - RefractiveIndexOfAir_Edlen(633.0)) < 2e-6);
    }
    SECTION("350nm")
    {
        REQUIRE(fabs(1.000280826 - RefractiveIndexOfAir_Edlen(350.0)) < 1e-5);
    }
    SECTION("400nm")
    {
        REQUIRE(fabs(1.000277516 - RefractiveIndexOfAir_Edlen(400)) < 5e-6);
    }
}

