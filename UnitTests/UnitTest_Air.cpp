#include "catch.hpp"
#include <SpectralEvaluation/Air.h>

TEST_CASE("NmVacuumToNmAir")
{
    SECTION("Red HeNe returns correct value.")
    {
        REQUIRE(Approx(632.815) == NmVacuumToNmAir(632.9915));
    }
}