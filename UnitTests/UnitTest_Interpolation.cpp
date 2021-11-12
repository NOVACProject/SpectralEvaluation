#include "catch.hpp"
#include <SpectralEvaluation/Interpolation.h>

using namespace novac;

TEST_CASE("Resample")
{
    SECTION("Same input and output - data is unchanged.")
    {
        std::vector<double> x = { 100.0, 101.0, 102.0, 103.0, 104.0, 105.0, 106.0, 107.0, 108.0, 109.0 };
        std::vector<double> y = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 };

        std::vector<double> result;
        novac::Resample(x, y, x, result);

        REQUIRE(result.size() == y.size());
        REQUIRE(result[0] == Approx(y[0]));
        REQUIRE(result[4] == Approx(y[4]));
        REQUIRE(result[9] == Approx(y[9]));
    }

    SECTION("Resample to double resolution - returns expected data.")
    {
        const std::vector<double> x = { -5.0, -4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0 };
        const std::vector<double> y = { 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0 };
        const std::vector<double> newX = 
            { -5.0, -4.5, -4.0, -3.5, -3.0, -2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 
               0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0 };

        std::vector<double> result;
        novac::Resample(x, y, newX, result);

        REQUIRE(result.size() == newX.size());
        REQUIRE(result.front() == Approx(y.front()));
        REQUIRE(result.back() == Approx(y.back()));
        REQUIRE(result[2] == Approx(y[1]));
        REQUIRE(result[8] == Approx(y[4]));
    }
}
