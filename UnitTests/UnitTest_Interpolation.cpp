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

TEST_CASE("GetFractionalIndex", "[Interpolation]")
{
    SECTION("Empty vector, returns Nan.")
    {
        std::vector<double> haystack;
        const double needle = 0.0;

        double result = GetFractionalIndex(haystack, needle);

        REQUIRE(std::isnan(result));
    }

    SECTION("Value equal to first - returns zero.")
    {
        std::vector<double> haystack{ 0.0, 1.0, 2.0 };
        const double needle = 0.0;

        double result = GetFractionalIndex(haystack, needle);

        REQUIRE(result == Approx(0.0));
    }

    SECTION("Value mid point between first and second - returns zero point five.")
    {
        std::vector<double> haystack{ 0.0, 10.0, 20.0 };
        const double needle = 5.0;

        double result = GetFractionalIndex(haystack, needle);

        REQUIRE(result == Approx(0.5));
    }

    SECTION("Value equal to last - returns last index.")
    {
        std::vector<double> haystack{ 0.0, 10.0, 20.0 };
        const double needle = 20.0;

        double result = GetFractionalIndex(haystack, needle);

        REQUIRE(result == Approx(2.0));
    }

    SECTION("Value smaller than first - returns NaN.")
    {
        std::vector<double> haystack{ 0.0, 10.0, 20.0 };
        const double needle = -5.0;

        double result = GetFractionalIndex(haystack, needle);

        REQUIRE(std::isnan(result));
    }

    SECTION("Value larger than last - returns NaN.")
    {
        std::vector<double> haystack{ 0.0, 10.0, 20.0 };
        const double needle = 30.0;

        double result = GetFractionalIndex(haystack, needle);

        REQUIRE(std::isnan(result));
    }
}