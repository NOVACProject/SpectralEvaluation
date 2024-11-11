#include "catch.hpp"
#include <SpectralEvaluation/VectorUtils.h>
#include <SpectralEvaluation/Evaluation/BasicMath.h>

TEST_CASE("Min", "[Min][VectorUtils]")
{
    SECTION("Expected value")
    {
        std::vector<double> values{ 2, 3, 4, 5, 6, 7, 8, 9, 1, 8, 7, 6, 5, 4, 3 };

        size_t minIdx = 0U;
        const double result = Min(values, minIdx);

        REQUIRE(minIdx == 8);
        REQUIRE(result == Approx(1.0));
    }

    SECTION("Empty vector returns zero")
    {
        std::vector<double> values;

        size_t minIdx = 1U;
        const double result = Min(values, minIdx);

        REQUIRE(minIdx == 0);
        REQUIRE(result == Approx(0.0));
    }

    SECTION("Max of vector range - returns expected value")
    {
        std::vector<double> values{ -1, 3, 4, 5, 6, 7, 8, 9, 1, 8, 7, 6, 5, 4, -3 };

        // Find the minimum value, excluding the first and the last
        const double result = Min(begin(values) + 1, end(values) - 1);

        // Find the maximum value, excluding the first and the last
        REQUIRE(result == Approx(1.0));
    }
}

TEST_CASE("Max", "[Max][VectorUtils]")
{
    SECTION("Expected value")
    {
        std::vector<double> values{ 2, 3, 4, 5, 6, 7, 8, 9, 1, 8, 7, 6, 5, 4, 3 };

        size_t maxIdx = 0U;
        const double result = Max(values, maxIdx);

        REQUIRE(maxIdx == 7);
        REQUIRE(result == Approx(9.0));
    }

    SECTION("Empty vector returns zero")
    {
        std::vector<double> values;

        size_t maxIdx = 1U;
        const double result = Max(values, maxIdx);

        REQUIRE(maxIdx == 0);
        REQUIRE(result == Approx(0.0));
    }

    SECTION("Max of vector range - returns expected value")
    {
        std::vector<double> values{ 10, 2, 4, 5, 6, 7, 8, 9, 1, 8, 7, 6, 5, 4, 11 };

        // Find the maximum value, excluding the first and the last
        const double result = Max(begin(values) + 1, end(values) - 1);

        REQUIRE(result == Approx(9.0));
    }
}

TEST_CASE("Normalize", "[Normalize][VectorUtils]")
{
    std::vector<double> values{ 2, 3, 4, 5, 6, 7, 8, 9, 8, 7, 6, 5, 4, 3 };
    std::vector<double> result;

    Normalize(values, result);

    SECTION("Result has correct length")
    {
        REQUIRE(result.size() == values.size());
    }

    SECTION("Minimum value is zero")
    {
        REQUIRE(0.0 == result[0]);
        REQUIRE(0.0 == Min(result));
    }

    SECTION("Maximum value is one")
    {
        REQUIRE(1.0 == result[7]);
        REQUIRE(1.0 == Max(result));
    }
}

TEST_CASE("Normalize with empty input", "[Normalize][VectorUtils]")
{
    std::vector<double> input, result;

    Normalize(input, result);

    REQUIRE(0 == result.size());
}

TEST_CASE("NormalizeArea", "[NormalizeArea][VectorUtils]")
{
    std::vector<double> values{ 2, 3, 4, 5, 6, 7, 8, 9, 8, 7, 6, 5, 4, 3 };
    std::vector<double> result;

    NormalizeArea(values, result);

    SECTION("Result has correct length")
    {
        REQUIRE(result.size() == values.size());
    }

    SECTION("Minimum value is zero")
    {
        REQUIRE(0.0 == result[0]);
        REQUIRE(0.0 == Min(result));
    }

    SECTION("Sum of all values is one")
    {
        REQUIRE(1.0 == Sum(result));
    }
}

TEST_CASE("NormalizeArea with empty input", "[Normalize][VectorUtils]")
{
    std::vector<double> input, result;

    NormalizeArea(input, result);

    REQUIRE(0 == result.size());
}


TEST_CASE("Reverse", "[Reverse][VectorUtils]")
{
    SECTION("Expected contents for odd length vector")
    {
        // There are 9 values here, with the '5' being in the center
        std::vector<double> values{ 1, 2, 3, 4, 5, 6, 7, 8, 9 };

        Reverse(values);

        REQUIRE(9 == values.size());
        REQUIRE(9 == Approx(values[0]));
        REQUIRE(8 == Approx(values[1]));
        REQUIRE(7 == Approx(values[2]));
        REQUIRE(6 == Approx(values[3]));
        REQUIRE(5 == Approx(values[4]));
        REQUIRE(4 == Approx(values[5]));
        REQUIRE(3 == Approx(values[6]));
        REQUIRE(2 == Approx(values[7]));
        REQUIRE(1 == Approx(values[8]));
    }

    SECTION("Expected contents for even length vector")
    {
        // There are 8 values here, and hence no value is in the center
        std::vector<double> values{ 1, 2, 3, 4, 5, 6, 7, 8 };

        Reverse(values);

        REQUIRE(8 == values.size());
        REQUIRE(8 == Approx(values[0]));
        REQUIRE(7 == Approx(values[1]));
        REQUIRE(6 == Approx(values[2]));
        REQUIRE(5 == Approx(values[3]));
        REQUIRE(4 == Approx(values[4]));
        REQUIRE(3 == Approx(values[5]));
        REQUIRE(2 == Approx(values[6]));
        REQUIRE(1 == Approx(values[7]));
    }

    SECTION("Empty vector returns")
    {
        std::vector<double> values;

        Reverse(values);

        REQUIRE(0 == values.size());
    }
}

TEST_CASE("FindValue - Constantly increasing vector", "[FindValue][VectorUtils]")
{
    const std::vector<double> values{ 1, 2, 3, 4 };

    SECTION("Finds values at integer index points.")
    {
        REQUIRE(0.0 == Approx(FindValue(values, 1.0, 0, 4)));
        REQUIRE(2.0 == Approx(FindValue(values, 3.0, 0, 4)));
        REQUIRE(3.0 == Approx(FindValue(values, 4.0, 0, 4)));
    }

    SECTION("Finds values midway between points.")
    {
        REQUIRE(0.5 == Approx(FindValue(values, 1.5, 0, 4)));
        REQUIRE(2.5 == Approx(FindValue(values, 3.5, 0, 4)));
    }

    SECTION("Finds values quarter between points.")
    {
        REQUIRE(0.25 == Approx(FindValue(values, 1.25, 0, 4)));
        REQUIRE(2.25 == Approx(FindValue(values, 3.25, 0, 4)));
    }

    SECTION("Finds values three-quarters between points.")
    {
        REQUIRE(0.75 == Approx(FindValue(values, 1.75, 0, 4)));
        REQUIRE(2.75 == Approx(FindValue(values, 3.75, 0, 4)));
    }
}

TEST_CASE("FindValue - Constantly decreasing vector", "[FindValue][VectorUtils]")
{
    const std::vector<double> values{ 4, 3, 2, 1 };

    SECTION("Finds values at integer index points.")
    {
        REQUIRE(0.0 == Approx(FindValue(values, 4.0, 0, 4)));
        REQUIRE(2.0 == Approx(FindValue(values, 2.0, 0, 4)));
        REQUIRE(3.0 == Approx(FindValue(values, 1.0, 0, 4)));
    }

    SECTION("Finds values midway between points.")
    {
        REQUIRE(2.5 == Approx(FindValue(values, 1.5, 0, 4)));
        REQUIRE(0.5 == Approx(FindValue(values, 3.5, 0, 4)));
    }

    SECTION("Finds values quarter between points.")
    {
        REQUIRE(2.25 == Approx(FindValue(values, 1.75, 0, 4)));
        REQUIRE(0.25 == Approx(FindValue(values, 3.75, 0, 4)));
    }

    SECTION("Finds values three-quarters between points.")
    {
        REQUIRE(2.75 == Approx(FindValue(values, 1.25, 0, 4)));
        REQUIRE(0.75 == Approx(FindValue(values, 3.25, 0, 4)));
    }
}

TEST_CASE("FindValue - invalid ranges", "[FindValue][VectorUtils]")
{
    const std::vector<double> values{ 4, 3, 2, 1 };

    SECTION("Startidx above StopIdx.")
    {
        REQUIRE(-1.0 == Approx(FindValue(values, 4.0, 3, 1)));
    }

    SECTION("Startidx out of bounds.")
    {
        REQUIRE(-1.0 == Approx(FindValue(values, 4.0, 6, 9)));
    }
}


TEST_CASE("GetAt", "[GetAt][VectorUtils]")
{
    const std::vector<double> values{ 1, 2, 3, 4 };

    SECTION("Finds values at integer index points.")
    {
        REQUIRE(1.0 == Approx(GetAt(values, 0.0)));
        REQUIRE(3.0 == Approx(GetAt(values, 2.0)));
        REQUIRE(4.0 == Approx(GetAt(values, 3.0)));
    }

    SECTION("Finds values at quarter index points.")
    {
        REQUIRE(1.25 == Approx(GetAt(values, 0.25)));
        REQUIRE(3.25 == Approx(GetAt(values, 2.25)));
    }
}

TEST_CASE("FindNLowest", "[FindNLowest][VectorUtils]")
{
    std::vector<double> input;
    std::vector<double> result;

    SECTION("Empty input, returns empty result vector")
    {
        result.resize(10);
        FindNLowest(input, 1, result);
        REQUIRE(result.size() == 0);
    }

    SECTION("Typical case - returns n lowest values (n = 3)")
    {
        input = { 3, 7, 2, 9, 4, 12, 3, 7 };
        FindNLowest(input, 3, result);
        REQUIRE(3 == result.size());
        REQUIRE(2 == result[0]);
        REQUIRE(3 == result[1]);
        REQUIRE(3 == result[2]);
    }

    SECTION("Typical case - returns n lowest values (n = 5)")
    {
        input = { 3, 7, 2, 9, 4, 12, 3, 7 };
        FindNLowest(input, 5, result);
        REQUIRE(5 == result.size());
        REQUIRE(2 == result[0]);
        REQUIRE(3 == result[1]);
        REQUIRE(3 == result[2]);
        REQUIRE(4 == result[3]);
        REQUIRE(7 == result[4]);
    }

    SECTION("N > input.size(), returns all values sorted")
    {
        input = { 3, 7, 2, 9, 4, 12, 3, 7 };
        FindNLowest(input, 12, result);
        REQUIRE(8 == result.size());
        REQUIRE(2 == result[0]);
        REQUIRE(3 == result[1]);
        REQUIRE(3 == result[2]);
        REQUIRE(4 == result[3]);
        REQUIRE(7 == result[4]);
        REQUIRE(7 == result[5]);
        REQUIRE(9 == result[6]);
        REQUIRE(12 == result[7]);
    }
}

TEST_CASE("HighPassBinomial - Constant vector", "[HighPassBinomial][VectorUtils]")
{
    CBasicMath math;

    SECTION("Constant value vector zeroes out after one iteration")
    {
        std::vector<double> values(128, 5.0);
        std::vector<double> values2(128, 5.0);

        // Act
        ::HighPassBinomial(values, 1);
        math.HighPassBinomial(values2.data(), static_cast<int>(values2.size()), 1);

        // Assert, the two vectors should now have the same value.
        for (size_t ii = 0; ii < values.size(); ++ii)
        {
            REQUIRE(Approx(values[ii]) == values2[ii]);
        }
    }

    SECTION("Constant value vector zeroes out after 500 iterations")
    {
        std::vector<double> values(128, 5.0);
        std::vector<double> values2(128, 5.0);

        // Act
        ::HighPassBinomial(values, 500);
        math.HighPassBinomial(values2.data(), static_cast<int>(values2.size()), 500);

        // Assert, the two vectors should now have the same value.
        for (size_t ii = 0; ii < values.size(); ++ii)
        {
            REQUIRE(Approx(values[ii]) == values2[ii]);
        }
    }

    SECTION("Odd sized Constant value vector zeroes out after 500 iterations")
    {
        std::vector<double> values(61, 5.0);
        std::vector<double> values2(61, 5.0);

        // Act
        ::HighPassBinomial(values, 1);
        math.HighPassBinomial(values2.data(), static_cast<int>(values2.size()), 500);

        // Assert, the two vectors should now have the same value.
        for (size_t ii = 0; ii < values.size(); ++ii)
        {
            REQUIRE(Approx(values[ii]) == values2[ii]);
        }
    }
}


TEST_CASE("HighPassBinomial - Random vector matches result from BasicMath", "[HighPassBinomial][VectorUtils]")
{
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<> dist(-10.0, 10.0); // distribution in range [0, +10.0]

    std::vector<double> values;
    std::vector<double> values2;

    CBasicMath math;

    SECTION("Even sized value vector: one iteration")
    {
        for (size_t ii = 0; ii < 128; ++ii)
        {
            const double v = dist(rng);
            values.push_back(v);
            values2.push_back(v);
        }

        // Act
        ::HighPassBinomial(values, 1);
        math.HighPassBinomial(values2.data(), static_cast<int>(values2.size()), 1);

        // Assert, the two vectors should now have the same value.
        for (size_t ii = 0; ii < values.size(); ++ii)
        {
            REQUIRE(Approx(values[ii]) == values2[ii]);
        }
    }
    SECTION("Even sized value vector: 500 iterations")
    {
        for (size_t ii = 0; ii < 128; ++ii)
        {
            const double v = dist(rng);
            values.push_back(v);
            values2.push_back(v);
        }

        // Act
        ::HighPassBinomial(values, 500);
        math.HighPassBinomial(values2.data(), static_cast<int>(values2.size()), 500);

        // Assert, the two vectors should now have the same value.
        for (size_t ii = 0; ii < values.size(); ++ii)
        {
            REQUIRE(Approx(values[ii]) == values2[ii]);
        }
    }

    SECTION("Odd sized value vector: 500 iterations")
    {
        for (size_t ii = 0; ii < 61; ++ii)
        {
            const double v = dist(rng);
            values.push_back(v);
            values2.push_back(v);
        }

        // Act
        ::HighPassBinomial(values, 500);
        math.HighPassBinomial(values2.data(), static_cast<int>(values2.size()), 500);

        // Assert, the two vectors should now have the same value.
        for (size_t ii = 0; ii < values.size(); ++ii)
        {
            REQUIRE(Approx(values[ii]) == values2[ii]);
        }
    }
}


TEST_CASE("LowPassBinomial - Constant vector", "[LowPassBinomial][VectorUtils]")
{
    CBasicMath math;

    SECTION("Constant value vector zeroes out after one iteration")
    {
        std::vector<double> values(128, 5.0);
        std::vector<double> values2(128, 5.0);

        // Act
        ::LowPassBinomial(values, 1);
        math.LowPassBinomial(values2.data(), static_cast<int>(values2.size()), 1);

        // Assert, the two vectors should now have the same value.
        for (size_t ii = 0; ii < values.size(); ++ii)
        {
            REQUIRE(Approx(values[ii]) == values2[ii]);
        }
    }

    SECTION("Constant value vector zeroes out after 500 iterations")
    {
        std::vector<double> values(128, 5.0);
        std::vector<double> values2(128, 5.0);

        // Act
        ::LowPassBinomial(values, 500);
        math.LowPassBinomial(values2.data(), static_cast<int>(values2.size()), 500);

        // Assert, the two vectors should now have the same value.
        for (size_t ii = 0; ii < values.size(); ++ii)
        {
            REQUIRE(Approx(values[ii]) == values2[ii]);
        }
    }

    SECTION("Odd sized Constant value vector zeroes out after 500 iterations")
    {
        std::vector<double> values(61, 5.0);
        std::vector<double> values2(61, 5.0);

        // Act
        ::LowPassBinomial(values, 1);
        math.LowPassBinomial(values2.data(), static_cast<int>(values2.size()), 500);

        // Assert, the two vectors should now have the same value.
        for (size_t ii = 0; ii < values.size(); ++ii)
        {
            REQUIRE(Approx(values[ii]) == values2[ii]);
        }
    }
}


TEST_CASE("LowPassBinomial - Random vector matches result from BasicMath", "[LowPassBinomial][VectorUtils]")
{
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_real_distribution<> dist(-10.0, 10.0); // distribution in range [0, +10.0]

    std::vector<double> values;
    std::vector<double> values2;

    CBasicMath math;

    SECTION("Even sized value vector: one iteration")
    {
        for (size_t ii = 0; ii < 128; ++ii)
        {
            const double v = dist(rng);
            values.push_back(v);
            values2.push_back(v);
        }

        // Act
        ::LowPassBinomial(values, 1);
        math.LowPassBinomial(values2.data(), static_cast<int>(values2.size()), 1);

        // Assert, the two vectors should now have the same value.
        for (size_t ii = 0; ii < values.size(); ++ii)
        {
            REQUIRE(Approx(values[ii]) == values2[ii]);
        }
    }

    SECTION("Even sized value vector: 500 iterations")
    {
        for (size_t ii = 0; ii < 128; ++ii)
        {
            const double v = dist(rng);
            values.push_back(v);
            values2.push_back(v);
        }

        // Act
        ::LowPassBinomial(values, 500);
        math.LowPassBinomial(values2.data(), static_cast<int>(values2.size()), 500);

        // Assert, the two vectors should now have the same value.
        for (size_t ii = 0; ii < values.size(); ++ii)
        {
            REQUIRE(Approx(values[ii]) == values2[ii]);
        }
    }

    SECTION("Odd sized value vector: 500 iterations")
    {
        for (size_t ii = 0; ii < 61; ++ii)
        {
            const double v = dist(rng);
            values.push_back(v);
            values2.push_back(v);
        }

        // Act
        ::LowPassBinomial(values, 500);
        math.LowPassBinomial(values2.data(), static_cast<int>(values2.size()), 500);

        // Assert, the two vectors should now have the same value.
        for (size_t ii = 0; ii < values.size(); ++ii)
        {
            REQUIRE(Approx(values[ii]) == values2[ii]);
        }
    }
}
