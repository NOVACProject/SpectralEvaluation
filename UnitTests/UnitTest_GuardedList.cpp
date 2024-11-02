#include "catch.hpp"
#include <SpectralEvaluation/ThreadUtils.h>

namespace novac
{

TEST_CASE("GuardedList initial state is as expected", "[GuardedList]")
{
    GuardedList<double> sut;

    SECTION("Initially empty")
    {
        REQUIRE(0 == sut.Size());
    }

    SECTION("PopFront from empty list return false and no item")
    {
        double value = 999.0;

        // Act
        const bool result = sut.PopFront(value);

        // Assert
        REQUIRE(false == result);
        REQUIRE(Approx(value) == 999.0); // value not modified
    }

    SECTION("Copy to vector, clears destination vector")
    {
        std::vector<double> destination(5, 3.14);

        // Act
        sut.CopyTo(destination);

        // Assert
        REQUIRE(destination.size() == 0);
    }
}

TEST_CASE("GuardedList with some added values behaves as expected", "[GuardedList]")
{
    GuardedList<double> sut;
    sut.AddItem(2.0);
    sut.AddItem(4.0);
    sut.AddItem(8.0);

    SECTION("Correct size")
    {
        REQUIRE(3 == sut.Size());
    }

    SECTION("Clear removes all items")
    {
        // Act
        sut.Clear();

        // Assert
        REQUIRE(0 == sut.Size());
    }

    SECTION("PopFront return true and gives first item")
    {
        double value = 999.0;

        // Act
        const bool result = sut.PopFront(value);

        // Assert
        REQUIRE(true == result);
        REQUIRE(Approx(value) == 2.0); // first value added

        REQUIRE(sut.Size() == 2); // first element removed
    }

    SECTION("PopFront multiple times removes one item at a time until no more")
    {
        double value1 = 999.0;
        double value2 = 999.0;
        double value3 = 999.0;
        double value4 = 999.0;

        // Act
        const bool result1 = sut.PopFront(value1);
        const bool result2 = sut.PopFront(value2);
        const bool result3 = sut.PopFront(value3);
        const bool result4 = sut.PopFront(value4);

        // Assert
        REQUIRE(true == result1);
        REQUIRE(Approx(value1) == 2.0); // first value added

        REQUIRE(true == result2);
        REQUIRE(Approx(value2) == 4.0); // second value added

        REQUIRE(true == result3);
        REQUIRE(Approx(value3) == 8.0); // third value added

        REQUIRE(false == result4);
        REQUIRE(Approx(value4) == 999.0); // not modified

        REQUIRE(sut.Size() == 0); // all elements removed
    }

    SECTION("Copy to vector gives expected destination vector")
    {
        std::vector<double> destination(5, 3.14);

        // Act
        sut.CopyTo(destination);

        // Assert
        REQUIRE(destination.size() == 3);
        REQUIRE(Approx(destination[0]) == 2.0);
        REQUIRE(Approx(destination[1]) == 4.0);
        REQUIRE(Approx(destination[2]) == 8.0);

        REQUIRE(sut.Size() == 3); // not changed.
    }
}


} // namespace novac
