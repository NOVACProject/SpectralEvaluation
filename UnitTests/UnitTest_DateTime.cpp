#include "catch.hpp"
#include <SpectralEvaluation/DateTime.h>
#include <sstream>

namespace novac
{

    TEST_CASE("CDateTime - Default constructor - all values set to zero", "[DateTime]")
    {
        CDateTime sut;

        REQUIRE(0 == sut.year);
        REQUIRE(0 == sut.month);
        REQUIRE(0 == sut.day);
        REQUIRE(0 == sut.hour);
        REQUIRE(0 == sut.minute);
        REQUIRE(0 == sut.second);
        REQUIRE(0 == sut.millisecond);
    }

    TEST_CASE("CDateTime - Copy constructor copy copies all values", "[DateTime]")
    {
        CDateTime sut;
        sut.year = 2009;
        sut.month = 12;
        sut.day = 14;
        sut.hour = 16;
        sut.minute = 34;
        sut.second = 23;
        sut.millisecond = 32;

        CDateTime copy = CDateTime{ sut };

        REQUIRE(copy.year == sut.year);
        REQUIRE(copy.month == sut.month);
        REQUIRE(copy.day == sut.day);
        REQUIRE(copy.hour == sut.hour);
        REQUIRE(copy.minute == sut.minute);
        REQUIRE(copy.second == sut.second);
        REQUIRE(copy.millisecond == sut.millisecond);
    }

    TEST_CASE("CDateTime - Copy assignment operator copies all values", "[DateTime]")
    {
        CDateTime sut;
        sut.year = 2009;
        sut.month = 12;
        sut.day = 14;
        sut.hour = 16;
        sut.minute = 34;
        sut.second = 23;
        sut.millisecond = 32;
        CDateTime copy;

        copy = sut;

        REQUIRE(copy.year == sut.year);
        REQUIRE(copy.month == sut.month);
        REQUIRE(copy.day == sut.day);
        REQUIRE(copy.hour == sut.hour);
        REQUIRE(copy.minute == sut.minute);
        REQUIRE(copy.second == sut.second);
        REQUIRE(copy.millisecond == sut.millisecond);
    }

    TEST_CASE("CDateTime - Constructor with parameters, sets correct values.", "[DateTime]")
    {
        CDateTime sut{ 2010, 11, 3, 22, 50, 56 };

        REQUIRE(2010 == sut.year);
        REQUIRE(11 == sut.month);
        REQUIRE(3 == sut.day);
        REQUIRE(22 == sut.hour);
        REQUIRE(50 == sut.minute);
        REQUIRE(56 == sut.second);
        REQUIRE(0 == sut.millisecond);
    }

    TEST_CASE("CDateTime - Increment positive number of seconds.", "[DateTime]")
    {
        SECTION("Common time")
        {
            CDateTime sut{ 2010, 11, 3, 22, 50, 12 };
            sut.Increment(5);

            REQUIRE(2010 == sut.year);
            REQUIRE(11 == sut.month);
            REQUIRE(3 == sut.day);
            REQUIRE(22 == sut.hour);
            REQUIRE(50 == sut.minute);
            REQUIRE(17 == sut.second);
            REQUIRE(0 == sut.millisecond);
        }

        SECTION("Minute wrap around")
        {
            CDateTime sut{ 2010, 11, 3, 22, 50, 56 };
            sut.Increment(5);

            REQUIRE(2010 == sut.year);
            REQUIRE(11 == sut.month);
            REQUIRE(3 == sut.day);
            REQUIRE(22 == sut.hour);
            REQUIRE(51 == sut.minute);
            REQUIRE(1 == sut.second);
            REQUIRE(0 == sut.millisecond);
        }

        SECTION("Hour wrap around")
        {
            CDateTime sut{ 2010, 11, 3, 22, 59, 57 };
            sut.Increment(5);

            REQUIRE(2010 == sut.year);
            REQUIRE(11 == sut.month);
            REQUIRE(3 == sut.day);
            REQUIRE(23 == sut.hour);
            REQUIRE(0 == sut.minute);
            REQUIRE(2 == sut.second);
            REQUIRE(0 == sut.millisecond);
        }

        SECTION("Day wrap around")
        {
            CDateTime sut{ 2010, 11, 3, 23, 59, 57 };
            sut.Increment(5);

            REQUIRE(2010 == sut.year);
            REQUIRE(11 == sut.month);
            REQUIRE(4 == sut.day);
            REQUIRE(0 == sut.hour);
            REQUIRE(0 == sut.minute);
            REQUIRE(2 == sut.second);
            REQUIRE(0 == sut.millisecond);
        }

        SECTION("Month wrap around")
        {
            CDateTime sut{ 2010, 12, 31, 23, 59, 57 };
            sut.Increment(10);

            REQUIRE(2011 == sut.year);
            REQUIRE(1 == sut.month);
            REQUIRE(1 == sut.day);
            REQUIRE(0 == sut.hour);
            REQUIRE(0 == sut.minute);
            REQUIRE(7 == sut.second);
            REQUIRE(0 == sut.millisecond);
        }

        SECTION("One day - Month wrap around")
        {
            CDateTime sut{ 2010, 12, 31, 23, 59, 57 };
            sut.Increment(23 * 3600); // 23 hours

            REQUIRE(2011 == sut.year);
            REQUIRE(1 == sut.month);
            REQUIRE(1 == sut.day);
            REQUIRE(22 == sut.hour);
            REQUIRE(59 == sut.minute);
            REQUIRE(57 == sut.second);
            REQUIRE(0 == sut.millisecond);
        }
    }

    TEST_CASE("CDateTime - Increment with negative number of seconds.", "[DateTime]")
    {
        SECTION("Common time")
        {
            CDateTime sut{ 2010, 11, 3, 22, 50, 12 };
            sut.Increment(-5);

            REQUIRE(2010 == sut.year);
            REQUIRE(11 == sut.month);
            REQUIRE(3 == sut.day);
            REQUIRE(22 == sut.hour);
            REQUIRE(50 == sut.minute);
            REQUIRE(7 == sut.second);
            REQUIRE(0 == sut.millisecond);
        }

        SECTION("Minute wrap around")
        {
            CDateTime sut{ 2010, 11, 3, 22, 50, 04 };
            sut.Increment(-5);

            REQUIRE(2010 == sut.year);
            REQUIRE(11 == sut.month);
            REQUIRE(3 == sut.day);
            REQUIRE(22 == sut.hour);
            REQUIRE(49 == sut.minute);
            REQUIRE(59 == sut.second);
            REQUIRE(0 == sut.millisecond);
        }

        SECTION("Hour wrap around")
        {
            CDateTime sut{ 2010, 11, 3, 22, 00, 02 };
            sut.Increment(-5);

            REQUIRE(2010 == sut.year);
            REQUIRE(11 == sut.month);
            REQUIRE(3 == sut.day);
            REQUIRE(21 == sut.hour);
            REQUIRE(59 == sut.minute);
            REQUIRE(57 == sut.second);
            REQUIRE(0 == sut.millisecond);
        }

        SECTION("Day wrap around")
        {
            CDateTime sut{ 2010, 11, 3, 0, 0, 2 };
            sut.Increment(-5);

            REQUIRE(2010 == sut.year);
            REQUIRE(11 == sut.month);
            REQUIRE(2 == sut.day);
            REQUIRE(23 == sut.hour);
            REQUIRE(59 == sut.minute);
            REQUIRE(57 == sut.second);
            REQUIRE(0 == sut.millisecond);
        }

        SECTION("Month wrap around")
        {
            CDateTime sut{ 2010, 12, 1, 0, 0, 2 };
            sut.Increment(-10);

            REQUIRE(2010 == sut.year);
            REQUIRE(11 == sut.month);
            REQUIRE(30 == sut.day);
            REQUIRE(23 == sut.hour);
            REQUIRE(59 == sut.minute);
            REQUIRE(52 == sut.second);
            REQUIRE(0 == sut.millisecond);
        }
    }

    TEST_CASE("CDateTime - Decrement with positive number of seconds.", "[DateTime]")
    {
        SECTION("Common time")
        {
            CDateTime sut{ 2010, 11, 3, 22, 50, 12 };
            sut.Decrement(5);

            REQUIRE(2010 == sut.year);
            REQUIRE(11 == sut.month);
            REQUIRE(3 == sut.day);
            REQUIRE(22 == sut.hour);
            REQUIRE(50 == sut.minute);
            REQUIRE(7 == sut.second);
            REQUIRE(0 == sut.millisecond);
        }

        SECTION("Minute wrap around")
        {
            CDateTime sut{ 2010, 11, 3, 22, 50, 2 };
            sut.Decrement(5);

            REQUIRE(2010 == sut.year);
            REQUIRE(11 == sut.month);
            REQUIRE(3 == sut.day);
            REQUIRE(22 == sut.hour);
            REQUIRE(49 == sut.minute);
            REQUIRE(57 == sut.second);
            REQUIRE(0 == sut.millisecond);
        }

        SECTION("Hour wrap around")
        {
            CDateTime sut{ 2010, 11, 3, 22, 0, 2 };
            sut.Decrement(5);

            REQUIRE(2010 == sut.year);
            REQUIRE(11 == sut.month);
            REQUIRE(3 == sut.day);
            REQUIRE(21 == sut.hour);
            REQUIRE(59 == sut.minute);
            REQUIRE(57 == sut.second);
            REQUIRE(0 == sut.millisecond);
        }

        SECTION("Day wrap around")
        {
            CDateTime sut{ 2010, 11, 3, 0, 0, 2 };
            sut.Decrement(5);

            REQUIRE(2010 == sut.year);
            REQUIRE(11 == sut.month);
            REQUIRE(2 == sut.day);
            REQUIRE(23 == sut.hour);
            REQUIRE(59 == sut.minute);
            REQUIRE(57 == sut.second);
            REQUIRE(0 == sut.millisecond);
        }

        SECTION("Month wrap around")
        {
            CDateTime sut{ 2010, 11, 1, 0, 0, 2 };
            sut.Decrement(10);

            REQUIRE(2010 == sut.year);
            REQUIRE(10 == sut.month);
            REQUIRE(31 == sut.day);
            REQUIRE(23 == sut.hour);
            REQUIRE(59 == sut.minute);
            REQUIRE(52 == sut.second);
            REQUIRE(0 == sut.millisecond);
        }
    }

    TEST_CASE("CDateTime - Decrement with negative number of seconds.", "[DateTime]")
    {
        SECTION("Common time")
        {
            CDateTime sut{ 2010, 11, 3, 22, 50, 12 };
            sut.Decrement(-5);

            REQUIRE(2010 == sut.year);
            REQUIRE(11 == sut.month);
            REQUIRE(3 == sut.day);
            REQUIRE(22 == sut.hour);
            REQUIRE(50 == sut.minute);
            REQUIRE(17 == sut.second);
            REQUIRE(0 == sut.millisecond);
        }

        SECTION("Minute wrap around")
        {
            CDateTime sut{ 2010, 11, 3, 22, 50, 56 };
            sut.Decrement(-5);

            REQUIRE(2010 == sut.year);
            REQUIRE(11 == sut.month);
            REQUIRE(3 == sut.day);
            REQUIRE(22 == sut.hour);
            REQUIRE(51 == sut.minute);
            REQUIRE(1 == sut.second);
            REQUIRE(0 == sut.millisecond);
        }

        SECTION("Hour wrap around")
        {
            CDateTime sut{ 2010, 11, 3, 22, 59, 57 };
            sut.Decrement(-5);

            REQUIRE(2010 == sut.year);
            REQUIRE(11 == sut.month);
            REQUIRE(3 == sut.day);
            REQUIRE(23 == sut.hour);
            REQUIRE(0 == sut.minute);
            REQUIRE(2 == sut.second);
            REQUIRE(0 == sut.millisecond);
        }

        SECTION("Day wrap around")
        {
            CDateTime sut{ 2010, 11, 3, 23, 59, 57 };
            sut.Decrement(-5);

            REQUIRE(2010 == sut.year);
            REQUIRE(11 == sut.month);
            REQUIRE(4 == sut.day);
            REQUIRE(0 == sut.hour);
            REQUIRE(0 == sut.minute);
            REQUIRE(2 == sut.second);
            REQUIRE(0 == sut.millisecond);
        }

        SECTION("Month wrap around")
        {
            CDateTime sut{ 2010, 12, 31, 23, 59, 57 };
            sut.Decrement(-10);

            REQUIRE(2011 == sut.year);
            REQUIRE(1 == sut.month);
            REQUIRE(1 == sut.day);
            REQUIRE(0 == sut.hour);
            REQUIRE(0 == sut.minute);
            REQUIRE(7 == sut.second);
            REQUIRE(0 == sut.millisecond);
        }
    }

    TEST_CASE("CDateTime - Print to string.", "[DateTime]")
    {
        CDateTime sut{ 2010, 6, 6, 22, 5, 12 };
        std::stringstream msg;

        // Act
        msg << sut;

        // Assert.
        REQUIRE("2010.06.06T22:05:12" == msg.str());
    }
}
