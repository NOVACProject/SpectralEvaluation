#include "catch.hpp"
#include <SpectralEvaluation/GPSData.h>

TEST_CASE("GpsMath::Distance", "[GpsData]")
{
    SECTION("Same source and destionation return zero distance (doubles)")
    {
        REQUIRE(Approx(0.00) == novac::GpsMath::Distance(57.6898586, 11.9725077, 57.6898586, 11.9725077));
    }

    SECTION("Same source and destionation return zero distance (structs)")
    {
        novac::CGPSData pos(57.6898586, 11.9725077, 0.0);
        REQUIRE(Approx(0.00) == novac::GpsMath::Distance(pos, pos));
    }

    SECTION("Distance from Eifel tower to Arc de Triomphe (double)")
    {
        const double sourceLat = 48.8583905296204;
        const double sourceLon = 2.2944259643554688;
        const double destinationLat = 48.87375689312557;
        const double destinationLon = 2.2950005773730364;

        const double fromEifelTowerToArcDeTriomphe = novac::GpsMath::Distance(sourceLat, sourceLon, destinationLat, destinationLon);

        REQUIRE(Approx(1708.08) == fromEifelTowerToArcDeTriomphe);
    }

    SECTION("Distance from Eifel tower to Arc de Triomphe (struct)")
    {
        const novac::CGPSData source(48.8583905296204, 2.2944259643554688, 100);
        const novac::CGPSData destination(48.87375689312557, 2.2950005773730364, 100);

        const double fromEifelTowerToArcDeTriomphe = novac::GpsMath::Distance(source, destination);

        REQUIRE(Approx(1708.08) == fromEifelTowerToArcDeTriomphe);
    }

    SECTION("From A to B equals distance from B to A (double)")
    {
        const double sourceLat = 48.8583701;
        const double sourceLon = 2.2944813;
        const double destinationLat = 48.87375689312557;
        const double destinationLon = 2.2950005773730364;

        const double distanceFromSourceToDestination = novac::GpsMath::Distance(sourceLat, sourceLon, destinationLat, destinationLon);
        const double distanceFromDestinationToSource = novac::GpsMath::Distance(destinationLat, destinationLon, sourceLat, sourceLon);

        REQUIRE(Approx(distanceFromSourceToDestination) == distanceFromDestinationToSource);
    }

    SECTION("From A to B equals distance from B to A (struct)")
    {
        const novac::CGPSData source(48.8583701, 2.2944813, 24);
        const novac::CGPSData destination(48.87375689312557, 2.2950005773730364, 24);

        const double distanceFromSourceToDestination = novac::GpsMath::Distance(source, destination);
        const double distanceFromDestinationToSource = novac::GpsMath::Distance(destination, source);

        REQUIRE(Approx(distanceFromSourceToDestination) == distanceFromDestinationToSource);
    }
}

TEST_CASE("GpsMath::Bearing", "[GpsData]")
{
    SECTION("Same source and destionation return zero bearing (double)")
    {
        REQUIRE(Approx(0.00) == novac::GpsMath::Bearing(57.6898586, 11.9725077, 57.6898586, 11.9725077));
    }

    SECTION("Same source and destionation return zero bearing (struct)")
    {
        const novac::CGPSData pos(57.6898586, 11.9725077, 50);
        REQUIRE(Approx(0.00) == novac::GpsMath::Bearing(pos, pos));
    }

    SECTION("Bearing straight north is zero (double)")
    {
        const double sourceLat = 48.8583905296204;
        const double sourceLon = 2.2944259643554688;

        const double bearing = novac::GpsMath::Bearing(sourceLat, sourceLon, sourceLat + 0.01, sourceLon);

        REQUIRE(std::abs(bearing) < 0.01);
    }

    SECTION("Bearing straight north is zero (struct)")
    {
        const novac::CGPSData source(48.8583905296204, 2.2944259643554688, 20);
        const novac::CGPSData destination(48.8583905296204 + 0.01, 2.2944259643554688, 20);

        const double bearing = novac::GpsMath::Bearing(source, destination);

        REQUIRE(std::abs(bearing) < 0.01);
    }

    SECTION("Bearing straight south is 180 (double)")
    {
        const double sourceLat = 48.8583905296204;
        const double sourceLon = 2.2944259643554688;

        const double bearing = novac::GpsMath::Bearing(sourceLat, sourceLon, sourceLat - 0.01, sourceLon);

        REQUIRE(std::abs(bearing - 180.0) < 0.01);
    }

    SECTION("Bearing straight south is 180 (struct)")
    {
        const novac::CGPSData source(48.8583905296204, 2.2944259643554688, 20);
        const novac::CGPSData destination(48.8583905296204 - 0.01, 2.2944259643554688, 20);

        const double bearing = novac::GpsMath::Bearing(source, destination);

        REQUIRE(std::abs(bearing - 180.0) < 0.01);
    }

    SECTION("Bearing straight west is 270")
    {
        const double sourceLat = 48.8583905296204;
        const double sourceLon = 2.2944259643554688;

        const double bearing = novac::GpsMath::Bearing(sourceLat, sourceLon, sourceLat, sourceLon - 0.01);

        REQUIRE(std::abs(bearing - 270.0) < 0.01);
    }

    SECTION("Bearing straight east is 90")
    {
        const double sourceLat = 48.8583905296204;
        const double sourceLon = 2.2944259643554688;

        const double bearing = novac::GpsMath::Bearing(sourceLat, sourceLon, sourceLat, sourceLon + 0.01);

        REQUIRE(std::abs(bearing - 90.0) < 0.01);
    }


    SECTION("Bearing from Eifel tower to Arc de Triomphe (double)")
    {
        const double sourceLat = 48.8583905296204;
        const double sourceLon = 2.2944259643554688;
        const double destinationLat = 48.87375689312557;
        const double destinationLon = 2.2950005773730364;

        const double fromEifelTowerToArcDeTriomphe = novac::GpsMath::Bearing(sourceLat, sourceLon, destinationLat, destinationLon);

        REQUIRE(std::abs(fromEifelTowerToArcDeTriomphe - 1.409) < 0.001);
    }

    SECTION("Bearing from Eifel tower to Arc de Triomphe (struct)")
    {
        const novac::CGPSData source(48.8583905296204, 2.2944259643554688, 20);
        const novac::CGPSData destination(48.87375689312557, 2.2950005773730364, 20);

        const double fromEifelTowerToArcDeTriomphe = novac::GpsMath::Bearing(source, destination);

        REQUIRE(std::abs(fromEifelTowerToArcDeTriomphe - 1.409) < 0.001);
    }

    SECTION("From A to B equals 180 minus bearing from B to A")
    {
        const double sourceLat = 48.8583701;
        const double sourceLon = 2.2944813;
        const double destinationLat = 48.87375689312557;
        const double destinationLon = 2.2950005773730364;

        const double bearingFromSourceToDestination = novac::GpsMath::Bearing(sourceLat, sourceLon, destinationLat, destinationLon);
        const double bearingFromDestinationToSource = novac::GpsMath::Bearing(destinationLat, destinationLon, sourceLat, sourceLon);

        REQUIRE(Approx(180.0 + bearingFromSourceToDestination) == bearingFromDestinationToSource);
    }
}

TEST_CASE("GpsData::DoubleToAngle", "[GpsData]")
{
    SECTION("Integer degree with zero minutes")
    {
        // Act
        const double result = novac::CGPSData::DoubleToAngle(5700.00);

        // Assert
        REQUIRE(Approx(result) == 57.00);
    }

    SECTION("Negative integer degree with zero minutes")
    {
        // Act
        const double result = novac::CGPSData::DoubleToAngle(-5700.00);

        // Assert
        REQUIRE(Approx(result) == -57.00);
    }

    SECTION("Degrees with 30 minutes")
    {
        // Act, 57 degrees and 30 minutes equals 57 and a half degrees
        const double result = novac::CGPSData::DoubleToAngle(5730.00);

        // Assert
        REQUIRE(Approx(result) == 57.50);
    }

    SECTION("Negative degrees with 30 minutes")
    {
        // Act, 57 degrees and 30 minutes equals 57 and a half degrees
        const double result = novac::CGPSData::DoubleToAngle(-5730.00);

        // Assert
        REQUIRE(Approx(result) == -57.50);
    }

    SECTION("Degrees with 15 minutes")
    {
        // Act, 48 degrees and 15 minutes equals 48 and a quarter degrees
        const double result = novac::CGPSData::DoubleToAngle(4815.00);

        // Assert
        REQUIRE(Approx(result) == 48.25);
    }

    SECTION("Negative degrees with 15 minutes")
    {
        // Act, 48 degrees and 15 minutes equals 48 and a quarter degrees
        const double result = novac::CGPSData::DoubleToAngle(-4815.00);

        // Assert
        REQUIRE(Approx(result) == -48.25);
    }
}