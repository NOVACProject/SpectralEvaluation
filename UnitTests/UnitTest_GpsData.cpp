#include "catch.hpp"
#include <SpectralEvaluation/GPSData.h>

TEST_CASE("GpsMath::Distance", "[GpsData]")
{
    SECTION("Same source and destionation return zero distance")
    {
        REQUIRE(Approx(0.00) == novac::GpsMath::Distance(57.6898586, 11.9725077, 57.6898586, 11.9725077));
    }

    SECTION("Distance from Eifel tower to Arc de Triomphe")
    {
        const double sourceLat = 48.8583905296204;
        const double sourceLon = 2.2944259643554688;
        const double destinationLat = 48.87375689312557;
        const double destinationLon = 2.2950005773730364;

        const double fromEifelTowerToArcDeTriomphe = novac::GpsMath::Distance(sourceLat, sourceLon, destinationLat, destinationLon);

        REQUIRE(Approx(1708.08) == fromEifelTowerToArcDeTriomphe);
    }

    SECTION("From A to B equals distance from B to A")
    {
        const double sourceLat = 48.8583701;
        const double sourceLon = 2.2944813;
        const double destinationLat = 48.87375689312557;
        const double destinationLon = 2.2950005773730364;

        const double distanceFromSourceToDestination = novac::GpsMath::Distance(sourceLat, sourceLon, destinationLat, destinationLon);
        const double distanceFromDestinationToSource = novac::GpsMath::Distance(destinationLat, destinationLon, sourceLat, sourceLon);

        REQUIRE(Approx(distanceFromSourceToDestination) == distanceFromDestinationToSource);
    }
}


TEST_CASE("GpsMath::Bearing", "[GpsData]")
{
    SECTION("Same source and destionation return zero bearing")
    {
        REQUIRE(Approx(0.00) == novac::GpsMath::Bearing(57.6898586, 11.9725077, 57.6898586, 11.9725077));
    }

    SECTION("Bearing straight north is zero")
    {
        const double sourceLat = 48.8583905296204;
        const double sourceLon = 2.2944259643554688;

        const double bearing = novac::GpsMath::Bearing(sourceLat, sourceLon, sourceLat + 0.01, sourceLon);

        REQUIRE(std::abs(bearing) < 0.01);
    }

    SECTION("Bearing straight south is 180")
    {
        const double sourceLat = 48.8583905296204;
        const double sourceLon = 2.2944259643554688;

        const double bearing = novac::GpsMath::Bearing(sourceLat, sourceLon, sourceLat - 0.01, sourceLon);

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


    SECTION("Bearing from Eifel tower to Arc de Triomphe")
    {
        const double sourceLat = 48.8583905296204;
        const double sourceLon = 2.2944259643554688;
        const double destinationLat = 48.87375689312557;
        const double destinationLon = 2.2950005773730364;

        const double fromEifelTowerToArcDeTriomphe = novac::GpsMath::Bearing(sourceLat, sourceLon, destinationLat, destinationLon);

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