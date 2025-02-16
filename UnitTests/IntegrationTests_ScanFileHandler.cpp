#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/StringUtils.h>
#include "catch.hpp"
#include "TestData.h"

namespace novac
{

TEST_CASE("ScanFileHandler CheckScanFile", "[ScanFileHandler][IntegrationTests]")
{
    novac::ConsoleLog log;
    novac::LogContext context;
    CScanFileHandler sut(log);
    std::string file = TestData::GetMeasuredSpectrumName_I2J8549();
    sut.CheckScanFile(context, file);

    SECTION("Sets filename")
    {
        REQUIRE(EqualsIgnoringCase(sut.GetFileName(), file));
    }

    SECTION("Reads device")
    {
        REQUIRE(EqualsIgnoringCase(sut.m_device, "I2J8549"));
    }

    SECTION("Reads channel")
    {
        REQUIRE(sut.m_channel == 0);
    }

    SECTION("Reads start time")
    {
        REQUIRE(sut.m_startTime.year == 2017);
        REQUIRE(sut.m_startTime.month == 2);
        REQUIRE(sut.m_startTime.day == 16);

        REQUIRE(sut.m_startTime.hour == 12);
        REQUIRE(sut.m_startTime.minute == 30);
        REQUIRE(sut.m_startTime.second == 44);
    }

    SECTION("Reads stop time")
    {
        REQUIRE(sut.m_stopTime.year == 2017);
        REQUIRE(sut.m_stopTime.month == 2);
        REQUIRE(sut.m_stopTime.day == 16);

        REQUIRE(sut.m_stopTime.hour == 12);
        REQUIRE(sut.m_stopTime.minute == 30);
        REQUIRE(sut.m_stopTime.second == 56);
    }
}

TEST_CASE("ScanFileHandler GetSpectrumLength", "[ScanFileHandler][IntegrationTests]")
{
    novac::ConsoleLog log;
    novac::LogContext context;
    CScanFileHandler sut(log);
    std::string file = TestData::GetMeasuredSpectrumName_I2J8549();
    sut.CheckScanFile(context, file);

    SECTION("Gets length of spectra")
    {
        REQUIRE(sut.GetSpectrumLength() == 2048);
    }
}

TEST_CASE("ScanFileHandler GetSpectrumNumInFile", "[ScanFileHandler][IntegrationTests]")
{
    novac::ConsoleLog log;
    novac::LogContext context;
    CScanFileHandler sut(log);
    std::string file = TestData::GetMeasuredSpectrumName_I2J8549();
    sut.CheckScanFile(context, file);

    SECTION("Gets number of spectra")
    {
        REQUIRE(sut.GetSpectrumNumInFile() == 53);
    }
}

TEST_CASE("ScanFileHandler GetDark", "[ScanFileHandler][IntegrationTests]")
{
    novac::ConsoleLog log;
    novac::LogContext context;
    CScanFileHandler sut(log);
    std::string file = TestData::GetMeasuredSpectrumName_I2J8549();
    sut.CheckScanFile(context, file);

    SECTION("Reads the dark spectrum")
    {
        CSpectrum spec;
        int returnCode = sut.GetDark(spec);

        REQUIRE(returnCode == 0);
        REQUIRE(spec.ScanIndex() == 1);
    }
}

TEST_CASE("ScanFileHandler GetSky", "[ScanFileHandler][IntegrationTests]")
{
    novac::ConsoleLog log;
    novac::LogContext context;
    CScanFileHandler sut(log);
    std::string file = TestData::GetMeasuredSpectrumName_I2J8549();
    sut.CheckScanFile(context, file);

    SECTION("Reads the sky spectrum")
    {
        CSpectrum spec;
        int returnCode = sut.GetSky(spec);

        REQUIRE(returnCode == 0);
        REQUIRE(spec.ScanIndex() == 0);
    }
}

TEST_CASE("ScanFileHandler GetSpectrum", "[ScanFileHandler][IntegrationTests]")
{
    novac::ConsoleLog log;
    novac::LogContext context;
    CScanFileHandler sut(log);
    std::string file = TestData::GetMeasuredSpectrumName_I2J8549();
    sut.CheckScanFile(context, file);

    SECTION("Reads all spectra in the file (INCLUDING sky and dark) one at a time")
    {
        CSpectrum spec;
        int counter = 0;
        while (0 != sut.GetSpectrum(context, spec, counter))
        {
            ++counter;
            REQUIRE(spec.ScanIndex() >= 0);
        }

        REQUIRE(counter == 53);
    }
}

TEST_CASE("ScanFileHandler GetNextSpectrum", "[ScanFileHandler][IntegrationTests]")
{
    novac::ConsoleLog log;
    novac::LogContext context;
    CScanFileHandler sut(log);
    std::string file = TestData::GetMeasuredSpectrumName_I2J8549();
    sut.CheckScanFile(context, file);

    SECTION("Reads all spectra in the file (except sky and dark) one at a time")
    {
        CSpectrum spec;
        int counter = 0;
        while (0 != sut.GetNextSpectrum(context, spec))
        {
            ++counter;
            REQUIRE(spec.ScanIndex() >= 0);
        }

        REQUIRE(counter == 51);
    }

    SECTION("ResetCounter makes it possible to read all spectra again")
    {
        CSpectrum spec;
        int firstCounter = 0;
        while (0 != sut.GetNextSpectrum(context, spec))
        {
            ++firstCounter;
        }

        // Act
        sut.ResetCounter();

        // Assert by verifying that we can read the same number of spectra as above.
        int secondCounter = 0;
        while (0 != sut.GetNextSpectrum(context, spec))
        {
            ++secondCounter;
        }

        REQUIRE(firstCounter == secondCounter);
    }
}
}