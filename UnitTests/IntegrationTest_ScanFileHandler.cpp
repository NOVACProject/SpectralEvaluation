#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/StringUtils.h>
#include "catch.hpp"

namespace novac
{
    std::string GetFileName()
    {
#ifdef _MSC_VER
        return std::string("../testData/I2J8549_170216_1230_0.pak");
#else
        return std::string("testData/I2J8549_170216_1230_0.pak");
#endif // _MSC_VER
    }

    TEST_CASE("ScanFileHandler CheckScanFile", "[ScanFileHandler][IntegrationTests]")
    {
        ::FileHandler::CScanFileHandler sut;
        std::string file = GetFileName();
        sut.CheckScanFile(file);

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
        }

    }

    TEST_CASE("ScanFileHandler GetSpectrumLength", "[ScanFileHandler][IntegrationTests]")
    {
        ::FileHandler::CScanFileHandler sut;
        std::string file = GetFileName();
        sut.CheckScanFile(file);

        SECTION("Gets length of spectra")
        {
            REQUIRE(sut.GetSpectrumLength() == 2048);
        }
    }

    TEST_CASE("ScanFileHandler GetSpectrumNumInFile", "[ScanFileHandler][IntegrationTests]")
    {
        ::FileHandler::CScanFileHandler sut;
        std::string file = GetFileName();
        sut.CheckScanFile(file);

        SECTION("Gets number of spectra")
        {
            REQUIRE(sut.GetSpectrumNumInFile() == 53);
        }
    }

    TEST_CASE("ScanFileHandler GetDark", "[ScanFileHandler][IntegrationTests]")
    {
        ::FileHandler::CScanFileHandler sut;
        std::string file = GetFileName();
        sut.CheckScanFile(file);

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
        ::FileHandler::CScanFileHandler sut;
        std::string file = GetFileName();
        sut.CheckScanFile(file);

        SECTION("Reads the sky spectrum")
        {
            CSpectrum spec;
            int returnCode = sut.GetSky(spec);

            REQUIRE(returnCode == 0);
            REQUIRE(spec.ScanIndex() == 0);
        }
    }

    TEST_CASE("ScanFileHandler GetNextSpectrum", "[ScanFileHandler][IntegrationTests]")
    {
        ::FileHandler::CScanFileHandler sut;
        std::string file = GetFileName();
        sut.CheckScanFile(file);

        SECTION("Reads all spectra in the file (except sky and dark) one at a time")
        {
            CSpectrum spec;
            int counter = 0;
            while (0 != sut.GetNextSpectrum(spec))
            {
                ++counter;
                REQUIRE(spec.ScanIndex() >= 0);
            }

            REQUIRE(counter == 51);
        }
    }
}