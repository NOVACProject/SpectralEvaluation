#include <SpectralEvaluation/File/ScanEvaluationLogFileHandler.h>
#include <SpectralEvaluation/StringUtils.h>
#include "catch.hpp"
#include "TestData.h"

namespace novac
{
TEST_CASE("CScanEvaluationLogFileHandler ReadEvaluationLog with one scan in file", "[CScanEvaluationLogFileHandler][IntegrationTests]")
{
    CScanEvaluationLogFileHandler sut;
    std::string file = TestData::GetBrORatioEvaluationFile1();

    // Act
    sut.ReadEvaluationLog(file);

    SECTION("Sets m_evaluationLog")
    {
        REQUIRE(EqualsIgnoringCase(sut.m_evaluationLog, file));
    }

    SECTION("Reads one scan")
    {
        REQUIRE(1 == sut.m_scan.size());
    }

    SECTION("Scan has correct number of spectra")
    {
        // There are 51 positions in the file, but sky and dark should be removed...
        REQUIRE(51 == sut.m_scan[0].m_spec.size());
    }

    SECTION("Reads expected scan properties")
    {
        const auto& firstEvaluatedSpectrum = sut.m_scan[0].m_spec[0];

        // There should be three species evaluated for (SO2, O3 and Ring)
        REQUIRE(3 == firstEvaluatedSpectrum.m_referenceResult.size());
        REQUIRE("so2" == firstEvaluatedSpectrum.m_referenceResult[0].m_specieName);
        REQUIRE("o3" == firstEvaluatedSpectrum.m_referenceResult[1].m_specieName);
        REQUIRE("ring" == firstEvaluatedSpectrum.m_referenceResult[2].m_specieName);
    }

    SECTION("Reads correct evaluation results for first spectrum")
    {
        const auto& firstEvaluatedSpectrum = sut.m_scan[0].m_spec[0];

        // There should be three species evaluated for (SO2, O3 and Ring)
        REQUIRE(3 == firstEvaluatedSpectrum.m_referenceResult.size());
        REQUIRE(Approx(1.56e+18) == firstEvaluatedSpectrum.m_referenceResult[0].m_column);
        REQUIRE(Approx(-4.71e+18) == firstEvaluatedSpectrum.m_referenceResult[1].m_column);
        REQUIRE(Approx(-1.53e+15) == firstEvaluatedSpectrum.m_referenceResult[2].m_column);

        REQUIRE(Approx(9.17e+17) == firstEvaluatedSpectrum.m_referenceResult[0].m_columnError);
        REQUIRE(Approx(2.32e+19) == firstEvaluatedSpectrum.m_referenceResult[1].m_columnError);
        REQUIRE(Approx(7.33e+14) == firstEvaluatedSpectrum.m_referenceResult[2].m_columnError);
    }

    SECTION("Reads correct meta data")
    {
        const auto& firstScan = sut.m_scan[0];

        REQUIRE(firstScan.m_spec.size() == firstScan.m_specInfo.size());

        REQUIRE("D2J2124" == firstScan.m_specInfo[0].m_device);
        REQUIRE("D2J2124" == firstScan.m_specInfo[50].m_device);

        REQUIRE(CDateTime(2016, 3, 31, 16, 8, 44) == firstScan.m_skySpecInfo.m_startTime);
        REQUIRE(CDateTime(2016, 3, 31, 16, 9, 42) == firstScan.m_specInfo[0].m_startTime);
        REQUIRE(CDateTime(2016, 3, 31, 16, 14, 58) == firstScan.m_specInfo[50].m_startTime);

        REQUIRE(0.15F == firstScan.m_specInfo[0].m_peakIntensity);
        REQUIRE(0.78F == firstScan.m_specInfo[50].m_peakIntensity);

        REQUIRE(0.09F == firstScan.m_specInfo[0].m_fitIntensity);
        REQUIRE(0.19F == firstScan.m_specInfo[50].m_fitIntensity);
    }

    SECTION("Reads correct evaluation results for last spectrum")
    {
        const auto& lastEvaluatedSpectrum = sut.m_scan[0].m_spec[50];

        // There should be three species evaluated for (SO2, O3 and Ring)
        REQUIRE(3 == lastEvaluatedSpectrum.m_referenceResult.size());
        REQUIRE(Approx(-1.23e+18) == lastEvaluatedSpectrum.m_referenceResult[0].m_column);
        REQUIRE(Approx(4.26e+18) == lastEvaluatedSpectrum.m_referenceResult[1].m_column);
        REQUIRE(Approx(-1.33e+11) == lastEvaluatedSpectrum.m_referenceResult[2].m_column);

        REQUIRE(Approx(9.66e+16) == lastEvaluatedSpectrum.m_referenceResult[0].m_columnError);
        REQUIRE(Approx(2.44e+18) == lastEvaluatedSpectrum.m_referenceResult[1].m_columnError);
        REQUIRE(Approx(7.72e+13) == lastEvaluatedSpectrum.m_referenceResult[2].m_columnError);
    }

    SECTION("Sets correct evaluation judgement on spectra")
    {
        const auto& firstScan = sut.m_scan[0];

        REQUIRE(true == firstScan.m_spec[0].IsBad());
        REQUIRE(true == firstScan.m_spec[1].IsBad());
        REQUIRE(true == firstScan.m_spec[2].IsBad());
        REQUIRE(true == firstScan.m_spec[3].IsBad());
        REQUIRE(true == firstScan.m_spec[4].IsBad());
        REQUIRE(true == firstScan.m_spec[5].IsBad());
        REQUIRE(true == firstScan.m_spec[6].IsBad());
        REQUIRE(true == firstScan.m_spec[7].IsBad());
        REQUIRE(true == firstScan.m_spec[8].IsBad());
        REQUIRE(true == firstScan.m_spec[9].IsBad());
        REQUIRE(true == firstScan.m_spec[10].IsBad());
        REQUIRE(false == firstScan.m_spec[11].IsBad());
        REQUIRE(false == firstScan.m_spec[12].IsBad());
        REQUIRE(false == firstScan.m_spec[13].IsBad());
    }
}

TEST_CASE("CScanEvaluationLogFileHandler ReadEvaluationLog with multiple scans in one file", "[CScanEvaluationLogFileHandler][IntegrationTests]")
{
    CScanEvaluationLogFileHandler sut;
    const std::string file = TestData::GetEvaluationLogfile1();

    // Act
    sut.ReadEvaluationLog(file);

    SECTION("Sets m_evaluationLog")
    {
        REQUIRE(EqualsIgnoringCase(sut.m_evaluationLog, file));
    }

    SECTION("Reads 69 scans")
    {
        REQUIRE(69 == sut.m_scan.size());
    }

    SECTION("The scans are sorted in increasing time order")
    {
        for (size_t idx = 1; idx < sut.m_scan.size(); ++idx)
        {
            // allowing for equal times here, as there are some duplicates in the file
            REQUIRE(sut.m_scan[idx].m_skySpecInfo.m_startTime >= sut.m_scan[idx - 1].m_skySpecInfo.m_startTime);
        }
    }
}

TEST_CASE("CScanEvaluationLogFileHandler ReadEvaluationLog - ReEvaluationLog from NovacProgram", "[CScanEvaluationLogFileHandler][IntegrationTests]")
{
    CScanEvaluationLogFileHandler sut;
    const std::string file = TestData::GetEvaluationLogfile2();

    // Act
    sut.ReadEvaluationLog(file);

    SECTION("Sets m_evaluationLog")
    {
        REQUIRE(EqualsIgnoringCase(sut.m_evaluationLog, file));
    }

    SECTION("Reads 78 scans")
    {
        REQUIRE(78 == sut.m_scan.size());
    }

    SECTION("The scans are sorted in increasing time order")
    {
        for (size_t idx = 1; idx < sut.m_scan.size(); ++idx)
        {
            // allowing for equal times here, as there are some duplicates in the file
            REQUIRE(sut.m_scan[idx].m_skySpecInfo.m_startTime >= sut.m_scan[idx - 1].m_skySpecInfo.m_startTime);
        }
    }
}
}