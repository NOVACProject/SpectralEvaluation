#include <SpectralEvaluation/Evaluation/PlumeSpectrumSelector.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h>
#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>
#include <SpectralEvaluation/File/ScanEvaluationLogFileHandler.h>
#include "catch.hpp"
#include "TestData.h"

namespace novac
{
    TEST_CASE("PlumeSpectrumSelector returns nullptr if input not readable", "[PlumeSpectrumSelector][IntegrationTest]")
    {
        novac::CScanFileHandler fileHandler;
        novac::BasicScanEvaluationResult evaluationResult;
        novac::CPlumeInScanProperty plumeInScanProperties;
        PlumeSpectrumSelector sut;

        SECTION("Invalid scan file")
        {
            const auto result = sut.CreatePlumeSpectra(fileHandler, evaluationResult, plumeInScanProperties);

            REQUIRE(nullptr == result);
        }

        SECTION("Invalid evaluation result")
        {
            REQUIRE(true == fileHandler.CheckScanFile(TestData::GetBrORatioScanFile1()));

            const auto result = sut.CreatePlumeSpectra(fileHandler, evaluationResult, plumeInScanProperties);

            REQUIRE(nullptr == result);
        }
    }

    // happy case where there is a good plume and we should be able to extract in plume and out of plume spectra
    TEST_CASE("PlumeSpectrumSelector returns expected in and out of plume for good scan - scan file 1", "[PlumeSpectrumSelector][IntegrationTest]")
    {
        novac::CScanFileHandler fileHandler;
        novac::CScanEvaluationLogFileHandler evaluationFileHandler;
        novac::BasicScanEvaluationResult evaluationResult;
        PlumeSpectrumSelector sut;

        SECTION("BrO Scan file 1")
        {
            const bool scanFileIsOk = fileHandler.CheckScanFile(TestData::GetBrORatioScanFile1());
            REQUIRE(scanFileIsOk); // check assumption on the setup

            const bool evaluationFileIsOk = evaluationFileHandler.ReadEvaluationLog(TestData::GetBrORatioEvaluationFile1());
            REQUIRE(evaluationFileIsOk); // check assumption on the setup
            REQUIRE(evaluationFileHandler.m_scan.size() == 1); // check assumption on the setup

            novac::CPlumeInScanProperty plumeInScanProperties;
            REQUIRE(true == novac::FindPlume(evaluationFileHandler.m_scan[0], 0, plumeInScanProperties));
            novac::CalculatePlumeOffset(evaluationFileHandler.m_scan[0], 0, plumeInScanProperties);
            REQUIRE(true == novac::CalculatePlumeCompleteness(evaluationFileHandler.m_scan[0], 0, plumeInScanProperties));

            const auto result = sut.CreatePlumeSpectra(fileHandler, evaluationFileHandler.m_scan[0], plumeInScanProperties);

            // Assert
            REQUIRE(result != nullptr);
            REQUIRE(result->inPlumeSpectrum != nullptr);
            REQUIRE(result->referenceSpectrum != nullptr);
            REQUIRE(result->darkSpectrum != nullptr);
        }
    }

    // handling the case where there isn't a good or visible plume in the scan


}
