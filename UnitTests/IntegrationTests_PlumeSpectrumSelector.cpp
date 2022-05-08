#include <SpectralEvaluation/Evaluation/PlumeSpectrumSelector.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h>
#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>
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


}
