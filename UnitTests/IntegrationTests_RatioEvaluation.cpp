#include <SpectralEvaluation/Configuration/DarkSettings.h>
#include <SpectralEvaluation/Evaluation/RatioEvaluation.h>
#include <SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h>
#include <SpectralEvaluation/File/FitWindowFileHandler.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/File/ScanEvaluationLogFileHandler.h>
#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>
#include "catch.hpp"
#include "TestData.h"

using namespace novac;

TEST_CASE("RatioEvaluation - IntegrationTest with good scan - scan file 1", "[RatioEvaluation][IntegrationTest]")
{
    RatioEvaluationSettings settings;
    Configuration::CDarkSettings darkSettings;

    // Prepare the test by reading in the .pak-file and the evaluation result and calculate the plume-properties from the result.
    novac::CScanFileHandler fileHandler;
    const bool scanFileIsOk = fileHandler.CheckScanFile(TestData::GetBrORatioScanFile1());
    REQUIRE(scanFileIsOk); // check assumption on the setup

    novac::CScanEvaluationLogFileHandler evaluationFileHandler;
    const bool evaluationFileIsOk = evaluationFileHandler.ReadEvaluationLog(TestData::GetBrORatioEvaluationFile1());
    REQUIRE(evaluationFileIsOk); // check assumption on the setup
    REQUIRE(evaluationFileHandler.m_scan.size() == 1); // check assumption on the setup

    novac::CPlumeInScanProperty plumeInScanProperties;
    novac::CalculatePlumeOffset(evaluationFileHandler.m_scan[0], 0, plumeInScanProperties);
    REQUIRE(true == novac::CalculatePlumeCompleteness(evaluationFileHandler.m_scan[0], 0, plumeInScanProperties));

    // Read in the fit windows to use
    CFitWindowFileHandler fitWindowFileHandler;
    auto allWindows = fitWindowFileHandler.ReadFitWindowFile(TestData::GetBrORatioFitWindowFileSO2());
    REQUIRE(allWindows.size() == 1);
    auto so2FitWindow = allWindows.front();

    allWindows = fitWindowFileHandler.ReadFitWindowFile(TestData::GetBrORatioFitWindowFileBrO());
    REQUIRE(allWindows.size() == 1);
    auto broFitWindow = allWindows.front();

    // Read in the references
    REQUIRE(true == ReadReferences(so2FitWindow));
    REQUIRE(true == ReadReferences(broFitWindow));

    // Setup the sut
    RatioEvaluation sut{ settings, darkSettings };
    sut.SetupFirstResult(evaluationFileHandler.m_scan[0], plumeInScanProperties);
    sut.SetupFitWindows(so2FitWindow, std::vector<CFitWindow>{ broFitWindow});

    std::string errorMessage;

    // Act
    const auto result = sut.Run(fileHandler, &errorMessage);

    // Assert
    REQUIRE(result.size() == 1);
    REQUIRE(errorMessage.empty());

    // TODO: Verify that the result is indeed correct!

}