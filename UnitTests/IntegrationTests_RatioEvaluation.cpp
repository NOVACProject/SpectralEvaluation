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

TEST_CASE("RatioEvaluation - IntegrationTest with good scan - scan file 1", "[RatioEvaluation][Ratios][IntegrationTest]")
{
    Configuration::RatioEvaluationSettings settings;
    Configuration::CDarkSettings darkSettings;

    // Prepare the test by reading in the .pak-file and the evaluation result and calculate the plume-properties from the result.
    novac::ConsoleLog log;
    novac::LogContext context;
    CScanFileHandler fileHandler(log);
    const bool scanFileIsOk = fileHandler.CheckScanFile(context, TestData::GetBrORatioScanFile1());
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
    so2FitWindow.fitLow = 439;
    so2FitWindow.fitHigh = 592;

    allWindows = fitWindowFileHandler.ReadFitWindowFile(TestData::GetBrORatioFitWindowFileBrO());
    REQUIRE(allWindows.size() == 1);
    auto broFitWindow = allWindows.front();
    broFitWindow.fitLow = 641;
    broFitWindow.fitHigh = 939;

    // Read in the references
    ReadReferences(so2FitWindow);
    ReadReferences(broFitWindow);

    SECTION("Polynomial fit")
    {
        // Setup the sut
        RatioEvaluation sut{ settings, darkSettings, log };
        sut.SetupFirstResult(evaluationFileHandler.m_scan[0], plumeInScanProperties);
        sut.SetupFitWindows(so2FitWindow, std::vector<CFitWindow>{ broFitWindow});

        // Act
        const auto ratios = sut.Run(context, fileHandler);

        // Assert
        REQUIRE(ratios.size() == 1);

        // Verify that the result is indeed correct!
        const auto& result = ratios.front();
        REQUIRE(result.minorSpecieName == "BrO");
        REQUIRE(result.minorResult == Approx(1.47e14).margin(2e13));

        REQUIRE(result.majorSpecieName == "SO2");
        REQUIRE(result.majorResult == Approx(2.0e18).margin(1e17));

        REQUIRE(result.ratio == Approx(7.2e-5).margin(1e-6));
        REQUIRE(result.error == Approx(1.4e-5).margin(1e-6));
    }

    SECTION("HP500 fit")
    {
        // Change the settings to use HP500 filtering and make sure to filter the references 
        // (but keep the unit in molec/cm2 as this makes it possible to have the same values in the assertions here as in the test above).
        so2FitWindow.fitType = FIT_TYPE::FIT_HP_DIV;
        broFitWindow.fitType = FIT_TYPE::FIT_HP_DIV;
        for (int refIdx = 0; refIdx < so2FitWindow.nRef; ++refIdx)
        {
            HighPassFilter(*so2FitWindow.ref[refIdx].m_data, false);
        }
        for (int refIdx = 0; refIdx < broFitWindow.nRef; ++refIdx)
        {
            HighPassFilter(*broFitWindow.ref[refIdx].m_data, false);
        }

        // Setup the sut
        RatioEvaluation sut{ settings, darkSettings, log };
        sut.SetupFirstResult(evaluationFileHandler.m_scan[0], plumeInScanProperties);
        sut.SetupFitWindows(so2FitWindow, std::vector<CFitWindow>{ broFitWindow});

        // Act
        const auto ratios = sut.Run(context, fileHandler);

        // Assert
        REQUIRE(ratios.size() == 1);

        // Verify that the result is indeed correct!
        const auto& result = ratios.front();
        REQUIRE(result.minorSpecieName == "BrO");
        REQUIRE(result.minorResult == Approx(1.66e14).margin(2e13));

        REQUIRE(result.majorSpecieName == "SO2");
        REQUIRE(result.majorResult == Approx(2.0e18).margin(1e17));

        REQUIRE(result.ratio == Approx(8.0e-5).margin(1e-6));
        REQUIRE(result.error == Approx(1.8e-5).margin(1e-6));
    }
}
