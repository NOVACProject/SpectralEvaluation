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

    SECTION("Polynomial fit")
    {
        // Setup the sut
        RatioEvaluation sut{ settings, darkSettings };
        sut.SetupFirstResult(evaluationFileHandler.m_scan[0], plumeInScanProperties);
        sut.SetupFitWindows(so2FitWindow, std::vector<CFitWindow>{ broFitWindow});

        // Act
        const auto ratios = sut.Run(fileHandler);

        // Assert
        REQUIRE(ratios.size() == 1);

        // Verify that the result is indeed correct!
        const auto& result = ratios.front();
        REQUIRE(result.minorSpecieName == "BrO");
        REQUIRE(std::abs(result.minorResult - 2.1e14) < 2e13);

        REQUIRE(result.majorSpecieName == "SO2");
        REQUIRE(std::abs(result.majorResult - 2.2e18) < 1e17);

        REQUIRE(std::abs(result.ratio - 9.6e-5) < 1e-5);
        REQUIRE(std::abs(result.error - 2e-5) < 1e-5);
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
        RatioEvaluation sut{ settings, darkSettings };
        sut.SetupFirstResult(evaluationFileHandler.m_scan[0], plumeInScanProperties);
        sut.SetupFitWindows(so2FitWindow, std::vector<CFitWindow>{ broFitWindow});

        // Act
        const auto ratios = sut.Run(fileHandler);

        // Assert
        REQUIRE(ratios.size() == 1);

        // Verify that the result is indeed correct!
        const auto& result = ratios.front();
        REQUIRE(result.minorSpecieName == "BrO");
        REQUIRE(std::abs(result.minorResult - 2.8e14) < 2e13); // TODO: This is not exactly the same result as above. Why??

        REQUIRE(result.majorSpecieName == "SO2");
        REQUIRE(std::abs(result.majorResult - 2.2e18) < 1e17);

        REQUIRE(std::abs(result.ratio - 1.3e-4) < 1e-5); // TODO: This is not exactly the same result as above. Why??
        REQUIRE(std::abs(result.error - 3e-5) < 1e-5); // TODO: This is not exactly the same result as above. Why??
    }
}
