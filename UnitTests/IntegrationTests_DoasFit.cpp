#include <SpectralEvaluation/Evaluation/DoasFit.h>
#include <SpectralEvaluation/Evaluation/DoasFitPreparation.h>
#include <SpectralEvaluation/File/FitWindowFileHandler.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>
#include "catch.hpp"
#include "TestData.h"

using namespace novac;

TEST_CASE("DoasFit - IntegrationTest with good scan - scan file 1", "[DoasFit][IntegrationTest]")
{
    // Prepare the test by reading in the .pak-file and the evaluation result and calculate the plume-properties from the result.
    novac::ConsoleLog log;
    novac::LogContext context;
    novac::CScanFileHandler fileHandler(log);
    const bool scanFileIsOk = fileHandler.CheckScanFile(context, TestData::GetBrORatioScanFile1());
    REQUIRE(scanFileIsOk); // check assumption on the setup

    // Read in the fit window to use
    CFitWindowFileHandler fitWindowFileHandler;
    auto allWindows = fitWindowFileHandler.ReadFitWindowFile(TestData::GetBrORatioFitWindowFileSO2());
    REQUIRE(allWindows.size() == 1);
    auto so2FitWindow = allWindows.front();

    // Read in the references
    ReadReferences(so2FitWindow);

    // Read in the spectra to use and dark-correct them
    CSpectrum darkSpectrum;
    fileHandler.GetDark(darkSpectrum);
    CSpectrum skySpectrum;
    fileHandler.GetSky(skySpectrum);
    skySpectrum.Sub(darkSpectrum);
    DoasFitPreparation::RemoveOffset(skySpectrum);
    CSpectrum measuredSpectrum;
    fileHandler.GetSpectrum(context, 42, measuredSpectrum);
    measuredSpectrum.Sub(darkSpectrum);

    SECTION("Polynomial fit")
    {
        // Setup the DOAS Fit
        DoasFit sut;
        auto filteredSkySpectrum = DoasFitPreparation::PrepareSkySpectrum(skySpectrum, so2FitWindow.fitType);
        AddAsSky(so2FitWindow, filteredSkySpectrum, SHIFT_TYPE::SHIFT_FREE);
        auto filteredMeasuredData = DoasFitPreparation::PrepareMeasuredSpectrum(measuredSpectrum, skySpectrum, so2FitWindow.fitType);
        sut.Setup(so2FitWindow);

        // Act
        DoasResult result;
        sut.Run(filteredMeasuredData.data(), filteredMeasuredData.size(), result);

        // Assert
        // Verify the setup
        REQUIRE(result.fitLow == so2FitWindow.fitLow);
        REQUIRE(result.fitHigh == so2FitWindow.fitHigh);
        REQUIRE(result.referenceResult.size() == 3);
        REQUIRE(static_cast<int>(result.polynomialCoefficients.size()) == so2FitWindow.polyOrder + 1);
        REQUIRE(static_cast<int>(result.polynomialValues.size()) == so2FitWindow.fitHigh - so2FitWindow.fitLow);

        // Verify that the result is indeed correct!
        REQUIRE(std::abs(result.chiSquare - 0.01) < 0.005);
        REQUIRE(std::abs(result.delta - 0.04) < 0.02);

        REQUIRE(result.referenceResult[0].name == "SO2");
        REQUIRE(std::abs(result.referenceResult[0].column + 2.2e18) < 1e17); // the SO2 column result is negative here
        REQUIRE(std::abs(result.referenceResult[0].columnError - 1.0e17) < 1e16);

        REQUIRE(result.referenceResult[1].name == "O3");
        REQUIRE(std::abs(result.referenceResult[1].column - 1.1e17) < 1e16);
        REQUIRE(std::abs(result.referenceResult[1].columnError - 2.2e17) < 1e16);

        REQUIRE(result.referenceResult[2].name == "sky");
        REQUIRE(std::abs(result.referenceResult[2].column + 1) < std::numeric_limits<double>::epsilon());
        REQUIRE(std::abs(result.referenceResult[2].columnError) < std::numeric_limits<double>::epsilon());
    }

    SECTION("HP_DIV fit")
    {
        // Change the settings to use HP500 filtering and make sure to filter the references
        // (but keep the unit in molec/cm2 as this makes it possible to have the same values in the assertions here as in the test above).
        so2FitWindow.fitType = FIT_TYPE::FIT_HP_DIV;
        for (size_t refIdx = 0; refIdx < so2FitWindow.nRef; ++refIdx)
        {
            HighPassFilter(*so2FitWindow.ref[refIdx].m_data, CrossSectionUnit::cm2_molecule);
        }

        // Setup the DOAS Fit
        DoasFit sut;
        auto filteredMeasuredData = DoasFitPreparation::PrepareMeasuredSpectrum(measuredSpectrum, skySpectrum, so2FitWindow.fitType);
        sut.Setup(so2FitWindow);

        // Act
        DoasResult result;
        sut.Run(filteredMeasuredData.data(), filteredMeasuredData.size(), result);

        // Assert
        // Verify the setup
        REQUIRE(result.fitLow == so2FitWindow.fitLow);
        REQUIRE(result.fitHigh == so2FitWindow.fitHigh);
        REQUIRE(result.referenceResult.size() == 2);
        REQUIRE(static_cast<int>(result.polynomialCoefficients.size()) == so2FitWindow.polyOrder + 1);
        REQUIRE(static_cast<int>(result.polynomialValues.size()) == so2FitWindow.fitHigh - so2FitWindow.fitLow);

        // Verify that the result is indeed correct!
        REQUIRE(std::abs(result.chiSquare - 0.01) < 0.005);
        REQUIRE(std::abs(result.delta - 0.04) < 0.02);

        REQUIRE(result.referenceResult[0].name == "SO2");
        REQUIRE(std::abs(result.referenceResult[0].column + 2.2e18) < 1e17); // the SO2 column result is negative here
        REQUIRE(std::abs(result.referenceResult[0].columnError - 1.0e17) < 1e16);

        REQUIRE(result.referenceResult[1].name == "O3");
        REQUIRE(std::abs(result.referenceResult[1].column - 1.1e17) < 3e16);
        REQUIRE(std::abs(result.referenceResult[1].columnError - 2.2e17) < 1e16);
    }

    SECTION("HP_SUB fit")
    {
        // Change the settings to use HP500 filtering and make sure to filter the references
        // (but keep the unit in molec/cm2 as this makes it possible to have the same values in the assertions here as in the test above).
        so2FitWindow.fitType = FIT_TYPE::FIT_HP_SUB;
        for (size_t refIdx = 0; refIdx < so2FitWindow.nRef; ++refIdx)
        {
            HighPassFilter(*so2FitWindow.ref[refIdx].m_data, novac::CrossSectionUnit::cm2_molecule);
        }

        // Setup the DOAS Fit
        DoasFit sut;
        auto filteredSkySpectrum = DoasFitPreparation::PrepareSkySpectrum(skySpectrum, so2FitWindow.fitType);
        AddAsSky(so2FitWindow, filteredSkySpectrum, SHIFT_TYPE::SHIFT_FREE);
        auto filteredMeasuredData = DoasFitPreparation::PrepareMeasuredSpectrum(measuredSpectrum, skySpectrum, so2FitWindow.fitType);
        sut.Setup(so2FitWindow);

        // Act
        DoasResult result;
        sut.Run(filteredMeasuredData.data(), filteredMeasuredData.size(), result);

        // Assert
        // Verify the setup
        REQUIRE(result.fitLow == so2FitWindow.fitLow);
        REQUIRE(result.fitHigh == so2FitWindow.fitHigh);
        REQUIRE(result.referenceResult.size() == 3);
        REQUIRE(static_cast<int>(result.polynomialCoefficients.size()) == so2FitWindow.polyOrder + 1);
        REQUIRE(static_cast<int>(result.polynomialValues.size()) == so2FitWindow.fitHigh - so2FitWindow.fitLow);

        // Verify that the result is indeed correct!
        REQUIRE(std::abs(result.chiSquare - 0.01) < 0.005);
        REQUIRE(std::abs(result.delta - 0.04) < 0.02);

        REQUIRE(result.referenceResult[0].name == "SO2");
        REQUIRE(std::abs(result.referenceResult[0].column + 2.2e18) < 1e17); // the SO2 column result is negative here
        REQUIRE(std::abs(result.referenceResult[0].columnError - 1.0e17) < 1e16);

        REQUIRE(result.referenceResult[1].name == "O3");
        REQUIRE(std::abs(result.referenceResult[1].column - 1.1e17) < 6e16);
        REQUIRE(std::abs(result.referenceResult[1].columnError - 2.2e17) < 1e16);

        REQUIRE(result.referenceResult[2].name == "sky");
        REQUIRE(std::abs(result.referenceResult[2].column - 1) < std::numeric_limits<double>::epsilon());
        REQUIRE(std::abs(result.referenceResult[2].columnError) < std::numeric_limits<double>::epsilon());
    }
}
