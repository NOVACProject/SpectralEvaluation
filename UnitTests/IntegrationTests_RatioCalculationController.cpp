#include <SpectralEvaluation/DialogControllers/RatioCalculationController.h>
#include <SpectralEvaluation/File/FitWindowFileHandler.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/File/ScanEvaluationLogFileHandler.h>

#include "catch.hpp"
#include "TestData.h"

using namespace novac;

// region helper methods

RatioCalculationFitSetup GetSetupOfFitWindowsForTest()
{
    RatioCalculationFitSetup result;

    // Read in one fit window. This helps us with ready-to-use references..
    CFitWindowFileHandler fitWindowFileHandler;
    auto allWindows = fitWindowFileHandler.ReadFitWindowFile(TestData::GetBrORatioFitWindowFileSO2());
    REQUIRE(allWindows.size() == 1);
    result.so2Window = allWindows.front();

    allWindows = fitWindowFileHandler.ReadFitWindowFile(TestData::GetBrORatioFitWindowFileBrO());
    REQUIRE(allWindows.size() == 1);
    result.broWindow = allWindows.front();

    return result;
}

// endregion helper methods

TEST_CASE("RatioCalculationController - SetupFitWindows - SO2 not included SO2 window - throws invalid_argument", "[RatioCalculationController][Ratios]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    novac::ConsoleLog log;
    RatioCalculationController sut(log);
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_includeInMajor = false;
    so2Reference->m_includeInMinor = false;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(sut, StandardDoasSpecie::BRO);
    broReference->m_path = testWindows.broWindow.reference[0].m_path;
    broReference->m_includeInMajor = true;
    broReference->m_includeInMinor = true;
    GetReferenceFor(sut, StandardDoasSpecie::RING)->m_automaticallyCalculate = false;
    GetReferenceFor(sut, StandardDoasSpecie::RING_LAMBDA4)->m_automaticallyCalculate = false;


    // Act
    REQUIRE_THROWS(sut.SetupFitWindows());
}

TEST_CASE("RatioCalculationController - SetupFitWindows - BrO not included BrO window - throws invalid_argument", "[RatioCalculationController][Ratios]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    novac::ConsoleLog log;
    RatioCalculationController sut(log);
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = testWindows.so2Window.reference[0].m_path;
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = false;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(sut, StandardDoasSpecie::BRO);
    broReference->m_includeInMajor = true;
    broReference->m_includeInMinor = true;
    GetReferenceFor(sut, StandardDoasSpecie::RING)->m_automaticallyCalculate = false;
    GetReferenceFor(sut, StandardDoasSpecie::RING_LAMBDA4)->m_automaticallyCalculate = false;


    // Act
    REQUIRE_THROWS(sut.SetupFitWindows());
}

TEST_CASE("RatioCalculationController - SetupFitWindows - SO2 and BrO included in both windows", "[RatioCalculationController][Ratios]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    novac::ConsoleLog log;
    RatioCalculationController sut(log);
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = testWindows.so2Window.reference[0].m_path;
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = true;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(sut, StandardDoasSpecie::BRO);
    broReference->m_path = testWindows.broWindow.reference[0].m_path;
    broReference->m_includeInMajor = true;
    broReference->m_includeInMinor = true;
    GetReferenceFor(sut, StandardDoasSpecie::RING)->m_automaticallyCalculate = false;
    GetReferenceFor(sut, StandardDoasSpecie::RING_LAMBDA4)->m_automaticallyCalculate = false;

    sut.m_so2PolynomialOrder = 5;
    sut.m_broPolynomialOrder = 4;

    // Act
    const auto result = sut.SetupFitWindows();

    // Assert
    // SO2 window, notice the order of the references.
    REQUIRE(result->so2Window.name == "SO2");
    REQUIRE(result->so2Window.polyOrder == 5);
    REQUIRE(result->so2Window.fitType == novac::FIT_TYPE::FIT_POLY);
    REQUIRE(result->so2Window.reference.size() == 2);
    REQUIRE(result->so2Window.reference[0].m_path == so2Reference->m_path);
    REQUIRE(result->so2Window.reference[0].m_specieName == "SO2");
    REQUIRE(result->so2Window.reference[0].m_shiftOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->so2Window.reference[0].m_squeezeOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->so2Window.reference[1].m_path == broReference->m_path);
    REQUIRE(result->so2Window.reference[1].m_specieName == "BrO");
    REQUIRE(result->so2Window.reference[1].m_shiftOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->so2Window.reference[1].m_squeezeOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->so2Window.ringCalculation == novac::RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING);


    // BrO window, notice the different order of the references.
    REQUIRE(result->broWindow.name == "BrO");
    REQUIRE(result->broWindow.polyOrder == 4);
    REQUIRE(result->broWindow.fitType == novac::FIT_TYPE::FIT_POLY);
    REQUIRE(result->broWindow.reference.size() == 2);
    REQUIRE(result->broWindow.reference[0].m_path == broReference->m_path);
    REQUIRE(result->broWindow.reference[0].m_specieName == "BrO");
    REQUIRE(result->broWindow.reference[0].m_shiftOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->broWindow.reference[0].m_squeezeOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->broWindow.reference[1].m_path == so2Reference->m_path);
    REQUIRE(result->broWindow.reference[1].m_specieName == "SO2");
    REQUIRE(result->broWindow.reference[1].m_shiftOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->broWindow.reference[1].m_squeezeOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->broWindow.ringCalculation == novac::RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING);
}

TEST_CASE("RatioCalculationController - SetupFitWindows - Sets up fitLow and fitHigh from SO2 and BrO references", "[RatioCalculationController][Ratios]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    novac::ConsoleLog log;
    RatioCalculationController sut(log);
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = testWindows.so2Window.reference[0].m_path;
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = true;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(sut, StandardDoasSpecie::BRO);
    broReference->m_path = testWindows.broWindow.reference[0].m_path;
    broReference->m_includeInMajor = true;
    broReference->m_includeInMinor = true;
    GetReferenceFor(sut, StandardDoasSpecie::RING)->m_automaticallyCalculate = false;
    GetReferenceFor(sut, StandardDoasSpecie::RING_LAMBDA4)->m_automaticallyCalculate = false;


    // Act
    const auto result = sut.SetupFitWindows();

    // Assert
    // SO2 window
    REQUIRE(result->so2Window.name == "SO2");
    REQUIRE(result->so2Window.fitLow == 439);
    REQUIRE(result->so2Window.fitHigh == 592);

    // BrO window
    REQUIRE(result->broWindow.name == "BrO");
    REQUIRE(result->broWindow.fitLow == 641);
    REQUIRE(result->broWindow.fitHigh == 939);
}

TEST_CASE("RatioCalculationController - SetupFitWindows - Sets fitType to HP_DIV of both windows if settings says so", "[RatioCalculationController][Ratios]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    novac::ConsoleLog log;
    RatioCalculationController sut(log);
    sut.m_doasFitType = novac::FIT_TYPE::FIT_HP_DIV;
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = testWindows.so2Window.reference[0].m_path;
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = true;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(sut, StandardDoasSpecie::BRO);
    broReference->m_path = testWindows.broWindow.reference[0].m_path;
    broReference->m_includeInMajor = true;
    broReference->m_includeInMinor = true;
    GetReferenceFor(sut, StandardDoasSpecie::RING)->m_automaticallyCalculate = false;
    GetReferenceFor(sut, StandardDoasSpecie::RING_LAMBDA4)->m_automaticallyCalculate = false;


    // Act
    const auto result = sut.SetupFitWindows();

    // Assert
    // SO2 window
    REQUIRE(result->so2Window.name == "SO2");
    REQUIRE(result->so2Window.fitType == novac::FIT_TYPE::FIT_HP_DIV);

    // BrO window
    REQUIRE(result->broWindow.name == "BrO");
    REQUIRE(result->broWindow.fitType == novac::FIT_TYPE::FIT_HP_DIV);
}

TEST_CASE("RatioCalculationController - SetupFitWindows - SO2 included in only major window", "[RatioCalculationController][Ratios]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    novac::ConsoleLog log;
    RatioCalculationController sut(log);
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = testWindows.so2Window.reference[0].m_path;
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = false;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(sut, StandardDoasSpecie::BRO);
    broReference->m_path = testWindows.broWindow.reference[0].m_path;
    broReference->m_includeInMajor = false;
    broReference->m_includeInMinor = true;
    GetReferenceFor(sut, StandardDoasSpecie::RING)->m_automaticallyCalculate = false;
    GetReferenceFor(sut, StandardDoasSpecie::RING_LAMBDA4)->m_automaticallyCalculate = false;

    // Act
    const auto result = sut.SetupFitWindows();

    // Assert
    // SO2 window
    REQUIRE(result->so2Window.name == "SO2");
    REQUIRE(result->so2Window.polyOrder == 3);
    REQUIRE(result->so2Window.fitType == novac::FIT_TYPE::FIT_POLY);
    REQUIRE(result->so2Window.reference.size() == 1);
    REQUIRE(result->so2Window.reference[0].m_specieName == "SO2");
    REQUIRE(result->so2Window.reference[0].m_shiftOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->so2Window.reference[0].m_squeezeOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->so2Window.ringCalculation == novac::RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING);


    // BrO window
    REQUIRE(result->broWindow.name == "BrO");
    REQUIRE(result->broWindow.polyOrder == 3);
    REQUIRE(result->broWindow.fitType == novac::FIT_TYPE::FIT_POLY);
    REQUIRE(result->broWindow.reference.size() == 1);
    REQUIRE(result->broWindow.reference[0].m_path == broReference->m_path);
    REQUIRE(result->broWindow.reference[0].m_specieName == "BrO");
    REQUIRE(result->broWindow.reference[0].m_shiftOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->broWindow.reference[0].m_squeezeOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->broWindow.ringCalculation == novac::RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING);
}

TEST_CASE("RatioCalculationController - SetupFitWindows - Two references in each window", "[RatioCalculationController][Ratios]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    novac::ConsoleLog log;
    RatioCalculationController sut(log);
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = testWindows.so2Window.reference[0].m_path;
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = false;
    ReferenceForRatioCalculation* o3Reference = GetReferenceFor(sut, StandardDoasSpecie::O3);
    o3Reference->m_path = testWindows.so2Window.reference[1].m_path;
    o3Reference->m_includeInMajor = true;
    o3Reference->m_includeInMinor = true;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(sut, StandardDoasSpecie::BRO);
    broReference->m_path = testWindows.broWindow.reference[1].m_path;
    broReference->m_includeInMajor = false;
    broReference->m_includeInMinor = true;
    GetReferenceFor(sut, StandardDoasSpecie::RING)->m_automaticallyCalculate = false;
    GetReferenceFor(sut, StandardDoasSpecie::RING_LAMBDA4)->m_automaticallyCalculate = false;

    // Act
    const auto result = sut.SetupFitWindows();

    // Assert
    // SO2 window
    REQUIRE(result->so2Window.name == "SO2");
    REQUIRE(result->so2Window.polyOrder == 3);
    REQUIRE(result->so2Window.fitType == novac::FIT_TYPE::FIT_POLY);
    REQUIRE(result->so2Window.reference.size() == 2);
    REQUIRE(result->so2Window.reference[0].m_specieName == "SO2");
    REQUIRE(result->so2Window.reference[0].m_shiftOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->so2Window.reference[0].m_shiftValue == 0.0);
    REQUIRE(result->so2Window.reference[0].m_squeezeOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->so2Window.reference[0].m_squeezeValue == 1.0);
    REQUIRE(result->so2Window.reference[1].m_specieName == "O3");
    REQUIRE(result->so2Window.reference[1].m_shiftOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->so2Window.reference[1].m_squeezeOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->so2Window.ringCalculation == novac::RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING);


    // BrO window
    REQUIRE(result->broWindow.name == "BrO");
    REQUIRE(result->broWindow.polyOrder == 3);
    REQUIRE(result->broWindow.fitType == novac::FIT_TYPE::FIT_POLY);
    REQUIRE(result->broWindow.reference.size() == 2);
    REQUIRE(result->broWindow.reference[0].m_specieName == "BrO");
    REQUIRE(result->broWindow.reference[0].m_shiftOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->broWindow.reference[0].m_shiftValue == 0.0);
    REQUIRE(result->broWindow.reference[0].m_squeezeOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->broWindow.reference[0].m_squeezeValue == 1.0);
    REQUIRE(result->broWindow.reference[1].m_specieName == "O3");
    REQUIRE(result->broWindow.reference[1].m_shiftOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->broWindow.reference[1].m_squeezeOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->broWindow.ringCalculation == novac::RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING);
}

TEST_CASE("RatioCalculationController - SetupFitWindows - Calculated Ring included in both windows", "[RatioCalculationController][Ratios]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    novac::ConsoleLog log;
    RatioCalculationController sut(log);
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = testWindows.so2Window.reference[0].m_path;
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = false;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(sut, StandardDoasSpecie::BRO);
    broReference->m_path = testWindows.broWindow.reference[0].m_path;
    broReference->m_includeInMajor = false;
    broReference->m_includeInMinor = true;
    ReferenceForRatioCalculation* ringReference = GetReferenceFor(sut, StandardDoasSpecie::RING);
    ringReference->m_path = "";
    ringReference->m_automaticallyCalculate = true;
    ringReference->m_includeInMajor = true;
    ringReference->m_includeInMinor = true;
    GetReferenceFor(sut, StandardDoasSpecie::RING_LAMBDA4)->m_automaticallyCalculate = false;


    // Act
    const auto result = sut.SetupFitWindows();

    // Assert
    // SO2 window
    REQUIRE(result->so2Window.name == "SO2");
    REQUIRE(result->so2Window.polyOrder == 3);
    REQUIRE(result->so2Window.fitType == novac::FIT_TYPE::FIT_POLY);
    REQUIRE(result->so2Window.reference.size() == 1);
    REQUIRE(result->so2Window.ringCalculation == novac::RING_CALCULATION_OPTION::CALCULATE_RING);


    // BrO window
    REQUIRE(result->broWindow.name == "BrO");
    REQUIRE(result->broWindow.polyOrder == 3);
    REQUIRE(result->broWindow.fitType == novac::FIT_TYPE::FIT_POLY);
    REQUIRE(result->broWindow.reference.size() == 1);
    REQUIRE(result->broWindow.ringCalculation == novac::RING_CALCULATION_OPTION::CALCULATE_RING);
}

TEST_CASE("RatioCalculationController - SetupFitWindows - Calculated Ring and Ring*Lambda4 included in both windows", "[RatioCalculationController][Ratios]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    novac::ConsoleLog log;
    RatioCalculationController sut(log);
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = testWindows.so2Window.reference[0].m_path;
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = false;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(sut, StandardDoasSpecie::BRO);
    broReference->m_path = testWindows.broWindow.reference[0].m_path;
    broReference->m_includeInMajor = false;
    broReference->m_includeInMinor = true;
    ReferenceForRatioCalculation* ringReference = GetReferenceFor(sut, StandardDoasSpecie::RING);
    ringReference->m_path = "";
    ringReference->m_automaticallyCalculate = true;
    ringReference->m_includeInMajor = true;
    ringReference->m_includeInMinor = true;
    ReferenceForRatioCalculation* ringReference2 = GetReferenceFor(sut, StandardDoasSpecie::RING_LAMBDA4);
    ringReference2->m_path = "";
    ringReference2->m_automaticallyCalculate = true;
    ringReference2->m_includeInMajor = true;
    ringReference2->m_includeInMinor = true;

    // Act
    const auto result = sut.SetupFitWindows();

    // Assert
    // SO2 window
    REQUIRE(result->so2Window.name == "SO2");
    REQUIRE(result->so2Window.polyOrder == 3);
    REQUIRE(result->so2Window.fitType == novac::FIT_TYPE::FIT_POLY);
    REQUIRE(result->so2Window.reference.size() == 1);
    REQUIRE(result->so2Window.ringCalculation == novac::RING_CALCULATION_OPTION::CALCULATE_RING_X2);


    // BrO window
    REQUIRE(result->broWindow.name == "BrO");
    REQUIRE(result->broWindow.polyOrder == 3);
    REQUIRE(result->broWindow.fitType == novac::FIT_TYPE::FIT_POLY);
    REQUIRE(result->broWindow.reference.size() == 1);
    REQUIRE(result->broWindow.ringCalculation == novac::RING_CALCULATION_OPTION::CALCULATE_RING_X2);
}

TEST_CASE("RatioCalculationController - Evaluate", "[RatioCalculationController][Ratios]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    novac::ConsoleLog log;
    novac::LogContext context;
    RatioCalculationController sut(log);
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references for a good fit
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = testWindows.so2Window.reference[0].m_path;
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = false;
    ReferenceForRatioCalculation* o3Reference = GetReferenceFor(sut, StandardDoasSpecie::O3);
    o3Reference->m_path = testWindows.so2Window.reference[1].m_path;
    o3Reference->m_includeInMajor = true;
    o3Reference->m_includeInMinor = true;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(sut, StandardDoasSpecie::BRO);
    broReference->m_path = testWindows.broWindow.reference[0].m_path;
    broReference->m_includeInMajor = false;
    broReference->m_includeInMinor = true;
    ReferenceForRatioCalculation* ringReference = GetReferenceFor(sut, StandardDoasSpecie::RING);
    ringReference->m_automaticallyCalculate = true;
    ringReference->m_includeInMajor = true;
    ringReference->m_includeInMinor = true;
    ReferenceForRatioCalculation* ringReference2 = GetReferenceFor(sut, StandardDoasSpecie::RING_LAMBDA4);
    ringReference2->m_automaticallyCalculate = true;
    ringReference2->m_includeInMajor = true;
    ringReference2->m_includeInMinor = true;

    // Prepare the test by reading in the .pak-file and the evaluation result and calculate the plume-properties from the result.
    CScanFileHandler fileHandler(log);
    const bool scanFileIsOk = fileHandler.CheckScanFile(context, TestData::GetBrORatioScanFile1());
    REQUIRE(scanFileIsOk); // check assumption on the setup

    novac::CScanEvaluationLogFileHandler evaluationFileHandler;
    const bool evaluationFileIsOk = evaluationFileHandler.ReadEvaluationLog(TestData::GetBrORatioEvaluationFile1());
    REQUIRE(evaluationFileIsOk); // check assumption on the setup
    REQUIRE(evaluationFileHandler.m_scan.size() == 1); // check assumption on the setup

    // setup the fit windows
    auto fitWindowsetup = sut.SetupFitWindows();

    // Act
    const auto result = sut.EvaluateScan(context, fileHandler, evaluationFileHandler.m_scan[0], fitWindowsetup);

    // Assert
    // Verify that the result is indeed correct!
    REQUIRE(result.ratio.minorSpecieName == "BrO");
    REQUIRE(result.ratio.majorSpecieName == "SO2");

    REQUIRE(result.debugInfo.doasResults[0].referenceResult[0].column == Approx(2.0e18).margin(1e17));
    REQUIRE(result.debugInfo.doasResults[1].referenceResult[0].column == Approx(1.5e14).margin(2e13));

    REQUIRE(result.ratio.ratio == Approx(7.4e-5).margin(1e-6));
    REQUIRE(result.ratio.error == Approx(1.4e-5).margin(1e-6));

    // Also, verify that the result contains the additional informatino which we need in order to show the results to the user
    REQUIRE(result.filename == fileHandler.GetFileName());
    REQUIRE(result.debugInfo.errorMessage == "");
    REQUIRE(result.debugInfo.plumeSpectrumIndices.size() == 10);
    REQUIRE(result.debugInfo.plumeSpectrumIndices.front() == 14);
    REQUIRE(result.debugInfo.plumeSpectrumIndices.back() == 23);
    REQUIRE(result.debugInfo.outOfPlumeSpectrumIndices.size() == 10);
    REQUIRE(result.debugInfo.outOfPlumeSpectrumIndices.front() == 29);
    REQUIRE(result.debugInfo.outOfPlumeSpectrumIndices.back() == 50);
    REQUIRE(result.debugInfo.inPlumeSpectrum.NumSpectra() == 15 * static_cast<long>(result.debugInfo.plumeSpectrumIndices.size()));
    REQUIRE(result.debugInfo.outOfPlumeSpectrum.NumSpectra() == 15 * static_cast<long>(result.debugInfo.outOfPlumeSpectrumIndices.size()));

    // The DOAS results
    REQUIRE(result.debugInfo.doasResults.size() == 2);
    REQUIRE(result.debugInfo.doasResults[0].referenceResult.size() == 6); // SO2, O3, Ring, RingxLamda4, IntensityPoly, Sky
    REQUIRE(result.debugInfo.doasResults[0].fitLow == fitWindowsetup->so2Window.fitLow);
    REQUIRE(result.debugInfo.doasResults[0].fitHigh == fitWindowsetup->so2Window.fitHigh);
    REQUIRE(result.debugInfo.doasResults[1].referenceResult.size() == 6); // BrO, O3, Ring, RingxLamda4, IntensityPoly, Sky
    REQUIRE(result.debugInfo.doasResults[1].fitLow == fitWindowsetup->broWindow.fitLow);
    REQUIRE(result.debugInfo.doasResults[1].fitHigh == fitWindowsetup->broWindow.fitHigh);
}

TEST_CASE("RatioCalculationController - Evaluate - plume too wide", "[RatioCalculationController][Ratios]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    novac::ConsoleLog log;
    novac::LogContext context;
    RatioCalculationController sut(log);
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references for a good fit
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = testWindows.so2Window.reference[0].m_path;
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = false;
    ReferenceForRatioCalculation* o3Reference = GetReferenceFor(sut, StandardDoasSpecie::O3);
    o3Reference->m_path = testWindows.so2Window.reference[1].m_path;
    o3Reference->m_includeInMajor = true;
    o3Reference->m_includeInMinor = true;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(sut, StandardDoasSpecie::BRO);
    broReference->m_path = testWindows.broWindow.reference[0].m_path;
    broReference->m_includeInMajor = false;
    broReference->m_includeInMinor = true;
    ReferenceForRatioCalculation* ringReference = GetReferenceFor(sut, StandardDoasSpecie::RING);
    ringReference->m_automaticallyCalculate = true;
    ringReference->m_includeInMajor = true;
    ringReference->m_includeInMinor = true;
    ReferenceForRatioCalculation* ringReference2 = GetReferenceFor(sut, StandardDoasSpecie::RING_LAMBDA4);
    ringReference2->m_automaticallyCalculate = true;
    ringReference2->m_includeInMajor = true;
    ringReference2->m_includeInMinor = true;

    // Prepare the test by reading in the .pak-file and the evaluation result and calculate the plume-properties from the result.
    CScanFileHandler fileHandler(log);
    const bool scanFileIsOk = fileHandler.CheckScanFile(context, TestData::GetBrORatioScanFile2());
    REQUIRE(scanFileIsOk); // check assumption on the setup

    novac::CScanEvaluationLogFileHandler evaluationFileHandler;
    const bool evaluationFileIsOk = evaluationFileHandler.ReadEvaluationLog(TestData::GetBrORatioEvaluationFile2());
    REQUIRE(evaluationFileIsOk); // check assumption on the setup
    REQUIRE(evaluationFileHandler.m_scan.size() == 1); // check assumption on the setup

    // setup the fit windows
    auto fitWindowsetup = sut.SetupFitWindows();

    // Act
    const auto result = sut.EvaluateScan(context, fileHandler, evaluationFileHandler.m_scan[0], fitWindowsetup);

    // Assert
    REQUIRE(!result.debugInfo.errorMessage.empty());

    REQUIRE(std::abs(result.ratio.ratio) < std::numeric_limits<double>::epsilon());
}

TEST_CASE("RatioCalculationController - Evaluate without Ring", "[RatioCalculationController][Ratios]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    novac::ConsoleLog log;
    novac::LogContext context;
    RatioCalculationController sut(log);
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references for a good fit
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = testWindows.so2Window.reference[0].m_path;
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = false;
    ReferenceForRatioCalculation* o3Reference = GetReferenceFor(sut, StandardDoasSpecie::O3);
    o3Reference->m_path = testWindows.so2Window.reference[1].m_path;
    o3Reference->m_includeInMajor = true;
    o3Reference->m_includeInMinor = true;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(sut, StandardDoasSpecie::BRO);
    broReference->m_path = testWindows.broWindow.reference[0].m_path;
    broReference->m_includeInMajor = false;
    broReference->m_includeInMinor = true;
    ReferenceForRatioCalculation* ringReference = GetReferenceFor(sut, StandardDoasSpecie::RING);
    ringReference->m_automaticallyCalculate = false;
    ReferenceForRatioCalculation* ringReference2 = GetReferenceFor(sut, StandardDoasSpecie::RING_LAMBDA4);
    ringReference2->m_automaticallyCalculate = false;

    // Prepare the test by reading in the .pak-file and the evaluation result and calculate the plume-properties from the result.
    CScanFileHandler fileHandler(log);
    const bool scanFileIsOk = fileHandler.CheckScanFile(context, TestData::GetBrORatioScanFile1());
    REQUIRE(scanFileIsOk); // check assumption on the setup

    novac::CScanEvaluationLogFileHandler evaluationFileHandler;
    const bool evaluationFileIsOk = evaluationFileHandler.ReadEvaluationLog(TestData::GetBrORatioEvaluationFile1());
    REQUIRE(evaluationFileIsOk); // check assumption on the setup
    REQUIRE(evaluationFileHandler.m_scan.size() == 1); // check assumption on the setup

    // setup the fit windows
    auto fitWindowsetup = sut.SetupFitWindows();

    // Act
    const auto result = sut.EvaluateScan(context, fileHandler, evaluationFileHandler.m_scan[0], fitWindowsetup);

    // Assert
    // Verify that the result is indeed correct!
    REQUIRE(result.ratio.minorSpecieName == "BrO");
    REQUIRE(result.ratio.majorSpecieName == "SO2");
    REQUIRE(std::abs(result.ratio.ratio) < 1.2e-4); // During development there may be variations in the result. Just verify the range is ok for now.
    REQUIRE(std::abs(result.ratio.ratio) > 6e-5); // During development there may be variations in the result. Just verify the range is ok for now.
    REQUIRE(std::abs(result.ratio.error) < std::abs(result.ratio.ratio / 2.0));

    // Also, verify that the result contains the additional informatino which we need in order to show the results to the user
    REQUIRE(result.filename == fileHandler.GetFileName());
    REQUIRE(result.debugInfo.errorMessage == "");
    REQUIRE(static_cast<int>(result.debugInfo.plumeSpectrumIndices.size()) >= sut.m_ratioEvaluationSettings.minNumberOfSpectraInPlume);
    REQUIRE(static_cast<int>(result.debugInfo.outOfPlumeSpectrumIndices.size()) >= sut.m_ratioEvaluationSettings.numberOfSpectraOutsideOfPlume);
    REQUIRE(result.debugInfo.inPlumeSpectrum.NumSpectra() == 15 * static_cast<long>(result.debugInfo.plumeSpectrumIndices.size()));
    REQUIRE(result.debugInfo.outOfPlumeSpectrum.NumSpectra() == 15 * static_cast<long>(result.debugInfo.outOfPlumeSpectrumIndices.size()));

    // The DOAS results
    REQUIRE(result.debugInfo.doasResults.size() == 2);
    REQUIRE(result.debugInfo.doasResults[0].referenceResult.size() == 4); // SO2, O3, IntensityPoly, Sky
    REQUIRE(result.debugInfo.doasResults[0].referenceResult[0].name == "SO2");
    REQUIRE(result.debugInfo.doasResults[0].referenceResult[1].name == "O3");
    REQUIRE(result.debugInfo.doasResults[0].referenceResult[2].name == "sky");
    REQUIRE(result.debugInfo.doasResults[0].referenceResult[3].name == "offset");
    REQUIRE(result.debugInfo.doasResults[0].fitLow == fitWindowsetup->so2Window.fitLow);
    REQUIRE(result.debugInfo.doasResults[0].fitHigh == fitWindowsetup->so2Window.fitHigh);
    REQUIRE(result.debugInfo.doasResults[1].referenceResult.size() == 4); // BrO, O3, IntensityPoly, Sky
    REQUIRE(result.debugInfo.doasResults[1].referenceResult[0].name == "BrO");
    REQUIRE(result.debugInfo.doasResults[1].referenceResult[1].name == "O3");
    REQUIRE(result.debugInfo.doasResults[1].referenceResult[2].name == "sky");
    REQUIRE(result.debugInfo.doasResults[1].referenceResult[3].name == "offset");
    REQUIRE(result.debugInfo.doasResults[1].fitLow == fitWindowsetup->broWindow.fitLow);
    REQUIRE(result.debugInfo.doasResults[1].fitHigh == fitWindowsetup->broWindow.fitHigh);
}

TEST_CASE("RatioCalculationController - Evaluate without offset polynomial", "[RatioCalculationController][Ratios]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    novac::ConsoleLog log;
    novac::LogContext context;
    RatioCalculationController sut(log);
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references, with varied number of references in each window
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = testWindows.so2Window.reference[0].m_path;
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = false;
    ReferenceForRatioCalculation* o3Reference = GetReferenceFor(sut, StandardDoasSpecie::O3);
    o3Reference->m_path = testWindows.so2Window.reference[1].m_path;
    o3Reference->m_includeInMajor = true;
    o3Reference->m_includeInMinor = true;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(sut, StandardDoasSpecie::BRO);
    broReference->m_path = testWindows.broWindow.reference[0].m_path;
    broReference->m_includeInMajor = false;
    broReference->m_includeInMinor = true;
    ReferenceForRatioCalculation* ringReference = GetReferenceFor(sut, StandardDoasSpecie::RING);
    ringReference->m_automaticallyCalculate = true;
    ringReference->m_includeInMajor = true; // Notice that it is not possible to include only Ring_Lambda4 without 'Ring'...
    ringReference->m_includeInMinor = true;
    ReferenceForRatioCalculation* ringReference2 = GetReferenceFor(sut, StandardDoasSpecie::RING_LAMBDA4);
    ringReference2->m_automaticallyCalculate = true;
    ringReference2->m_includeInMajor = true;
    ringReference2->m_includeInMinor = false;

    // Prepare the test by reading in the .pak-file and the evaluation result and calculate the plume-properties from the result.
    CScanFileHandler fileHandler(log);
    const bool scanFileIsOk = fileHandler.CheckScanFile(context, TestData::GetBrORatioScanFile1());
    REQUIRE(scanFileIsOk); // check assumption on the setup

    novac::CScanEvaluationLogFileHandler evaluationFileHandler;
    const bool evaluationFileIsOk = evaluationFileHandler.ReadEvaluationLog(TestData::GetBrORatioEvaluationFile1());
    REQUIRE(evaluationFileIsOk); // check assumption on the setup
    REQUIRE(evaluationFileHandler.m_scan.size() == 1); // check assumption on the setup

    // setup the fit windows
    auto fitWindowsetup = sut.SetupFitWindows();

    // Remove the offset poly in the two fit windows
    fitWindowsetup->so2Window.includeIntensitySpacePolyominal = false;
    fitWindowsetup->broWindow.includeIntensitySpacePolyominal = false;


    // Act
    const auto result = sut.EvaluateScan(context, fileHandler, evaluationFileHandler.m_scan[0], fitWindowsetup);

    // Assert
    // Verify that the result is indeed correct!
    REQUIRE(result.ratio.minorSpecieName == "BrO");
    REQUIRE(result.ratio.majorSpecieName == "SO2");
    REQUIRE(std::abs(result.ratio.ratio) < 1.2e-4); // During development there may be variations in the result. Just verify the range is ok for now.
    REQUIRE(std::abs(result.ratio.ratio) > 6e-5); // During development there may be variations in the result. Just verify the range is ok for now.
    REQUIRE(std::abs(result.ratio.error) < std::abs(result.ratio.ratio / 2.0));

    // Also, verify that the result contains the additional informatino which we need in order to show the results to the user
    REQUIRE(result.filename == fileHandler.GetFileName());
    REQUIRE(result.debugInfo.errorMessage == "");
    REQUIRE(static_cast<int>(result.debugInfo.plumeSpectrumIndices.size()) >= sut.m_ratioEvaluationSettings.minNumberOfSpectraInPlume);
    REQUIRE(static_cast<int>(result.debugInfo.outOfPlumeSpectrumIndices.size()) >= sut.m_ratioEvaluationSettings.numberOfSpectraOutsideOfPlume);
    REQUIRE(result.debugInfo.inPlumeSpectrum.NumSpectra() == 15 * static_cast<long>(result.debugInfo.plumeSpectrumIndices.size()));
    REQUIRE(result.debugInfo.outOfPlumeSpectrum.NumSpectra() == 15 * static_cast<long>(result.debugInfo.outOfPlumeSpectrumIndices.size()));

    // The DOAS results
    REQUIRE(result.debugInfo.doasResults.size() == 2);
    REQUIRE(result.debugInfo.doasResults[0].referenceResult.size() == 5); // SO2, O3, IntensityPoly, Sky
    REQUIRE(result.debugInfo.doasResults[0].referenceResult[0].name == "SO2");
    REQUIRE(result.debugInfo.doasResults[0].referenceResult[1].name == "O3");
    REQUIRE(result.debugInfo.doasResults[0].referenceResult[2].name == "sky");
    REQUIRE(result.debugInfo.doasResults[0].referenceResult[3].name == "ring");
    REQUIRE(result.debugInfo.doasResults[0].referenceResult[4].name == "ringLambda4");
    REQUIRE(result.debugInfo.doasResults[0].fitLow == fitWindowsetup->so2Window.fitLow);
    REQUIRE(result.debugInfo.doasResults[0].fitHigh == fitWindowsetup->so2Window.fitHigh);
    REQUIRE(result.debugInfo.doasResults[1].referenceResult.size() == 4); // BrO, O3, IntensityPoly, Sky
    REQUIRE(result.debugInfo.doasResults[1].referenceResult[0].name == "BrO");
    REQUIRE(result.debugInfo.doasResults[1].referenceResult[1].name == "O3");
    REQUIRE(result.debugInfo.doasResults[1].referenceResult[2].name == "sky");
    REQUIRE(result.debugInfo.doasResults[1].referenceResult[3].name == "ring");
    REQUIRE(result.debugInfo.doasResults[1].fitLow == fitWindowsetup->broWindow.fitLow);
    REQUIRE(result.debugInfo.doasResults[1].fitHigh == fitWindowsetup->broWindow.fitHigh);
}

TEST_CASE("RatioCalculationController - DoInitialEvaluation", "[RatioCalculationController][Ratios]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    novac::ConsoleLog log;
    novac::LogContext context;
    RatioCalculationController sut(log);
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references for a good fit
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = testWindows.so2Window.reference[0].m_path;
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = false;
    ReferenceForRatioCalculation* o3Reference = GetReferenceFor(sut, StandardDoasSpecie::O3);
    o3Reference->m_path = testWindows.so2Window.reference[1].m_path;
    o3Reference->m_includeInMajor = true;
    o3Reference->m_includeInMinor = true;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(sut, StandardDoasSpecie::BRO);
    broReference->m_path = testWindows.broWindow.reference[0].m_path;
    broReference->m_includeInMajor = false;
    broReference->m_includeInMinor = true;
    ReferenceForRatioCalculation* ringReference = GetReferenceFor(sut, StandardDoasSpecie::RING);
    ringReference->m_automaticallyCalculate = true;
    ringReference->m_includeInMajor = true;
    ringReference->m_includeInMinor = true;
    ReferenceForRatioCalculation* ringReference2 = GetReferenceFor(sut, StandardDoasSpecie::RING_LAMBDA4);
    ringReference2->m_automaticallyCalculate = true;
    ringReference2->m_includeInMajor = true;
    ringReference2->m_includeInMinor = true;

    // Prepare the test by reading in the .pak-file and the evaluation result and calculate the plume-properties from the result.
    CScanFileHandler fileHandler(log);
    const bool scanFileIsOk = fileHandler.CheckScanFile(context, TestData::GetBrORatioScanFile1());
    REQUIRE(scanFileIsOk); // check assumption on the setup

    novac::CScanEvaluationLogFileHandler evaluationFileHandler;
    const bool evaluationFileIsOk = evaluationFileHandler.ReadEvaluationLog(TestData::GetBrORatioEvaluationFile1());
    REQUIRE(evaluationFileIsOk); // check assumption on the setup
    REQUIRE(evaluationFileHandler.m_scan.size() == 1); // check assumption on the setup

    // setup the fit windows
    auto fitWindowsetup = sut.SetupFitWindows();

    // Act
    const auto result = sut.DoInitialEvaluation(fileHandler, fitWindowsetup);

    // Assert
    REQUIRE(result.m_spec.size() == evaluationFileHandler.m_scan[0].m_spec.size());
    for (size_t idx = 0; idx < evaluationFileHandler.m_scan[0].m_spec.size(); ++idx)
    {
        // With the above setup, the DOAS fit should be better than the original and hence should the chi2 be smaller.
        REQUIRE(result.m_spec[idx].m_chiSquare < evaluationFileHandler.m_scan[0].m_spec[idx].m_chiSquare);

        REQUIRE(result.m_spec[idx].m_referenceResult.size() == 3);
    }
}

TEST_CASE("RatioCalculationController - Save and LoadSetup restores original data", "[RatioCalculationController][Ratios]")
{
    novac::ConsoleLog log;
    RatioCalculationController original(log);
    // Setup the references
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(original, StandardDoasSpecie::SO2);
    so2Reference->m_path = "This is a test ";
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = false;
    ReferenceForRatioCalculation* o3Reference = GetReferenceFor(original, StandardDoasSpecie::O3);
    o3Reference->m_path = "to verify that ";
    o3Reference->m_includeInMajor = true;
    o3Reference->m_includeInMinor = true;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(original, StandardDoasSpecie::BRO);
    broReference->m_path = "save and load does work";
    broReference->m_includeInMajor = false;
    broReference->m_includeInMinor = true;
    ReferenceForRatioCalculation* ringReference = GetReferenceFor(original, StandardDoasSpecie::RING);
    ringReference->m_automaticallyCalculate = true;
    ringReference->m_includeInMajor = true;
    ringReference->m_includeInMinor = true;
    ReferenceForRatioCalculation* ringReference2 = GetReferenceFor(original, StandardDoasSpecie::RING_LAMBDA4);
    ringReference2->m_automaticallyCalculate = true;
    ringReference2->m_includeInMajor = true;
    ringReference2->m_includeInMinor = true;
    original.m_so2FitRange = novac::WavelengthRange(305.1, 311.2);
    original.m_so2PolynomialOrder = 7;
    original.m_broFitRange = novac::WavelengthRange(333.3, 366.6);
    original.m_broPolynomialOrder = 9;
    original.m_doasFitType = novac::FIT_TYPE::FIT_HP_SUB;
    RatioCalculationController restored(log);
    const std::string temporaryFile = TestData::GetTemporaryConfigurationFileName();

    // Act
    original.SaveSetup(temporaryFile);
    restored.LoadSetup(temporaryFile);

    // Assert the setup of the two objects is identical
    REQUIRE(std::abs(original.m_so2FitRange.low - restored.m_so2FitRange.low) < 0.01);
    REQUIRE(std::abs(original.m_so2FitRange.high - restored.m_so2FitRange.high) < 0.01);
    REQUIRE(std::abs(original.m_broFitRange.low - restored.m_broFitRange.low) < 0.01);
    REQUIRE(std::abs(original.m_broFitRange.high - restored.m_broFitRange.high) < 0.01);
    REQUIRE(original.m_so2PolynomialOrder == restored.m_so2PolynomialOrder);
    REQUIRE(original.m_broPolynomialOrder == restored.m_broPolynomialOrder);
    REQUIRE(original.m_doasFitType == restored.m_doasFitType);

    for (size_t ii = 0; ii < original.m_references.size(); ++ii)
    {
        REQUIRE(original.m_references[ii].specie == restored.m_references[ii].specie);
        REQUIRE(original.m_references[ii].m_name == restored.m_references[ii].m_name);
        REQUIRE(original.m_references[ii].m_path == restored.m_references[ii].m_path);
        REQUIRE(original.m_references[ii].m_includeInMajor == restored.m_references[ii].m_includeInMajor);
        REQUIRE(original.m_references[ii].m_includeInMinor == restored.m_references[ii].m_includeInMinor);
        REQUIRE(original.m_references[ii].m_automaticallyCalculate == restored.m_references[ii].m_automaticallyCalculate);
    }
}