#include <SpectralEvaluation/DialogControllers/RatioCalculationController.h>
#include <SpectralEvaluation/File/FitWindowFileHandler.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/File/ScanEvaluationLogFileHandler.h>

#include "catch.hpp"
#include "TestData.h"

using namespace novac;

#pragma region helper methods

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

#pragma endregion

TEST_CASE("RatioCalculationController - SetupFitWindows - SO2 not included SO2 window - throws invalid_argument", "[RatioCalculationController]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    RatioCalculationController sut;
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_includeInMajor = false;
    so2Reference->m_includeInMinor = false;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(sut, StandardDoasSpecie::BRO);
    broReference->m_path = testWindows.broWindow.ref[0].m_path;
    broReference->m_includeInMajor = true;
    broReference->m_includeInMinor = true;
    GetReferenceFor(sut, StandardDoasSpecie::RING)->m_automaticallyCalculate = false;
    GetReferenceFor(sut, StandardDoasSpecie::RING_LAMBDA4)->m_automaticallyCalculate = false;


    // Act
    REQUIRE_THROWS(sut.SetupFitWindows());
}

TEST_CASE("RatioCalculationController - SetupFitWindows - BrO not included BrO window - throws invalid_argument", "[RatioCalculationController]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    RatioCalculationController sut;
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = testWindows.so2Window.ref[0].m_path;
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

TEST_CASE("RatioCalculationController - SetupFitWindows - SO2 and BrO included in both windows", "[RatioCalculationController]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    RatioCalculationController sut;
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = testWindows.so2Window.ref[0].m_path;
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = true;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(sut, StandardDoasSpecie::BRO);
    broReference->m_path = testWindows.broWindow.ref[0].m_path;
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
    REQUIRE(result->so2Window.nRef == 2);
    REQUIRE(result->so2Window.ref[0].m_path == so2Reference->m_path);
    REQUIRE(result->so2Window.ref[0].m_specieName == "SO2");
    REQUIRE(result->so2Window.ref[0].m_shiftOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->so2Window.ref[0].m_squeezeOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->so2Window.ref[1].m_path == broReference->m_path);
    REQUIRE(result->so2Window.ref[1].m_specieName == "BrO");
    REQUIRE(result->so2Window.ref[1].m_shiftOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->so2Window.ref[1].m_squeezeOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->so2Window.ringCalculation == novac::RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING);


    // BrO window, notice the different order of the references.
    REQUIRE(result->broWindow.name == "BrO");
    REQUIRE(result->broWindow.polyOrder == 4);
    REQUIRE(result->broWindow.fitType == novac::FIT_TYPE::FIT_POLY);
    REQUIRE(result->broWindow.nRef == 2);
    REQUIRE(result->broWindow.ref[0].m_path == broReference->m_path);
    REQUIRE(result->broWindow.ref[0].m_specieName == "BrO");
    REQUIRE(result->broWindow.ref[0].m_shiftOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->broWindow.ref[0].m_squeezeOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->broWindow.ref[1].m_path == so2Reference->m_path);
    REQUIRE(result->broWindow.ref[1].m_specieName == "SO2");
    REQUIRE(result->broWindow.ref[1].m_shiftOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->broWindow.ref[1].m_squeezeOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->broWindow.ringCalculation == novac::RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING);
}


TEST_CASE("RatioCalculationController - SetupFitWindows - Sets up fitLow and fitHigh from SO2 and BrO references", "[RatioCalculationController]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    RatioCalculationController sut;
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = testWindows.so2Window.ref[0].m_path;
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = true;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(sut, StandardDoasSpecie::BRO);
    broReference->m_path = testWindows.broWindow.ref[0].m_path;
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

TEST_CASE("RatioCalculationController - SetupFitWindows - SO2 included in only major window", "[RatioCalculationController]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    RatioCalculationController sut;
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = testWindows.so2Window.ref[0].m_path;
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = false;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(sut, StandardDoasSpecie::BRO);
    broReference->m_path = testWindows.broWindow.ref[0].m_path;
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
    REQUIRE(result->so2Window.nRef == 1);
    REQUIRE(result->so2Window.ref[0].m_specieName == "SO2");
    REQUIRE(result->so2Window.ref[0].m_shiftOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->so2Window.ref[0].m_squeezeOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->so2Window.ringCalculation == novac::RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING);


    // BrO window
    REQUIRE(result->broWindow.name == "BrO");
    REQUIRE(result->broWindow.polyOrder == 3);
    REQUIRE(result->broWindow.fitType == novac::FIT_TYPE::FIT_POLY);
    REQUIRE(result->broWindow.nRef == 1);
    REQUIRE(result->broWindow.ref[0].m_path == broReference->m_path);
    REQUIRE(result->broWindow.ref[0].m_specieName == "BrO");
    REQUIRE(result->broWindow.ref[0].m_shiftOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->broWindow.ref[0].m_squeezeOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->broWindow.ringCalculation == novac::RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING);
}

TEST_CASE("RatioCalculationController - SetupFitWindows - Two references in each window", "[RatioCalculationController]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    RatioCalculationController sut;
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = testWindows.so2Window.ref[0].m_path;
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = false;
    ReferenceForRatioCalculation* o3Reference = GetReferenceFor(sut, StandardDoasSpecie::O3);
    o3Reference->m_path = testWindows.so2Window.ref[1].m_path;
    o3Reference->m_includeInMajor = true;
    o3Reference->m_includeInMinor = true;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(sut, StandardDoasSpecie::BRO);
    broReference->m_path = testWindows.broWindow.ref[1].m_path;
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
    REQUIRE(result->so2Window.nRef == 2);
    REQUIRE(result->so2Window.ref[0].m_specieName == "SO2");
    REQUIRE(result->so2Window.ref[0].m_shiftOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->so2Window.ref[0].m_shiftValue == 0.0);
    REQUIRE(result->so2Window.ref[0].m_squeezeOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->so2Window.ref[0].m_squeezeValue == 1.0);
    REQUIRE(result->so2Window.ref[1].m_specieName == "O3");
    REQUIRE(result->so2Window.ref[1].m_shiftOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->so2Window.ref[1].m_squeezeOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->so2Window.ringCalculation == novac::RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING);


    // BrO window
    REQUIRE(result->broWindow.name == "BrO");
    REQUIRE(result->broWindow.polyOrder == 3);
    REQUIRE(result->broWindow.fitType == novac::FIT_TYPE::FIT_POLY);
    REQUIRE(result->broWindow.nRef == 2);
    REQUIRE(result->broWindow.ref[0].m_specieName == "BrO");
    REQUIRE(result->broWindow.ref[0].m_shiftOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->broWindow.ref[0].m_shiftValue == 0.0);
    REQUIRE(result->broWindow.ref[0].m_squeezeOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->broWindow.ref[0].m_squeezeValue == 1.0);
    REQUIRE(result->broWindow.ref[1].m_specieName == "O3");
    REQUIRE(result->broWindow.ref[1].m_shiftOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->broWindow.ref[1].m_squeezeOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->broWindow.ringCalculation == novac::RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING);
}

TEST_CASE("RatioCalculationController - SetupFitWindows - Calculated Ring included in both windows", "[RatioCalculationController]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    RatioCalculationController sut;
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = testWindows.so2Window.ref[0].m_path;
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = false;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(sut, StandardDoasSpecie::BRO);
    broReference->m_path = testWindows.broWindow.ref[0].m_path;
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
    REQUIRE(result->so2Window.nRef == 1);
    REQUIRE(result->so2Window.ringCalculation == novac::RING_CALCULATION_OPTION::CALCULATE_RING);


    // BrO window
    REQUIRE(result->broWindow.name == "BrO");
    REQUIRE(result->broWindow.polyOrder == 3);
    REQUIRE(result->broWindow.fitType == novac::FIT_TYPE::FIT_POLY);
    REQUIRE(result->broWindow.nRef == 1);
    REQUIRE(result->broWindow.ringCalculation == novac::RING_CALCULATION_OPTION::CALCULATE_RING);
}

TEST_CASE("RatioCalculationController - SetupFitWindows - Calculated Ring and Ring*Lambda4 included in both windows", "[RatioCalculationController]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    RatioCalculationController sut;
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = testWindows.so2Window.ref[0].m_path;
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = false;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(sut, StandardDoasSpecie::BRO);
    broReference->m_path = testWindows.broWindow.ref[0].m_path;
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
    REQUIRE(result->so2Window.nRef == 1);
    REQUIRE(result->so2Window.ringCalculation == novac::RING_CALCULATION_OPTION::CALCULATE_RING_X2);


    // BrO window
    REQUIRE(result->broWindow.name == "BrO");
    REQUIRE(result->broWindow.polyOrder == 3);
    REQUIRE(result->broWindow.fitType == novac::FIT_TYPE::FIT_POLY);
    REQUIRE(result->broWindow.nRef == 1);
    REQUIRE(result->broWindow.ringCalculation == novac::RING_CALCULATION_OPTION::CALCULATE_RING_X2);
}

TEST_CASE("RatioCalculationController - Evaluate", "[RatioCalculationController]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    RatioCalculationController sut;
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references for a good fit
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = testWindows.so2Window.ref[0].m_path;
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = false;
    ReferenceForRatioCalculation* o3Reference = GetReferenceFor(sut, StandardDoasSpecie::O3);
    o3Reference->m_path = testWindows.so2Window.ref[1].m_path;
    o3Reference->m_includeInMajor = true;
    o3Reference->m_includeInMinor = true;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(sut, StandardDoasSpecie::BRO);
    broReference->m_path = testWindows.broWindow.ref[0].m_path;
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
    novac::CScanFileHandler fileHandler;
    const bool scanFileIsOk = fileHandler.CheckScanFile(TestData::GetBrORatioScanFile1());
    REQUIRE(scanFileIsOk); // check assumption on the setup

    novac::CScanEvaluationLogFileHandler evaluationFileHandler;
    const bool evaluationFileIsOk = evaluationFileHandler.ReadEvaluationLog(TestData::GetBrORatioEvaluationFile1());
    REQUIRE(evaluationFileIsOk); // check assumption on the setup
    REQUIRE(evaluationFileHandler.m_scan.size() == 1); // check assumption on the setup

    // setup the fit windows
    auto fitWindowsetup = sut.SetupFitWindows();

    // Act
    const auto result = sut.EvaluateScan(fileHandler, evaluationFileHandler.m_scan[0], fitWindowsetup);

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
    REQUIRE(result.debugInfo.plumeSpectra.size() >= sut.m_ratioEvaluationSettings.minNumberOfSpectraInPlume);
    REQUIRE(result.debugInfo.outOfPlumeSpectra.size() >= sut.m_ratioEvaluationSettings.minNumberOfReferenceSpectra);
    REQUIRE(result.debugInfo.inPlumeSpectrum.NumSpectra() == 15 * result.debugInfo.plumeSpectra.size());
    REQUIRE(result.debugInfo.outOfPlumeSpectrum.NumSpectra() == 15 * result.debugInfo.outOfPlumeSpectra.size());

    // The DOAS results
    REQUIRE(result.debugInfo.doasResults.size() == 2);
    REQUIRE(result.debugInfo.doasResults[0].referenceResult.size() == 6); // SO2, O3, Ring, RingxLamda4, IntensityPoly, Sky
    REQUIRE(result.debugInfo.doasResults[0].fitLow == fitWindowsetup->so2Window.fitLow);
    REQUIRE(result.debugInfo.doasResults[0].fitHigh == fitWindowsetup->so2Window.fitHigh);
    REQUIRE(result.debugInfo.doasResults[1].referenceResult.size() == 6); // BrO, O3, Ring, RingxLamda4, IntensityPoly, Sky
    REQUIRE(result.debugInfo.doasResults[1].fitLow == fitWindowsetup->broWindow.fitLow);
    REQUIRE(result.debugInfo.doasResults[1].fitHigh == fitWindowsetup->broWindow.fitHigh);
}

TEST_CASE("RatioCalculationController - Evaluate without Ring", "[RatioCalculationController]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    RatioCalculationController sut;
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references for a good fit
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = testWindows.so2Window.ref[0].m_path;
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = false;
    ReferenceForRatioCalculation* o3Reference = GetReferenceFor(sut, StandardDoasSpecie::O3);
    o3Reference->m_path = testWindows.so2Window.ref[1].m_path;
    o3Reference->m_includeInMajor = true;
    o3Reference->m_includeInMinor = true;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(sut, StandardDoasSpecie::BRO);
    broReference->m_path = testWindows.broWindow.ref[0].m_path;
    broReference->m_includeInMajor = false;
    broReference->m_includeInMinor = true;
    ReferenceForRatioCalculation* ringReference = GetReferenceFor(sut, StandardDoasSpecie::RING);
    ringReference->m_automaticallyCalculate = false;
    ReferenceForRatioCalculation* ringReference2 = GetReferenceFor(sut, StandardDoasSpecie::RING_LAMBDA4);
    ringReference2->m_automaticallyCalculate = false;

    // Prepare the test by reading in the .pak-file and the evaluation result and calculate the plume-properties from the result.
    novac::CScanFileHandler fileHandler;
    const bool scanFileIsOk = fileHandler.CheckScanFile(TestData::GetBrORatioScanFile1());
    REQUIRE(scanFileIsOk); // check assumption on the setup

    novac::CScanEvaluationLogFileHandler evaluationFileHandler;
    const bool evaluationFileIsOk = evaluationFileHandler.ReadEvaluationLog(TestData::GetBrORatioEvaluationFile1());
    REQUIRE(evaluationFileIsOk); // check assumption on the setup
    REQUIRE(evaluationFileHandler.m_scan.size() == 1); // check assumption on the setup

    // setup the fit windows
    auto fitWindowsetup = sut.SetupFitWindows();

    // Act
    const auto result = sut.EvaluateScan(fileHandler, evaluationFileHandler.m_scan[0], fitWindowsetup);

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
    REQUIRE(result.debugInfo.plumeSpectra.size() >= sut.m_ratioEvaluationSettings.minNumberOfSpectraInPlume);
    REQUIRE(result.debugInfo.outOfPlumeSpectra.size() >= sut.m_ratioEvaluationSettings.minNumberOfReferenceSpectra);
    REQUIRE(result.debugInfo.inPlumeSpectrum.NumSpectra() == 15 * result.debugInfo.plumeSpectra.size());
    REQUIRE(result.debugInfo.outOfPlumeSpectrum.NumSpectra() == 15 * result.debugInfo.outOfPlumeSpectra.size());

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

TEST_CASE("RatioCalculationController - Evaluate without offset polynomial", "[RatioCalculationController]")
{
    const auto testWindows = GetSetupOfFitWindowsForTest();

    RatioCalculationController sut;
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references, with varied number of references in each window
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = testWindows.so2Window.ref[0].m_path;
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = false;
    ReferenceForRatioCalculation* o3Reference = GetReferenceFor(sut, StandardDoasSpecie::O3);
    o3Reference->m_path = testWindows.so2Window.ref[1].m_path;
    o3Reference->m_includeInMajor = true;
    o3Reference->m_includeInMinor = true;
    ReferenceForRatioCalculation* broReference = GetReferenceFor(sut, StandardDoasSpecie::BRO);
    broReference->m_path = testWindows.broWindow.ref[0].m_path;
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
    novac::CScanFileHandler fileHandler;
    const bool scanFileIsOk = fileHandler.CheckScanFile(TestData::GetBrORatioScanFile1());
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
    const auto result = sut.EvaluateScan(fileHandler, evaluationFileHandler.m_scan[0], fitWindowsetup);

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
    REQUIRE(result.debugInfo.plumeSpectra.size() >= sut.m_ratioEvaluationSettings.minNumberOfSpectraInPlume);
    REQUIRE(result.debugInfo.outOfPlumeSpectra.size() >= sut.m_ratioEvaluationSettings.minNumberOfReferenceSpectra);
    REQUIRE(result.debugInfo.inPlumeSpectrum.NumSpectra() == 15 * result.debugInfo.plumeSpectra.size());
    REQUIRE(result.debugInfo.outOfPlumeSpectrum.NumSpectra() == 15 * result.debugInfo.outOfPlumeSpectra.size());

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