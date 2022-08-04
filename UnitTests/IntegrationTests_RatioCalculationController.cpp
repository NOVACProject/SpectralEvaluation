#include <SpectralEvaluation/DialogControllers/RatioCalculationController.h>
#include <SpectralEvaluation/File/FitWindowFileHandler.h>

#include "catch.hpp"
#include "TestData.h"

using namespace novac;

#pragma region helper methods

ReferenceForRatioCalculation* GetReferenceFor(RatioCalculationController& controller, StandardDoasSpecie specie) {
    for (auto& ref : controller.m_references)
    {
        if (ref.specie == specie)
        {
            return &ref;
        }
    }

    return nullptr; // not found
}
#pragma endregion

TEST_CASE("RatioCalculationController - SetupFitWindows - SO2 included in both windows", "[RatioCalculationController]")
{
    // Read in one fit window. This helps us with ready-to-use references..
    CFitWindowFileHandler fitWindowFileHandler;
    auto allWindows = fitWindowFileHandler.ReadFitWindowFile(TestData::GetBrORatioFitWindowFileSO2());
    REQUIRE(allWindows.size() == 1);
    auto so2FitWindow = allWindows.front();

    RatioCalculationController sut;
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references
    ReferenceForRatioCalculation *so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = so2FitWindow.ref[0].m_path;
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = true;
    GetReferenceFor(sut, StandardDoasSpecie::RING)->m_automaticallyCalculate = false;
    GetReferenceFor(sut, StandardDoasSpecie::RING_LAMBDA4)->m_automaticallyCalculate = false;


    // Act
    const auto result = sut.SetupFitWindows();

    // Assert
    // SO2 window
    REQUIRE(result->m_so2Window.name == "SO2");
    REQUIRE(result->m_so2Window.polyOrder == 3);
    REQUIRE(result->m_so2Window.fitType == novac::FIT_TYPE::FIT_POLY);
    REQUIRE(result->m_so2Window.nRef == 1);
    REQUIRE(result->m_so2Window.ref[0].m_path == so2Reference->m_path);
    REQUIRE(result->m_so2Window.ref[0].m_specieName == "SO2");
    REQUIRE(result->m_so2Window.ref[0].m_shiftOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->m_so2Window.ref[0].m_squeezeOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->m_so2Window.ringCalculation == novac::RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING);


    // BrO window
    REQUIRE(result->m_broWindow.name == "BrO");
    REQUIRE(result->m_broWindow.polyOrder == 3);
    REQUIRE(result->m_broWindow.fitType == novac::FIT_TYPE::FIT_POLY);
    REQUIRE(result->m_broWindow.nRef == 1);
    REQUIRE(result->m_broWindow.ref[0].m_path == so2Reference->m_path);
    REQUIRE(result->m_broWindow.ref[0].m_specieName == "SO2");
    REQUIRE(result->m_broWindow.ref[0].m_shiftOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->m_broWindow.ref[0].m_squeezeOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->m_broWindow.ringCalculation == novac::RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING);
}

TEST_CASE("RatioCalculationController - SetupFitWindows - SO2 included in only major window", "[RatioCalculationController]")
{
    // Read in one fit window. This helps us with ready-to-use references..
    CFitWindowFileHandler fitWindowFileHandler;
    auto allWindows = fitWindowFileHandler.ReadFitWindowFile(TestData::GetBrORatioFitWindowFileSO2());
    REQUIRE(allWindows.size() == 1);
    auto so2FitWindow = allWindows.front();

    RatioCalculationController sut;
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references
    ReferenceForRatioCalculation* so2Reference = GetReferenceFor(sut, StandardDoasSpecie::SO2);
    so2Reference->m_path = so2FitWindow.ref[0].m_path;
    so2Reference->m_includeInMajor = true;
    so2Reference->m_includeInMinor = false;
    GetReferenceFor(sut, StandardDoasSpecie::RING)->m_automaticallyCalculate = false;
    GetReferenceFor(sut, StandardDoasSpecie::RING_LAMBDA4)->m_automaticallyCalculate = false;

    // Act
    const auto result = sut.SetupFitWindows();

    // Assert
    // SO2 window
    REQUIRE(result->m_so2Window.name == "SO2");
    REQUIRE(result->m_so2Window.polyOrder == 3);
    REQUIRE(result->m_so2Window.fitType == novac::FIT_TYPE::FIT_POLY);
    REQUIRE(result->m_so2Window.nRef == 1);
    REQUIRE(result->m_so2Window.ref[0].m_specieName == "SO2");
    REQUIRE(result->m_so2Window.ref[0].m_shiftOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->m_so2Window.ref[0].m_squeezeOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->m_so2Window.ringCalculation == novac::RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING);


    // BrO window
    REQUIRE(result->m_broWindow.name == "BrO");
    REQUIRE(result->m_broWindow.polyOrder == 3);
    REQUIRE(result->m_broWindow.fitType == novac::FIT_TYPE::FIT_POLY);
    REQUIRE(result->m_broWindow.nRef == 0);
    REQUIRE(result->m_broWindow.ref[0].m_shiftOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->m_broWindow.ref[0].m_squeezeOption == novac::SHIFT_TYPE::SHIFT_FIX);
    REQUIRE(result->m_broWindow.ringCalculation == novac::RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING);
}

TEST_CASE("RatioCalculationController - SetupFitWindows - Calculated Ring included in both windows", "[RatioCalculationController]")
{
    RatioCalculationController sut;
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references
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
    REQUIRE(result->m_so2Window.name == "SO2");
    REQUIRE(result->m_so2Window.polyOrder == 3);
    REQUIRE(result->m_so2Window.fitType == novac::FIT_TYPE::FIT_POLY);
    REQUIRE(result->m_so2Window.nRef == 0);
    REQUIRE(result->m_so2Window.ringCalculation == novac::RING_CALCULATION_OPTION::CALCULATE_RING);


    // BrO window
    REQUIRE(result->m_broWindow.name == "BrO");
    REQUIRE(result->m_broWindow.polyOrder == 3);
    REQUIRE(result->m_broWindow.fitType == novac::FIT_TYPE::FIT_POLY);
    REQUIRE(result->m_broWindow.nRef == 0);
    REQUIRE(result->m_broWindow.ringCalculation == novac::RING_CALCULATION_OPTION::CALCULATE_RING);
}

TEST_CASE("RatioCalculationController - SetupFitWindows - Calculated Ring and Ring*Lambda4 included in both windows", "[RatioCalculationController]")
{
    RatioCalculationController sut;
    REQUIRE(sut.m_references.size() == 5); // check assumption here.

    // Setup the references
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
    REQUIRE(result->m_so2Window.name == "SO2");
    REQUIRE(result->m_so2Window.polyOrder == 3);
    REQUIRE(result->m_so2Window.fitType == novac::FIT_TYPE::FIT_POLY);
    REQUIRE(result->m_so2Window.nRef == 0);
    REQUIRE(result->m_so2Window.ringCalculation == novac::RING_CALCULATION_OPTION::CALCULATE_RING_X2);


    // BrO window
    REQUIRE(result->m_broWindow.name == "BrO");
    REQUIRE(result->m_broWindow.polyOrder == 3);
    REQUIRE(result->m_broWindow.fitType == novac::FIT_TYPE::FIT_POLY);
    REQUIRE(result->m_broWindow.nRef == 0);
    REQUIRE(result->m_broWindow.ringCalculation == novac::RING_CALCULATION_OPTION::CALCULATE_RING_X2);
}


/*
TEST_CASE("RatioCalculationController - Evaluate", "[RatioCalculationController]")
{

    // Prepare the test by reading in the .pak-file and the evaluation result and calculate the plume-properties from the result.
    novac::CScanFileHandler fileHandler;
    const bool scanFileIsOk = fileHandler.CheckScanFile(TestData::GetBrORatioScanFile1());
    REQUIRE(scanFileIsOk); // check assumption on the setup

    novac::CScanEvaluationLogFileHandler evaluationFileHandler;
    const bool evaluationFileIsOk = evaluationFileHandler.ReadEvaluationLog(TestData::GetBrORatioEvaluationFile1());
    REQUIRE(evaluationFileIsOk); // check assumption on the setup
    REQUIRE(evaluationFileHandler.m_scan.size() == 1); // check assumption on the setup


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
    RatioCalculationController sut;
    sut.m_so2FitRange = WavelengthRange{ 314.8, 326.8 }; // TODO: decide on values here!
    sut.m_broFitRange = WavelengthRange{ 330.6, 352.8 };

    ReferenceForRatioCalculation so2Reference{StandardDoasSpecie::SO2, "SO2", so2FitWindow.ref[0].m_path, true, true, false};

    sut.m_references.push_back(so2Reference);

    std::string errorMessage;

    // Act
    const bool result = sut.EvaluateScan(fileHandler, evaluationFileHandler.m_scan[0]);

    // Assert
    REQUIRE(result.size() == 1);
    REQUIRE(errorMessage.empty());

    // Verify that the result is indeed correct!
    REQUIRE(result.front().minorSpecieName == "BrO");
    REQUIRE(result.front().majorSpecieName == "SO2");
    REQUIRE(std::abs(result.front().ratio) < 1.2e-4); // During development there may be variations in the result. Just verify the range is ok for now.
    REQUIRE(std::abs(result.front().ratio) > 6e-5); // During development there may be variations in the result. Just verify the range is ok for now.
    REQUIRE(std::abs(result.front().error) < std::abs(result.front().ratio / 2.0));

}
*/