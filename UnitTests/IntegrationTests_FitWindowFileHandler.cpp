#include <SpectralEvaluation/Evaluation/PlumeSpectrumSelector.h>
#include <SpectralEvaluation/File/FitWindowFileHandler.h>
#include "catch.hpp"
#include "TestData.h"

namespace novac
{
    TEST_CASE("FitWindowFileHandler can read .nfw file for SO2 from the NovacProgram", "[FitWindowFileHandler][IntegrationTest]")
    {
        const std::string& filename = TestData::GetBrORatioFitWindowFileSO2();

        CFitWindowFileHandler sut;

        const auto result = sut.ReadFitWindowFile(filename);

        REQUIRE(result.size() == 1);

        // Make sure that the contents is as expected
        REQUIRE(result.front().name == "SO2");
        REQUIRE(result.front().fitLow == 442);
        REQUIRE(result.front().fitHigh == 595);
        REQUIRE(result.front().polyOrder == 3);
        REQUIRE(result.front().includeIntensitySpacePolyominal == true);
        REQUIRE(result.front().ringCalculation == RING_CALCULATION_OPTION::CALCULATE_RING_X2);

        REQUIRE(result.front().fitType == FIT_TYPE::FIT_POLY);
        REQUIRE(result.front().channel == 0);
        REQUIRE(result.front().specLength == 2048);

        REQUIRE(result.front().findOptimalShift == 0);
        REQUIRE(result.front().UV == 1);
        REQUIRE(result.front().shiftSky == 0);
        REQUIRE(result.front().interlaceStep == 1);

        REQUIRE(result.front().fraunhoferRef.m_path == "../TestData/BrORatio/D2J2124_SolarSpec_Master.txt");

        // finally the references
        REQUIRE(result.front().nRef == 2);
        REQUIRE(result.front().ref[0].m_specieName == "SO2");
        REQUIRE(result.front().ref[0].m_path == "../TestData/BrORatio/D2J2124_SO2_Bogumil_293K_Master.txt");
        REQUIRE(result.front().ref[1].m_specieName == "O3");
        REQUIRE(result.front().ref[1].m_path == "../TestData/BrORatio/D2J2124_O3_Voigt_223K_Master.txt");
    }

    TEST_CASE("FitWindowFileHandler can read .nfw file for BrO from the NovacProgram", "[FitWindowFileHandler][IntegrationTest]")
    {
        const std::string& filename = TestData::GetBrORatioFitWindowFileBrO();

        CFitWindowFileHandler sut;

        const auto result = sut.ReadFitWindowFile(filename);

        REQUIRE(result.size() == 1);

        // Make sure that the contents is as expected
        REQUIRE(result.front().name == "BrO");
        REQUIRE(result.front().fitLow == 644);
        REQUIRE(result.front().fitHigh == 923);
        REQUIRE(result.front().polyOrder == 2);
        REQUIRE(result.front().includeIntensitySpacePolyominal == true);
        REQUIRE(result.front().ringCalculation == RING_CALCULATION_OPTION::CALCULATE_RING_X2);

        REQUIRE(result.front().fitType == FIT_TYPE::FIT_POLY);
        REQUIRE(result.front().channel == 0);
        REQUIRE(result.front().specLength == 2048);

        REQUIRE(result.front().findOptimalShift == 0);
        REQUIRE(result.front().UV == 1);
        REQUIRE(result.front().shiftSky == 0);
        REQUIRE(result.front().interlaceStep == 1);

        REQUIRE(result.front().fraunhoferRef.m_path == "../TestData/BrORatio/D2J2124_SolarSpec_Master.txt");

        // finally the references
        REQUIRE(result.front().nRef == 4);
        REQUIRE(result.front().ref[0].m_specieName == "BrO");
        REQUIRE(result.front().ref[0].m_path == "../TestData/BrORatio/D2J2124_BrO_Fleischmann_298K.txt");
        REQUIRE(result.front().ref[1].m_specieName == "SO2");
        REQUIRE(result.front().ref[1].m_path == "../TestData/BrORatio/D2J2124_SO2_Bogumil_293K_Master.txt");
        REQUIRE(result.front().ref[2].m_specieName == "CHO2");
        REQUIRE(result.front().ref[2].m_path == "../TestData/BrORatio/D2J2124_CH2O_MellerMoortgat_298K.txt");
        REQUIRE(result.front().ref[3].m_specieName == "O3");
        REQUIRE(result.front().ref[3].m_path == "../TestData/BrORatio/D2J2124_O3_Voigt_223K_Master.txt");
    }

    TEST_CASE("FitWindowFileHandler can read back one saved fit window", "[FitWindowFileHandler][IntegrationTest]")
    {
        const std::string filename = TestData::GetTemporaryFitWindowFileName();

        // Setup a fit-window to write to file. Add some odd values here such that we know we're not getting defaults back
        CFitWindow originalWindow;
        originalWindow.channel = 1;
        originalWindow.findOptimalShift = 1;
        originalWindow.fitHigh = 987;
        originalWindow.fitLow = 765;
        originalWindow.fitType = FIT_TYPE::FIT_HP_SUB;
        originalWindow.fraunhoferRef.m_path = "C://Temp//Some//File.txt";
        originalWindow.includeIntensitySpacePolyominal = true;
        originalWindow.interlaceStep = 3;
        originalWindow.name = "Some type of window";
        originalWindow.polyOrder = 9;
        originalWindow.ringCalculation = RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING;
        originalWindow.shiftSky = 2;
        originalWindow.skyShift = 3.141;
        originalWindow.skySqueeze = 0.98765;
        originalWindow.specLength = 1234;
        originalWindow.UV = 0;
        originalWindow.nRef = 6;

        // Act: save the fit window to file and then read it back again.
        CFitWindowFileHandler sut;
        sut.WriteFitWindow(originalWindow, filename, true);
        const auto readBackWindow = sut.ReadFitWindowFile(filename);

        // Assert that the read back fit window has the same contents as the original.
        REQUIRE(readBackWindow.size() == 1);
        REQUIRE(originalWindow.channel == readBackWindow.front().channel);
        REQUIRE(originalWindow.findOptimalShift == readBackWindow.front().findOptimalShift);
        REQUIRE(originalWindow.fitHigh == readBackWindow.front().fitHigh);
        REQUIRE(originalWindow.fitLow == readBackWindow.front().fitLow);
        REQUIRE(originalWindow.fitType == readBackWindow.front().fitType);
        REQUIRE(originalWindow.fraunhoferRef.m_path == readBackWindow.front().fraunhoferRef.m_path);
        REQUIRE(originalWindow.includeIntensitySpacePolyominal == readBackWindow.front().includeIntensitySpacePolyominal);
        REQUIRE(originalWindow.interlaceStep == readBackWindow.front().interlaceStep);
        REQUIRE(originalWindow.name == readBackWindow.front().name);
        REQUIRE(originalWindow.polyOrder == readBackWindow.front().polyOrder);
        REQUIRE(originalWindow.ringCalculation == readBackWindow.front().ringCalculation);
        REQUIRE(originalWindow.shiftSky == readBackWindow.front().shiftSky);
        REQUIRE(originalWindow.skyShift == readBackWindow.front().skyShift);
        REQUIRE(originalWindow.skySqueeze == readBackWindow.front().skySqueeze);
        REQUIRE(originalWindow.specLength == readBackWindow.front().specLength);
        REQUIRE(originalWindow.UV == readBackWindow.front().UV);

    }
}