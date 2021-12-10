#include <SpectralEvaluation/DialogControllers/NovacProgramWavelengthCalibrationController.h>
#include <SpectralEvaluation/Calibration/InstrumentCalibration.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/File/File.h>
#include <SpectralEvaluation/VectorUtils.h>
#include "catch.hpp"

namespace novac
{
    static std::string GetTestDataDirectory()
    {
#ifdef _MSC_VER
        return std::string("../TestData/");
#else
        return std::string("TestData/");
#endif // _MSC_VER 
    }

    static std::string GetMeasuredSpectrumName()
    {
        return GetTestDataDirectory() + std::string("I2J8549/I2J8549_170216_1230_0.pak");
    }

    static std::string GetInitialPixelToWavelengthCalibration()
    {
        return GetTestDataDirectory() + std::string("I2J8549/I2J8549.clb");
    }

    static std::string GetInitialInstrumentLineShapefile()
    {
        return GetTestDataDirectory() + std::string("I2J8549/I2J8549_302nm.slf");
    }

    static std::string GetSolarAtlasFile()
    {
        return GetTestDataDirectory() + std::string("SOLARFL_296-440nm.xs");
    }

    TEST_CASE(
        "NovacProgramWavelengthCalibrationController (measured I2J8549 spectrum) with no ILF fit",
        "[NovacProgramWavelengthCalibrationController][WavelengthCalibrationController][IntegrationTest][LongRunningTest]")
    {
        NovacProgramWavelengthCalibrationController sut;
        sut.m_inputSpectrumFile = GetMeasuredSpectrumName();
        sut.m_initialCalibrationFile = GetInitialPixelToWavelengthCalibration();
        sut.m_initialLineShapeFile = GetInitialInstrumentLineShapefile();
        sut.m_instrumentLineShapeFitOption = WavelengthCalibrationController::InstrumentLineShapeFitOption::None;
        sut.m_solarSpectrumFile = GetSolarAtlasFile();

        // Act
        sut.RunCalibration();
        REQUIRE(nullptr != sut.m_resultingCalibration);

        // Make sure that RunCalibration produces a correct pixel to wavelength polynomial.
        {
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial.size() == 4);
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial[0] == Approx(278.618).margin(0.01));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial[1] == Approx(0.0855182).margin(1e-6));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial[2] == Approx(-6.24709e-06).margin(1e-9));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial[3] == Approx(-3.76497e-10).margin(1e-14));
        }

        // Make sure that RunCalibration produces a correct pixel to wavelength mapping.
        {
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.size() == 2048);
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.front() == Approx(278.618).margin(0.01));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.back() == Approx(424.268).margin(0.01));
        }

        // Make sure that GetFinalCalibration keeps input instrument line shape.
        {
            CCrossSectionData initialInstrumentLineShape;
            ReadCrossSectionFile(GetInitialInstrumentLineShapefile(), initialInstrumentLineShape);
            
            const auto finalCalibration = sut.GetFinalCalibration();

            REQUIRE(finalCalibration->instrumentLineShape.size() == initialInstrumentLineShape.GetSize());
            REQUIRE(finalCalibration->instrumentLineShape.front() == initialInstrumentLineShape.m_crossSection.front() / Max(initialInstrumentLineShape.m_crossSection));
            REQUIRE(finalCalibration->instrumentLineShape.back() == initialInstrumentLineShape.m_crossSection.back() / Max(initialInstrumentLineShape.m_crossSection));

            REQUIRE(finalCalibration->instrumentLineShapeGrid.size() == initialInstrumentLineShape.GetSize());
            REQUIRE(finalCalibration->instrumentLineShapeGrid.front() == initialInstrumentLineShape.m_waveLength.front());
            REQUIRE(finalCalibration->instrumentLineShapeGrid.back() == initialInstrumentLineShape.m_waveLength.back());
        }
    }
}