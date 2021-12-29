#include <SpectralEvaluation/DialogControllers/NovacProgramWavelengthCalibrationController.h>
#include <SpectralEvaluation/Calibration/InstrumentCalibration.h>
#include <SpectralEvaluation/Calibration/InstrumentLineShape.h>
#include <SpectralEvaluation/Calibration/InstrumentLineShapeEstimation.h>
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

    static std::string GetMeasuredSpectrumName_I2J8549()
    {
        return GetTestDataDirectory() + std::string("I2J8549/I2J8549_170216_1230_0.pak");
    }

    static std::string GetInitialPixelToWavelengthCalibration_I2J8549()
    {
        return GetTestDataDirectory() + std::string("I2J8549/I2J8549.clb");
    }

    static std::string GetInitialInstrumentLineShapefile_I2J8549()
    {
        return GetTestDataDirectory() + std::string("I2J8549/I2J8549_302nm.slf");
    }

    static std::string GetMeasuredSpectrumName_2009175M1()
    {
        return GetTestDataDirectory() + std::string("2009175M1/2009175M1_211214_1817_0.pak");
    }

    static std::string GetInitialPixelToWavelengthCalibration_2009175M1()
    {
        return GetTestDataDirectory() + std::string("2009175M1/DD2J3040_MASTER_SO2_HP500_PPMM.txt");
    }

    static std::string GetInitialPixelToWavelengthCalibration_D2J2200()
    {
        return GetTestDataDirectory() + std::string("D2J2200/D2J2200_Master.clb");
    }

    static std::string GetSkySpectrumName_MAYP11440()
    {
        return GetTestDataDirectory() + std::string("MAYP11440/sky_0.STD");
    }

    static std::string GetDarkSpectrumName_MAYP11440()
    {
        return GetTestDataDirectory() + std::string("MAYP11440/dark_0.STD");
    }

    static std::string GetInitialPixelToWavelengthCalibration_MAYP11440()
    {
        return GetTestDataDirectory() + std::string("MAYP11440/MAYP11440_SO2_293K_Bogumil_334nm.txt");
    }

    static std::string GetMeasuredSpectrumName_FLMS14634()
    {
        return GetTestDataDirectory() + std::string("FLMS14634/00007_0.STD");
    }

    static std::string GetDarkSpectrumName_FLMS14634()
    {
        return GetTestDataDirectory() + std::string("FLMS14634/dark_0.STD");
    }

    static std::string GetInitialPixelToWavelengthCalibration_FLMS14634()
    {
        return GetTestDataDirectory() + std::string("FLMS14634/FLMS14634.clb");
    }

    static std::string GetInitialInstrumentLineShapefile_FLMS14634()
    {
        return GetTestDataDirectory() + std::string("FLMS14634/FLMS14634_302nm.slf");
    }

    static std::string GetMeasuredSpectrumName_I2J0093()
    {
        return GetTestDataDirectory() + std::string("I2P0093/00000_0.STD");
    }

    static std::string GetDarkSpectrumName_I2P0093()
    {
        return GetTestDataDirectory() + std::string("I2P0093/dark_0.STD");
    }

    static std::string GetInitialPixelToWavelengthCalibration_I2P0093()
    {
        return GetTestDataDirectory() + std::string("I2P0093/I2P0093_Master.clb");
    }

    static std::string GetInitialInstrumentLineShape_I2P0093()
    {
        return GetTestDataDirectory() + std::string("I2P0093/I2P0093_302nm_Master.slf");
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
        sut.m_inputSpectrumFile = GetMeasuredSpectrumName_I2J8549();
        sut.m_initialCalibrationFile = GetInitialPixelToWavelengthCalibration_I2J8549();
        sut.m_initialLineShapeFile = GetInitialInstrumentLineShapefile_I2J8549();
        sut.m_instrumentLineShapeFitOption = WavelengthCalibrationController::InstrumentLineShapeFitOption::None;
        sut.m_solarSpectrumFile = GetSolarAtlasFile();

        // Act
        sut.RunCalibration();
        REQUIRE(nullptr != sut.m_resultingCalibration);

        // Make sure we have a good quality result
        {
            REQUIRE(sut.m_calibrationDebug.inlierCorrespondencePixels.size() > 30);
        }

        // Make sure that RunCalibration produces a correct pixel to wavelength polynomial.
        {
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial.size() == 4);
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial[0] == Approx(278.6).margin(1.5));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial[1] == Approx(0.086).margin(1e-2));
            REQUIRE(std::abs(sut.m_resultingCalibration->pixelToWavelengthPolynomial[2]) < 1e-5);
            REQUIRE(std::abs(sut.m_resultingCalibration->pixelToWavelengthPolynomial[3]) < 1e-9);
        }

        // Make sure that RunCalibration produces a correct pixel to wavelength mapping.
        //  (notice that the uncertainty is highest for the first and the last pixels)
        {
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.size() == 2048);
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.front() == Approx(278.618).margin(0.7));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.back() == Approx(424.268).margin(0.7));
        }

        // Make sure that RunCalibration produces a highly accurate pixel to wavelength mapping in the DOAS range
        //  (notice that the uncertainty is smaller in the range where the intensity is good)
        {
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping[500] == Approx(319.764).margin(0.1));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping[700] == Approx(335.290).margin(0.1));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping[900] == Approx(350.249).margin(0.1));
        }

        // Make sure that GetFinalCalibration keeps input instrument line shape.
        {
            CCrossSectionData initialInstrumentLineShape;
            ReadCrossSectionFile(GetInitialInstrumentLineShapefile_I2J8549(), initialInstrumentLineShape);

            const auto finalCalibration = sut.GetFinalCalibration();

            REQUIRE(finalCalibration->instrumentLineShape.size() == initialInstrumentLineShape.GetSize());
            REQUIRE(finalCalibration->instrumentLineShape.front() == initialInstrumentLineShape.m_crossSection.front() / Max(initialInstrumentLineShape.m_crossSection));
            REQUIRE(finalCalibration->instrumentLineShape.back() == initialInstrumentLineShape.m_crossSection.back() / Max(initialInstrumentLineShape.m_crossSection));

            REQUIRE(finalCalibration->instrumentLineShapeGrid.size() == initialInstrumentLineShape.GetSize());
            REQUIRE(finalCalibration->instrumentLineShapeGrid.front() == initialInstrumentLineShape.m_waveLength.front());
            REQUIRE(finalCalibration->instrumentLineShapeGrid.back() == initialInstrumentLineShape.m_waveLength.back());
        }
    }

    TEST_CASE(
        "NovacProgramWavelengthCalibrationController (measured I2J8549 spectrum) with fitted instrument line shape ",
        "[NovacProgramWavelengthCalibrationController][WavelengthCalibrationController][IntegrationTest][LongRunningTest]")
    {
        NovacProgramWavelengthCalibrationController sut;
        sut.m_inputSpectrumFile = GetMeasuredSpectrumName_I2J8549();
        sut.m_initialCalibrationFile = GetInitialPixelToWavelengthCalibration_I2J8549();
        sut.m_initialLineShapeFile = "";
        sut.m_instrumentLineShapeFitOption = WavelengthCalibrationController::InstrumentLineShapeFitOption::SuperGaussian;
        sut.m_solarSpectrumFile = GetSolarAtlasFile();

        // Act
        sut.RunCalibration();
        REQUIRE(nullptr != sut.m_resultingCalibration);

        // Make sure we have a good quality result
        {
            REQUIRE(sut.m_calibrationDebug.inlierCorrespondencePixels.size() > 30);
        }

        // Make sure that RunCalibration produces a correct pixel to wavelength polynomial.
        {
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial.size() == 4);
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial[0] == Approx(278.6).margin(1.5));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial[1] == Approx(0.086).margin(1e-2));
            REQUIRE(std::abs(sut.m_resultingCalibration->pixelToWavelengthPolynomial[2]) < 1e-5);
            REQUIRE(std::abs(sut.m_resultingCalibration->pixelToWavelengthPolynomial[3]) < 1e-9);
        }

        // Make sure that RunCalibration produces a correct pixel to wavelength mapping.
        //  (notice that the uncertainty is highest for the first and the last pixels)
        {
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.size() == 2048);
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.front() == Approx(278.618).margin(0.7));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.back() == Approx(424.268).margin(0.7));
        }

        // Make sure that RunCalibration produces a highly accurate pixel to wavelength mapping in the DOAS range
        //  (notice that the uncertainty is smaller in the range where the intensity is good)
        {
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping[500] == Approx(319.764).margin(0.1));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping[700] == Approx(335.290).margin(0.1));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping[900] == Approx(350.249).margin(0.1));
        }

        // Make sure that RunCalibration produces a correct instrument line shape parametrization.
        {
            const auto* superGaussian = dynamic_cast<novac::SuperGaussianLineShape*>(sut.m_resultingCalibration->instrumentLineShapeParameter);

            REQUIRE(superGaussian != nullptr);
            REQUIRE(superGaussian->k == Approx(2.5).margin(0.6));
            REQUIRE(superGaussian->w == Approx(0.292).margin(0.05));

            CCrossSectionData actualInstrumentLineShape;
            ReadCrossSectionFile(GetInitialInstrumentLineShapefile_I2J8549(), actualInstrumentLineShape);
            REQUIRE(GetFwhm(actualInstrumentLineShape) == Approx(superGaussian->Fwhm()).margin(0.05));
        }
    }

    TEST_CASE(
        "NovacProgramWavelengthCalibrationController (measured I2P0093 spectrum) with fitted instrument line shape ",
        "[NovacProgramWavelengthCalibrationController][WavelengthCalibrationController][IntegrationTest][LongRunningTest]")
    {
        MobileDoasWavelengthCalibrationController sut;
        sut.m_inputSpectrumFile = GetMeasuredSpectrumName_I2J0093();
        sut.m_darkSpectrumFile = GetDarkSpectrumName_I2P0093();
        sut.m_initialCalibrationFile = GetInitialPixelToWavelengthCalibration_I2P0093();
        sut.m_initialLineShapeFile = GetInitialInstrumentLineShape_I2P0093();
        sut.m_instrumentLineShapeFitOption = WavelengthCalibrationController::InstrumentLineShapeFitOption::SuperGaussian;
        sut.m_solarSpectrumFile = GetSolarAtlasFile();

        // Act
        sut.RunCalibration();
        REQUIRE(nullptr != sut.m_resultingCalibration);

        // Make sure we have a good quality result
        {
            REQUIRE(sut.m_calibrationDebug.inlierCorrespondencePixels.size() > 30);
        }

        // Make sure that RunCalibration produces a correct pixel to wavelength polyno  mial.
        {
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial.size() == 4);
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial[0] == Approx(278.0).margin(1.5));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial[1] == Approx(0.086).margin(1e-2));
            REQUIRE(std::abs(sut.m_resultingCalibration->pixelToWavelengthPolynomial[2]) < 1e-5);
            REQUIRE(std::abs(sut.m_resultingCalibration->pixelToWavelengthPolynomial[3]) < 1e-9);
        }

        // Make sure that RunCalibration produces a correct pixel to wavelength mapping.
        //  (notice that the uncertainty is highest for the first and the last pixels)
        {
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.size() == 2048);
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.front() == Approx(278.0).margin(0.7));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.back() == Approx(421.6).margin(0.7));
        }

        // Make sure that RunCalibration produces a highly accurate pixel to wavelength mapping in the DOAS range
        //  (notice that the uncertainty is smaller in the range where the intensity is good)
        {
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping[500] == Approx(318.860).margin(0.1));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping[700] == Approx(334.222).margin(0.1));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping[900] == Approx(348.974).margin(0.1));
        }

        // Make sure that RunCalibration produces a correct instrument line shape parametrization.
        {
            const auto* superGaussian = dynamic_cast<novac::SuperGaussianLineShape*>(sut.m_resultingCalibration->instrumentLineShapeParameter);

            REQUIRE(superGaussian != nullptr);
            REQUIRE(superGaussian->k == Approx(2.16).margin(0.4));
            REQUIRE(superGaussian->w == Approx(0.34).margin(0.05));

            CCrossSectionData actualInstrumentLineShape;
            ReadCrossSectionFile(GetInitialInstrumentLineShape_I2P0093(), actualInstrumentLineShape);
            REQUIRE(GetFwhm(actualInstrumentLineShape) > 0.1 + superGaussian->Fwhm()); // the initial isn't that good here..
        }
    }

    TEST_CASE(
        "NovacProgramWavelengthCalibrationController (measured 2009175M1 spectrum) with no initial instrument line shape",
        "[NovacProgramWavelengthCalibrationController][WavelengthCalibrationController][IntegrationTest][LongRunningTest][AvaSpec]")
    {
        /** This is a test making sure that the wavelength calibration is able to handle a measured spectrum file with
        *   an initially unknown instrument line shape and a relatively large initial error in the pixel to wavelength mapping */

        NovacProgramWavelengthCalibrationController sut;
        sut.m_inputSpectrumFile = GetMeasuredSpectrumName_2009175M1();
        sut.m_initialCalibrationFile = GetInitialPixelToWavelengthCalibration_2009175M1();
        sut.m_instrumentLineShapeFitOption = WavelengthCalibrationController::InstrumentLineShapeFitOption::SuperGaussian;
        sut.m_instrumentLineShapeFitRegion = WavelengthRange(330.0, 350.0);
        sut.m_solarSpectrumFile = GetSolarAtlasFile();

        // Act
        sut.RunCalibration();
        REQUIRE(nullptr != sut.m_resultingCalibration);

        // Make sure we have a good quality result
        {
            REQUIRE(sut.m_calibrationDebug.inlierCorrespondencePixels.size() > 30);
        }

        // Make sure that RunCalibration produces a correct pixel to wavelength polynomial.
        {
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial.size() == 4);
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial[0] == Approx(266.0).margin(1.5));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial[1] == Approx(0.094).margin(2e-2));
            REQUIRE(std::abs(sut.m_resultingCalibration->pixelToWavelengthPolynomial[2]) < 1e-5);
            REQUIRE(std::abs(sut.m_resultingCalibration->pixelToWavelengthPolynomial[3]) < 2e-9);
        }

        // Make sure that RunCalibration produces a correct pixel to wavelength mapping
        //  (notice that the uncertainty is highest for the first and the last pixels)
        {
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.size() == 2048);
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.front() == Approx(266.0).margin(1.5));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.back() == Approx(440.3).margin(1.0));
        }

        // Make sure that RunCalibration produces a highly accurate pixel to wavelength mapping in the DOAS range
        //  (notice that the uncertainty is smaller in the range where the intensity is good)
        {
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping[500] == Approx(312.15).margin(0.25));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping[700] == Approx(330.13).margin(0.25));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping[900] == Approx(347.71).margin(0.25));
        }

        // Make sure that RunCalibration produces a correct instrument line shape parametrization.
        {
            const auto* superGaussian = dynamic_cast<novac::SuperGaussianLineShape*>(sut.m_resultingCalibration->instrumentLineShapeParameter);

            REQUIRE(superGaussian != nullptr);
            REQUIRE(superGaussian->k == Approx(4.5).margin(0.4));
            REQUIRE(superGaussian->w == Approx(0.209).margin(0.05));
        }
    }

    TEST_CASE(
    "NovacProgramWavelengthCalibrationController (measured 2009175M1 spectrum) with wrong initial instrument line shape",
    "[NovacProgramWavelengthCalibrationController][WavelengthCalibrationController][IntegrationTest][LongRunningTest][AvaSpec]")
    {
        /** This is a test making sure that the wavelength calibration is able to handle a measured spectrum file with
        *   an initial pixel-to-wavelength mapping from another device (and spectrometer model) and the instrument line shape is unknown. */

        NovacProgramWavelengthCalibrationController sut;
        sut.m_inputSpectrumFile = GetMeasuredSpectrumName_2009175M1();
        sut.m_initialCalibrationFile = GetInitialPixelToWavelengthCalibration_D2J2200();
        sut.m_instrumentLineShapeFitOption = WavelengthCalibrationController::InstrumentLineShapeFitOption::SuperGaussian;
        sut.m_instrumentLineShapeFitRegion = WavelengthRange(330.0, 350.0);
        sut.m_solarSpectrumFile = GetSolarAtlasFile();

        // Act
        sut.RunCalibration();
        REQUIRE(nullptr != sut.m_resultingCalibration);

        // Make sure we have a good quality result
        {
            REQUIRE(sut.m_calibrationDebug.inlierCorrespondencePixels.size() > 30);
        }

        // Make sure that RunCalibration produces a correct pixel to wavelength polynomial.
        {
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial.size() == 4);
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial[0] == Approx(265.5).margin(2.0));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial[1] == Approx(0.094).margin(2e-2));
            REQUIRE(std::abs(sut.m_resultingCalibration->pixelToWavelengthPolynomial[2]) < 1e-5);
            REQUIRE(std::abs(sut.m_resultingCalibration->pixelToWavelengthPolynomial[3]) < 1e-9);
        }

        // Make sure that RunCalibration produces a correct pixel to wavelength mapping
        //  (notice that the uncertainty is highest for the first and the last pixels)
        {
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.size() == 2048);
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.front() == Approx(265.5).margin(2.0));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.back() == Approx(440.3).margin(1.0));
        }

        // Make sure that RunCalibration produces a highly accurate pixel to wavelength mapping in the DOAS range
        //  (notice that the uncertainty is smaller in the range where the intensity is good)
        {
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping[500] == Approx(312.13).margin(0.25));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping[700] == Approx(330.13).margin(0.25));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping[900] == Approx(347.71).margin(0.25));
        }

        // Make sure that RunCalibration produces a correct instrument line shape parametrization.
        {
            const auto* superGaussian = dynamic_cast<novac::SuperGaussianLineShape*>(sut.m_resultingCalibration->instrumentLineShapeParameter);

            REQUIRE(superGaussian != nullptr);
            REQUIRE(superGaussian->k == Approx(4.5).margin(0.4));
            REQUIRE(superGaussian->w == Approx(0.209).margin(0.05));
        }
    }

    TEST_CASE(
        "MobileDoasWavelengthCalibrationController (measured MAYP11440 spectrum) with no initial instrument line shape",
        "[MobileDoasWavelengthCalibrationController][WavelengthCalibrationController][IntegrationTest][LongRunningTest][MayaPro]")
    {
        /** This is a test making sure that the wavelength calibration is able to handle a measured spectrum file with
        *   an initially unknown instrument line shape and a relatively large initial error in the pixel to wavelength mapping */

        MobileDoasWavelengthCalibrationController sut;
        sut.m_inputSpectrumFile = GetSkySpectrumName_MAYP11440();
        sut.m_darkSpectrumFile = GetDarkSpectrumName_MAYP11440();
        sut.m_initialCalibrationFile = GetInitialPixelToWavelengthCalibration_MAYP11440();
        sut.m_instrumentLineShapeFitOption = WavelengthCalibrationController::InstrumentLineShapeFitOption::SuperGaussian;
        sut.m_instrumentLineShapeFitRegion = WavelengthRange(330.0, 350.0);
        sut.m_solarSpectrumFile = GetSolarAtlasFile();

        // Act
        sut.RunCalibration();
        REQUIRE(nullptr != sut.m_resultingCalibration);

        // Make sure we have a good quality result
        {
            REQUIRE(sut.m_calibrationDebug.inlierCorrespondencePixels.size() > 30);
        }

        // Make sure that RunCalibration produces a correct pixel to wavelength polynomial.
        {
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial.size() == 4);
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial[0] == Approx(280.2).margin(0.3));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial[1] == Approx(0.052).margin(2e-2));
            REQUIRE(std::abs(sut.m_resultingCalibration->pixelToWavelengthPolynomial[2]) < 1e-5);
            REQUIRE(std::abs(sut.m_resultingCalibration->pixelToWavelengthPolynomial[3]) < 1e-9);
        }

        // Make sure that RunCalibration produces a correct pixel to wavelength mapping
        //  (notice that the uncertainty is highest for the first and the last pixels)
        {
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.size() == 2068);
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.front() == Approx(280.2).margin(0.7));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.back() == Approx(377.5).margin(0.7));
        }

        // Make sure that RunCalibration produces a highly accurate pixel to wavelength mapping in the DOAS range
        //  (notice that the uncertainty is smaller in the range where the intensity is good)
        {
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping[500] == Approx(305.82).margin(0.1));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping[700] == Approx(315.68).margin(0.1));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping[900] == Approx(325.34).margin(0.1));
        }

        // Make sure that RunCalibration produces a correct instrument line shape parametrization.
        {
            const auto* superGaussian = dynamic_cast<novac::SuperGaussianLineShape*>(sut.m_resultingCalibration->instrumentLineShapeParameter);

            REQUIRE(superGaussian != nullptr);
            REQUIRE(superGaussian->k == Approx(3.21).margin(0.4));
            REQUIRE(superGaussian->w == Approx(0.231).margin(0.05));
        }
    }

    TEST_CASE(
        "NovacProgramWavelengthCalibrationController (measured FLMS14634 spectrum) with fitted instrument line shape ",
        "[NovacProgramWavelengthCalibrationController][WavelengthCalibrationController][IntegrationTest][LongRunningTest][Flame]")
    {
        MobileDoasWavelengthCalibrationController sut;
        sut.m_inputSpectrumFile = GetMeasuredSpectrumName_FLMS14634();
        sut.m_darkSpectrumFile = GetDarkSpectrumName_FLMS14634();
        sut.m_initialCalibrationFile = GetInitialPixelToWavelengthCalibration_FLMS14634(); // notice that the initial error here is very large in parts of the spectrum.
        sut.m_initialLineShapeFile = GetInitialInstrumentLineShapefile_FLMS14634();
        sut.m_instrumentLineShapeFitOption = WavelengthCalibrationController::InstrumentLineShapeFitOption::SuperGaussian;
        sut.m_solarSpectrumFile = GetSolarAtlasFile();

        // Act
        sut.RunCalibration();
        REQUIRE(nullptr != sut.m_resultingCalibration);

        // Make sure we have a good quality result
        {
            // TODO: Find a way to improve this
            REQUIRE(sut.m_calibrationDebug.inlierCorrespondencePixels.size() > 15);
        }

        // Make sure that RunCalibration produces a correct pixel to wavelength polynomial.
        {
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial.size() == 4);
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial[0] == Approx(278.2).margin(1.5));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthPolynomial[1] == Approx(0.082).margin(1e-2));
            REQUIRE(std::abs(sut.m_resultingCalibration->pixelToWavelengthPolynomial[2]) < 1e-5);
            REQUIRE(std::abs(sut.m_resultingCalibration->pixelToWavelengthPolynomial[3]) < 2e-9);
        }

        // Make sure that RunCalibration produces a correct pixel to wavelength mapping.
        //  (notice that the uncertainty is highest for the first and the last pixels)
        {
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.size() == 2048);
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.front() == Approx(278.2).margin(1.5));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.back() == Approx(420.7).margin(1.5));
        }

        // Make sure that RunCalibration produces a highly accurate pixel to wavelength mapping in the DOAS range
        //  (notice that the uncertainty is smaller in the range where the intensity is good)
        {
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping[500] == Approx(317.75).margin(0.20));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping[700] == Approx(332.80).margin(0.20));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping[900] == Approx(347.50).margin(0.20));
        }

        // Make sure that RunCalibration produces a correct instrument line shape parametrization.
        {
            const auto* superGaussian = dynamic_cast<novac::SuperGaussianLineShape*>(sut.m_resultingCalibration->instrumentLineShapeParameter);

            REQUIRE(superGaussian != nullptr);
            REQUIRE(superGaussian->k == Approx(1.86).margin(0.4));
            REQUIRE(superGaussian->w == Approx(0.35).margin(0.05));
        }
    }
}