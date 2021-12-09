#include <SpectralEvaluation/DialogControllers/InstrumentLineshapeCalibrationController.h>
#include <SpectralEvaluation/Calibration/InstrumentCalibration.h>
#include <SpectralEvaluation/Calibration/InstrumentLineShape.h>
#include <SpectralEvaluation/File/File.h>
#include <SpectralEvaluation/Calibration/InstrumentLineShapeEstimation.h>
#include "catch.hpp"

namespace novac
{
    static std::string GetTestFileName()
    {
#ifdef _MSC_VER
        return std::string("../TestData/hglampnov152021.std");
#else
        return std::string("./TestData/hglampnov152021.std");
#endif // _MSC_VER
    }

    static std::string GetDarkFileName()
    {
#ifdef _MSC_VER
        return std::string("../TestData/hglampnov152021_dark.std");
#else
        return std::string("./TestData/hglampnov152021_dark.std");
#endif // _MSC_VER
    }

    static std::string GetInstrumentCalibrationStdFileName()
    {
#ifdef _MSC_VER
        return std::string("../TestData/InstrumentCalibration.std");
#else
        return std::string("TestData/InstrumentCalibration.std");
#endif // _MSC_VER
    }


    TEST_CASE(
        "InstrumentLineshapeCalibrationController with simple mercury spectrum file with no wavelength calibration in file",
        "[InstrumentLineshapeCalibrationController][IntegrationTest]")
    {
        InstrumentLineshapeCalibrationController sut;
        sut.m_inputSpectrumPath = GetTestFileName();
        sut.m_darkSpectrumPath = GetDarkFileName();

        SECTION("Update - reads dark corrected mercury spectrum and locates emission lines.")
        {
            sut.m_readWavelengthCalibrationFromFile = true; // read the calibration as it is in the file.

            sut.Update();

            // The spectrum should be read in and dark-corrected  (i.e. difference between first value in input-spectrum and dark-spectrum).
            REQUIRE(sut.m_inputSpectrum.size() == 2048);
            REQUIRE(sut.m_inputSpectrum.front() == Approx(39.180567858));

            // There is no inherent calibration in the file, hence the sut should produce one.
            REQUIRE(sut.m_inputSpectrumContainsWavelength == false);
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.size() == 2048);
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.front() == 0);
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.back() == 2047);

            // The sut should have located the emission lines in the measurement at the correct pixels.
            //  Since there is no wavelength calibration, the wavelengths should equal the pixel values.
            REQUIRE(sut.m_peaksFound.size() == 5);
            REQUIRE(sut.m_peaksFound[0].pixel == Approx(81.38));
            REQUIRE(sut.m_peaksFound[0].wavelength == sut.m_peaksFound[0].pixel);
            REQUIRE(sut.m_peaksFound[1].pixel == Approx(168.7630));
            REQUIRE(sut.m_peaksFound[1].wavelength == sut.m_peaksFound[1].pixel);
            REQUIRE(sut.m_peaksFound[2].pixel == Approx(234.63414));
            REQUIRE(sut.m_peaksFound[2].wavelength == sut.m_peaksFound[2].pixel);
            REQUIRE(sut.m_peaksFound[3].pixel == Approx(634.02));
            REQUIRE(sut.m_peaksFound[3].wavelength == sut.m_peaksFound[3].pixel);
            REQUIRE(sut.m_peaksFound[4].pixel == Approx(1690.24349));
            REQUIRE(sut.m_peaksFound[4].wavelength == sut.m_peaksFound[4].pixel);
        }

        SECTION("Update - reads dark corrected mercury spectrum and locates emission lines.")
        {
            sut.m_readWavelengthCalibrationFromFile = false; // i.e. perform a calibration 

            sut.Update();

            // The spectrum should be read in and dark-corrected  (i.e. difference between first value in input-spectrum and dark-spectrum).
            REQUIRE(sut.m_inputSpectrum.size() == 2048);
            REQUIRE(sut.m_inputSpectrum.front() == Approx(39.180567858));

            // There is no inherent calibration in the file, hence the sut should produce one.
            REQUIRE(sut.m_inputSpectrumContainsWavelength == false);
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.size() == 2048);
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.front() == Approx(282.5505124798));
            REQUIRE(sut.m_resultingCalibration->pixelToWavelengthMapping.back() == Approx(421.7909430782));

            // The sut should have located the emission lines in the measurement at the correct pixels
            //  and allocated a correct wavelength value to them (based on known mercury emission lines).
            REQUIRE(sut.m_peaksFound.size() == 5);
            REQUIRE(sut.m_peaksFound[0].pixel == Approx(81.38));
            REQUIRE(sut.m_peaksFound[0].wavelength == Approx(289.4));
            REQUIRE(sut.m_peaksFound[1].pixel == Approx(168.7630));
            REQUIRE(sut.m_peaksFound[1].wavelength == Approx(296.69571));
            REQUIRE(sut.m_peaksFound[2].pixel == Approx(234.63414));
            REQUIRE(sut.m_peaksFound[2].wavelength == Approx(302.1498));
            REQUIRE(sut.m_peaksFound[3].pixel == Approx(634.02));
            REQUIRE(sut.m_peaksFound[3].wavelength == Approx(334.148));
            REQUIRE(sut.m_peaksFound[4].pixel == Approx(1690.24349));
            REQUIRE(sut.m_peaksFound[4].wavelength == Approx(404.6563));
        }

        SECTION("SaveResultAsStd - No function fitted - saves the measured instrument line shape to a file")
        {
            sut.m_readWavelengthCalibrationFromFile = false; // i.e. perform a calibration 
            sut.Update();   // read the files and locate the emission lines
            sut.FitFunctionToLineShape(2, InstrumentLineshapeCalibrationController::LineShapeFunction::None);

            // Act
            sut.SaveResultAsStd(2, GetInstrumentCalibrationStdFileName()); // save peak number 2, the 302nm, to file.

            // Now read in the calibration again from file and make sure that we get the expected result back. 
            InstrumentCalibration readCalibration;
            bool readSuccessfully = ReadInstrumentCalibration(GetInstrumentCalibrationStdFileName(), readCalibration);

            REQUIRE(readSuccessfully);

            REQUIRE(readCalibration.pixelToWavelengthMapping.size() == 2048);
            REQUIRE(readCalibration.pixelToWavelengthMapping.front() == Approx(282.5505));
            REQUIRE(readCalibration.pixelToWavelengthMapping.back() == Approx(421.7909));

            REQUIRE(readCalibration.instrumentLineShapeParameter == nullptr);
            REQUIRE(readCalibration.instrumentLineShape.size() == 33);
            REQUIRE(readCalibration.instrumentLineShapeGrid.size() == 33);
            REQUIRE(readCalibration.instrumentLineShapeCenter == Approx(302.190));
            REQUIRE(GetFwhm(readCalibration.instrumentLineShapeGrid, readCalibration.instrumentLineShape) == Approx(0.70).margin(0.01));
        }

        SECTION("SaveResultAsStd - Super Gaussian function fitted - saves the gaussian instrument line shape to a file")
        {
            sut.m_readWavelengthCalibrationFromFile = false; // i.e. perform a calibration 
            sut.Update();   // read the files and locate the emission lines
            sut.FitFunctionToLineShape(2, InstrumentLineshapeCalibrationController::LineShapeFunction::SuperGauss);

            // Act
            sut.SaveResultAsStd(2, GetInstrumentCalibrationStdFileName()); // save peak number 2, the 302nm, to file.

            // Now read in the calibration again from file and make sure that we get the expected result back. 
            InstrumentCalibration readCalibration;
            bool readSuccessfully = ReadInstrumentCalibration(GetInstrumentCalibrationStdFileName(), readCalibration);

            REQUIRE(readSuccessfully);

            REQUIRE(readCalibration.pixelToWavelengthMapping.size() == 2048);
            REQUIRE(readCalibration.pixelToWavelengthMapping.front() == Approx(282.5505));
            REQUIRE(readCalibration.pixelToWavelengthMapping.back() == Approx(421.7909));

            REQUIRE(readCalibration.instrumentLineShapeParameter != nullptr);
            REQUIRE(readCalibration.instrumentLineShapeParameter->Type() == novac::InstrumentLineShapeFunctionType::SuperGaussian);
            const auto* superGaussian = static_cast<SuperGaussianLineShape*>(readCalibration.instrumentLineShapeParameter);
            REQUIRE(superGaussian->Fwhm() == Approx(0.70).margin(0.01));

            REQUIRE(readCalibration.instrumentLineShape.size() == 33);
            REQUIRE(readCalibration.instrumentLineShapeGrid.size() == 33);
            REQUIRE(readCalibration.instrumentLineShapeCenter == Approx(302.150));
            REQUIRE(GetFwhm(readCalibration.instrumentLineShapeGrid, readCalibration.instrumentLineShape) == Approx(0.70).margin(0.01));
        }
    }
}
