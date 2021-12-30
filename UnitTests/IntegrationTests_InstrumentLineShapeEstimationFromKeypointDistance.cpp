#include "catch.hpp"
#include <string.h>
#include <SpectralEvaluation/Calibration/InstrumentLineShapeEstimation.h>
#include <SpectralEvaluation/Calibration/FraunhoferSpectrumGeneration.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/File/STDFile.h>

namespace novac
{
    std::vector<double> GetPixelToWavelengthMappingFromFile(const std::string& clbFile);

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

    static std::string GetMeasuredSpectrumName_2009175M1()
    {
        return GetTestDataDirectory() + std::string("2009175M1/2009175M1_211214_1817_0.pak");
    }

    static std::string GetInitialPixelToWavelengthCalibration_2009175M1()
    {
        return GetTestDataDirectory() + std::string("2009175M1/DD2J3040_MASTER_SO2_HP500_PPMM.txt");
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

    static std::string GetSolarAtlasFile()
    {
        return GetTestDataDirectory() + std::string("SOLARFL_296-440nm.xs");
    }

    static void GetMeasuredPakFileSpectrum(const std::string& fileName, CSpectrum& measuredSpectrum)
    {
        novac::CScanFileHandler pakFileHandler;
        if (!pakFileHandler.CheckScanFile(fileName))
        {
            throw std::invalid_argument("Test failed, cannot read provided measured spectrum file.");
        }

        if (pakFileHandler.GetSky(measuredSpectrum))
        {
            throw std::invalid_argument("Cannot read the provided input spectrum file");
        }

        CSpectrum darkSpectrum;
        if (pakFileHandler.GetDark(darkSpectrum))
        {
            throw std::invalid_argument("Cannot read the provided input spectrum file");
        }

        measuredSpectrum.Sub(darkSpectrum);
    }

    static void GetMeasuredSpectrum(const std::string& fileName, const std::string& darkSpectrumFileName, CSpectrum& measuredSpectrum)
    {
        if (!CSTDFile::ReadSpectrum(measuredSpectrum, fileName))
        {
            throw std::invalid_argument("Failed to read spectrum for test.");
        }

        CSpectrum darkSpectrum;
        if (!CSTDFile::ReadSpectrum(darkSpectrum, darkSpectrumFileName))
        {
            throw std::invalid_argument("Failed to read dark spectrum for test.");
        }

        measuredSpectrum.Sub(darkSpectrum);
    }

    TEST_CASE(
        "InstrumentLineShapeEstimationFromKeypointDistance (measured I2J8549 spectrum)",
        "[InstrumentLineShapeEstimationFromKeypointDistance][IntegrationTest][LongRunningTest]")
    {
        std::vector<std::pair<std::string, double>> noCrossSections;
        FraunhoferSpectrumGeneration fraunhoferSpectrumGenerator{ GetSolarAtlasFile(), noCrossSections };

        const auto pixelToWavelengthMapping = GetPixelToWavelengthMappingFromFile(GetInitialPixelToWavelengthCalibration_I2J8549());
        InstrumentLineShapeEstimationFromKeypointDistance sut{ pixelToWavelengthMapping };
        const double actualFwhm = 0.35;
        CSpectrum measuredSpectrum;
        GetMeasuredPakFileSpectrum(GetMeasuredSpectrumName_I2J8549(), measuredSpectrum);

        // Act
        CCrossSectionData estimatedLineShape;
        double estimatedFwhm = 0.0;
        sut.EstimateInstrumentLineShape(fraunhoferSpectrumGenerator, measuredSpectrum, estimatedLineShape, estimatedFwhm);

        // Assert
        REQUIRE(estimatedFwhm == Approx(actualFwhm).margin(0.50 * actualFwhm)); // 50% margin
    }

    TEST_CASE(
        "InstrumentLineShapeEstimationFromKeypointDistance (measured 2009175M1 spectrum)",
        "[InstrumentLineShapeEstimationFromKeypointDistance][IntegrationTest][LongRunningTest]")
    {
        std::vector<std::pair<std::string, double>> noCrossSections;
        FraunhoferSpectrumGeneration fraunhoferSpectrumGenerator{ GetSolarAtlasFile(), noCrossSections };

        // Notice that here it is known that the initial pixel-to-wavelength calibration is quite far off...
        const auto pixelToWavelengthMapping = GetPixelToWavelengthMappingFromFile(GetInitialPixelToWavelengthCalibration_2009175M1());
        InstrumentLineShapeEstimationFromKeypointDistance sut{ pixelToWavelengthMapping };

        const double actualFwhm = 0.48;
        CSpectrum measuredSpectrum;
        GetMeasuredPakFileSpectrum(GetMeasuredSpectrumName_2009175M1(), measuredSpectrum);

        // Act
        CCrossSectionData estimatedLineShape;
        double estimatedFwhm = 0.0;
        sut.EstimateInstrumentLineShape(fraunhoferSpectrumGenerator, measuredSpectrum, estimatedLineShape, estimatedFwhm);

        // Assert
        REQUIRE(estimatedFwhm == Approx(actualFwhm).margin(0.50 * actualFwhm)); // 50% margin
    }

    TEST_CASE(
        "InstrumentLineShapeEstimationFromKeypointDistance (measured MAYP11440 spectrum)",
        "[InstrumentLineShapeEstimationFromKeypointDistance][IntegrationTest][LongRunningTest]")
    {
        std::vector<std::pair<std::string, double>> noCrossSections;
        FraunhoferSpectrumGeneration fraunhoferSpectrumGenerator{ GetSolarAtlasFile(), noCrossSections };

        // Notice that here it is known that the initial pixel-to-wavelength calibration is quite far off...
        const auto pixelToWavelengthMapping = GetPixelToWavelengthMappingFromFile(GetInitialPixelToWavelengthCalibration_MAYP11440());
        InstrumentLineShapeEstimationFromKeypointDistance sut{ pixelToWavelengthMapping };

        const double actualFwhm = 0.42;
        CSpectrum measuredSpectrum;
        GetMeasuredSpectrum(GetSkySpectrumName_MAYP11440(), GetDarkSpectrumName_MAYP11440(), measuredSpectrum);

        // Act
        CCrossSectionData estimatedLineShape;
        double estimatedFwhm = 0.0;
        sut.EstimateInstrumentLineShape(fraunhoferSpectrumGenerator, measuredSpectrum, estimatedLineShape, estimatedFwhm);

        // Assert
        REQUIRE(estimatedFwhm == Approx(actualFwhm).margin(0.50 * actualFwhm)); // 50% margin
    }
}
