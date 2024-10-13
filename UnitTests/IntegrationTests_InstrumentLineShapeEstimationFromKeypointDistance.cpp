#include "catch.hpp"
#include "TestData.h"
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
static void GetMeasuredPakFileSpectrum(const std::string& fileName, CSpectrum& measuredSpectrum)
{
    novac::ConsoleLog log;
    novac::LogContext context("file", fileName);
    CScanFileHandler pakFileHandler(log);
    if (!pakFileHandler.CheckScanFile(context, fileName))
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
    FraunhoferSpectrumGeneration fraunhoferSpectrumGenerator{ TestData::GetSolarAtlasFile(), noCrossSections };

    const auto pixelToWavelengthMapping = GetPixelToWavelengthMappingFromFile(TestData::GetInitialPixelToWavelengthCalibration_I2J8549());
    InstrumentLineShapeEstimationFromKeypointDistance sut{ pixelToWavelengthMapping };
    const double actualFwhm = 0.35;
    CSpectrum measuredSpectrum;
    GetMeasuredPakFileSpectrum(TestData::GetMeasuredSpectrumName_I2J8549(), measuredSpectrum);

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
    FraunhoferSpectrumGeneration fraunhoferSpectrumGenerator{ TestData::GetSolarAtlasFile(), noCrossSections };

    // Notice that here it is known that the initial pixel-to-wavelength calibration is quite far off...
    const auto pixelToWavelengthMapping = GetPixelToWavelengthMappingFromFile(TestData::GetInitialPixelToWavelengthCalibration_2009175M1());
    InstrumentLineShapeEstimationFromKeypointDistance sut{ pixelToWavelengthMapping };

    const double actualFwhm = 0.48;
    CSpectrum measuredSpectrum;
    GetMeasuredPakFileSpectrum(TestData::GetMeasuredSpectrumName_2009175M1(), measuredSpectrum);

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
    FraunhoferSpectrumGeneration fraunhoferSpectrumGenerator{ TestData::GetSolarAtlasFile(), noCrossSections };

    // Notice that here it is known that the initial pixel-to-wavelength calibration is quite far off...
    const auto pixelToWavelengthMapping = GetPixelToWavelengthMappingFromFile(TestData::GetInitialPixelToWavelengthCalibration_MAYP11440());
    InstrumentLineShapeEstimationFromKeypointDistance sut{ pixelToWavelengthMapping };

    const double actualFwhm = 0.42;
    CSpectrum measuredSpectrum;
    GetMeasuredSpectrum(TestData::GetSkySpectrumName_MAYP11440(), TestData::GetDarkSpectrumName_MAYP11440(), measuredSpectrum);

    // Act
    CCrossSectionData estimatedLineShape;
    double estimatedFwhm = 0.0;
    sut.EstimateInstrumentLineShape(fraunhoferSpectrumGenerator, measuredSpectrum, estimatedLineShape, estimatedFwhm);

    // Assert
    REQUIRE(estimatedFwhm == Approx(actualFwhm).margin(0.50 * actualFwhm)); // 50% margin
}
}
