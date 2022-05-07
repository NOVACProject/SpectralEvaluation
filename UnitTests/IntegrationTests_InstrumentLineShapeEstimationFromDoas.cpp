#include "catch.hpp"
#include "TestData.h"
#include <string.h>
#include <SpectralEvaluation/Calibration/InstrumentLineShapeEstimation.h>
#include <SpectralEvaluation/Calibration/FraunhoferSpectrumGeneration.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Spectra/SpectrumUtils.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/File/STDFile.h>

namespace novac
{
    std::vector<double> GetPixelToWavelengthMappingFromFile(const std::string& clbFile);

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
        "InstrumentLineshapeEstimationFromDoas (measured I2J8549 spectrum)",
        "[InstrumentLineshapeEstimationFromDoas][IntegrationTest][LongRunningTest]")
    {
        std::vector<std::pair<std::string, double>> noCrossSections;
        FraunhoferSpectrumGeneration fraunhoferSpectrumGenerator{ TestData::GetSolarAtlasFile(), noCrossSections };
        const auto pixelToWavelengthMapping = GetPixelToWavelengthMappingFromFile(TestData::GetInitialPixelToWavelengthCalibration_I2J8549());
        novac::SuperGaussianLineShape initialLineshape{ 0.3, 2.0 }; // start with a pure Gaussian
        InstrumentLineshapeEstimationFromDoas sut{ pixelToWavelengthMapping, initialLineshape };
        InstrumentLineshapeEstimationFromDoas::LineShapeEstimationSettings settings;
        settings.startPixel = WavelengthToPixel(pixelToWavelengthMapping, 330.0);
        settings.endPixel = WavelengthToPixel(pixelToWavelengthMapping, 350.0);
        const double actualFwhm = 0.50;
        CSpectrum measuredSpectrum;
        GetMeasuredPakFileSpectrum(TestData::GetMeasuredSpectrumName_I2J8549(), measuredSpectrum);

        // Act
        const auto result = sut.EstimateInstrumentLineShape(measuredSpectrum, settings, fraunhoferSpectrumGenerator);

        // Assert
        REQUIRE(result.result.lineShape.Fwhm() == Approx(actualFwhm).margin(0.10 * actualFwhm)); // 20% margin
    }

    TEST_CASE(
        "InstrumentLineshapeEstimationFromDoas (measured 2009175M1 spectrum)",
        "[InstrumentLineshapeEstimationFromDoas][IntegrationTest][LongRunningTest]")
    {
        std::vector<std::pair<std::string, double>> noCrossSections;
        FraunhoferSpectrumGeneration fraunhoferSpectrumGenerator{ TestData::GetSolarAtlasFile(), noCrossSections };

        // Notice that here it is known that the initial pixel-to-wavelength calibration is quite far off...
        const auto pixelToWavelengthMapping = GetPixelToWavelengthMappingFromFile(TestData::GetInitialPixelToWavelengthCalibration_2009175M1());
        novac::SuperGaussianLineShape initialLineshape{ 0.3, 2.0 }; // start with a pure Gaussian
        InstrumentLineshapeEstimationFromDoas sut{ pixelToWavelengthMapping, initialLineshape };
        InstrumentLineshapeEstimationFromDoas::LineShapeEstimationSettings settings;
        settings.startPixel = WavelengthToPixel(pixelToWavelengthMapping, 330.0);
        settings.endPixel = WavelengthToPixel(pixelToWavelengthMapping, 350.0);

        const double actualFwhm = 0.5;
        CSpectrum measuredSpectrum;
        GetMeasuredPakFileSpectrum(TestData::GetMeasuredSpectrumName_2009175M1(), measuredSpectrum);

        // Act
        const auto result = sut.EstimateInstrumentLineShape(measuredSpectrum, settings, fraunhoferSpectrumGenerator);

        // Assert
        REQUIRE(result.result.lineShape.Fwhm() == Approx(actualFwhm).margin(0.15 * actualFwhm)); // 15% margin
    }

    TEST_CASE(
        "InstrumentLineshapeEstimationFromDoas (measured MAYP11440 spectrum)",
        "[InstrumentLineshapeEstimationFromDoas][IntegrationTest][LongRunningTest]")
    {
        std::vector<std::pair<std::string, double>> noCrossSections;
        FraunhoferSpectrumGeneration fraunhoferSpectrumGenerator{ TestData::GetSolarAtlasFile(), noCrossSections };

        // Notice that here it is known that the initial pixel-to-wavelength calibration is quite far off...
        const auto pixelToWavelengthMapping = GetPixelToWavelengthMappingFromFile(TestData::GetInitialPixelToWavelengthCalibration_MAYP11440());
        novac::SuperGaussianLineShape initialLineshape{ 0.3, 2.0 }; // start with a pure Gaussian
        InstrumentLineshapeEstimationFromDoas sut{ pixelToWavelengthMapping, initialLineshape };
        InstrumentLineshapeEstimationFromDoas::LineShapeEstimationSettings settings;
        settings.startPixel = WavelengthToPixel(pixelToWavelengthMapping, 330.0);
        settings.endPixel = WavelengthToPixel(pixelToWavelengthMapping, 350.0);

        const double actualFwhm = 0.42;
        CSpectrum measuredSpectrum;
        GetMeasuredSpectrum(TestData::GetSkySpectrumName_MAYP11440(), TestData::GetDarkSpectrumName_MAYP11440(), measuredSpectrum);

        // Act
        const auto result = sut.EstimateInstrumentLineShape(measuredSpectrum, settings, fraunhoferSpectrumGenerator);

        // Assert. Here we have used a 30% margin, instead of the 20% above as the results of this spectrum is 
        //  slightly worse results than the results above for some reason
        REQUIRE(result.result.lineShape.Fwhm() == Approx(actualFwhm).margin(0.30 * actualFwhm)); // 30% margin
    }
}
