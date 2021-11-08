#include "catch.hpp"
#include <string.h>
#include <SpectralEvaluation/Calibration/InstrumentLineShapeEstimation.h>
#include <SpectralEvaluation/Calibration/FraunhoferSpectrumGeneration.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>

using namespace novac;

// Functions defined elsewhere...
std::vector<double> GeneratePixelToWavelengthMapping(double lambdaMin, double lambdaMax, double lambdaStep);
std::string GetSolarAtlasFileName();

static double GaussianFwhmToSigma(double fwhm)
{
    return fwhm / (2.0 * std::sqrt(2.0 * std::log(2.0)));
}


TEST_CASE("EstimateInstrumentLineShape with measured spectrum being Gaussian of 0.5nm fwhm - returns correct fitted function", "[InstrumentLineshapeEstimationFromDoas]")
{
    std::vector<double> pixelToWavelengthMapping = GeneratePixelToWavelengthMapping(330.0, 350.0, 0.05);

    std::vector<std::pair<std::string, double>> noCrossSections;
    FraunhoferSpectrumGeneration fraunhoferSpectrumGenerator{ GetSolarAtlasFileName(), noCrossSections };

    const double actualFwhm = 0.5;
    CCrossSectionData actualnstrumentLineShape;
    CreateGaussian(GaussianFwhmToSigma(actualFwhm), 0.1, actualnstrumentLineShape);

    // Create the measured spectrum using the actual Fwhm
    auto measuredSpectrum = fraunhoferSpectrumGenerator.GetFraunhoferSpectrum(pixelToWavelengthMapping, actualnstrumentLineShape);

    SECTION("Correct Estimation of Instrument Line Shape with actual line shape as input")
    {
        // Create the SUT with the correct initial values, should return the same solution back...
        InstrumentLineshapeEstimationFromDoas sut{ pixelToWavelengthMapping, actualnstrumentLineShape };

        // Act
        InstrumentLineshapeEstimationFromDoas::LineShapeEstimationSettings settings;
        settings.startPixel = 10;
        settings.endPixel = pixelToWavelengthMapping.size() - 10;
        const auto output = sut.EstimateInstrumentLineShape(fraunhoferSpectrumGenerator, *measuredSpectrum, settings);

        // Assert
        REQUIRE(fabs(actualFwhm - output.result.lineShape.Fwhm()) < 0.01 * actualFwhm); // 1% margin
        REQUIRE(fabs(output.result.shift) < 0.01); // in pixels
    }

    SECTION("Correct Estimation of Instrument Line Shape with initial guess being Gaussian of 0.4nm fwhm")
    {
        const double initialGuessForFwhm = 0.4;
        CCrossSectionData initialGuessForInstrumentLineShape;
        CreateGaussian(GaussianFwhmToSigma(initialGuessForFwhm), 0.1, initialGuessForInstrumentLineShape);
        InstrumentLineshapeEstimationFromDoas sut{ pixelToWavelengthMapping, initialGuessForInstrumentLineShape };

        // Act
        InstrumentLineshapeEstimationFromDoas::LineShapeEstimationSettings settings;
        settings.startPixel = 10;
        settings.endPixel = pixelToWavelengthMapping.size() - 10;
        const auto output = sut.EstimateInstrumentLineShape(fraunhoferSpectrumGenerator, *measuredSpectrum, settings);

        // Assert
        REQUIRE(fabs(actualFwhm - output.result.lineShape.Fwhm()) < 0.01 * actualFwhm); // 1% margin
        REQUIRE(fabs(output.result.shift) < 0.01); // in pixels
    }

}