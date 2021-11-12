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
double GaussianFwhmToSigma(double fwhm);

TEST_CASE("EstimateInstrumentLineShape with measured spectrum being Gaussian of 0.5nm fwhm - returns correct fitted function", "[InstrumentLineshapeEstimationFromDoas]")
{
    std::vector<double> pixelToWavelengthMapping = GeneratePixelToWavelengthMapping(330.0, 350.0, 0.05);

    InstrumentLineshapeEstimationFromDoas::LineShapeEstimationSettings settings;
    settings.startPixel = pixelToWavelengthMapping.size() / 20; // 5% margin in either end
    settings.endPixel = pixelToWavelengthMapping.size() - pixelToWavelengthMapping.size() / 20;

    std::vector<std::pair<std::string, double>> noCrossSections;
    FraunhoferSpectrumGeneration fraunhoferSpectrumGenerator{ GetSolarAtlasFileName(), noCrossSections };

    const double actualFwhm = 0.5;
    CCrossSectionData actualnstrumentLineShape;
    CreateGaussian(GaussianFwhmToSigma(actualFwhm), 0.05, actualnstrumentLineShape);

    // Create the measured spectrum using the actual Fwhm
    auto measuredSpectrum = fraunhoferSpectrumGenerator.GetFraunhoferSpectrum(pixelToWavelengthMapping, actualnstrumentLineShape);

    SECTION("Correct Estimation of Instrument Line Shape with actual line shape as input")
    {
        // Create the SUT with the correct initial values, should return the same solution back...
        InstrumentLineshapeEstimationFromDoas sut{ pixelToWavelengthMapping, actualnstrumentLineShape };

        // Act
        const auto output = sut.EstimateInstrumentLineShape(fraunhoferSpectrumGenerator, *measuredSpectrum, settings);

        // Assert
        REQUIRE(fabs(actualFwhm - output.result.lineShape.Fwhm()) < 0.01 * actualFwhm); // 1% margin
        REQUIRE(fabs(2.0 - output.result.lineShape.k) < 0.03); // the function must be nearly Gaussian.
        REQUIRE(fabs(output.result.shift) < 0.1); // in pixels
    }

    SECTION("Correct Estimation of Instrument Line Shape with initial guess being Gaussian of 0.3nm fwhm")
    {
        const double initialGuessForFwhm = 0.3;
        CCrossSectionData initialGuessForInstrumentLineShape;
        CreateGaussian(GaussianFwhmToSigma(initialGuessForFwhm), 0.05, initialGuessForInstrumentLineShape);
        InstrumentLineshapeEstimationFromDoas sut{ pixelToWavelengthMapping, initialGuessForInstrumentLineShape };

        // Act
        const auto output = sut.EstimateInstrumentLineShape(fraunhoferSpectrumGenerator, *measuredSpectrum, settings);

        // Assert that the result is nearly gaussian with nearly the correct fwhm.
        // Since the initial fwhm guess is rather far off the tolerance can be quite large.
        REQUIRE(fabs(actualFwhm - output.result.lineShape.Fwhm()) < 0.02 * actualFwhm); // 2% margin
        REQUIRE(fabs(2.0 - output.result.lineShape.k) < 0.06);
        REQUIRE(fabs(output.result.shift) < 0.1); // in pixels
    }

    SECTION("Correct Estimation of Instrument Line Shape with initial guess being Gaussian of 0.7nm fwhm")
    {
        const double initialGuessForFwhm = 0.7;
        CCrossSectionData initialGuessForInstrumentLineShape;
        CreateGaussian(GaussianFwhmToSigma(initialGuessForFwhm), 0.05, initialGuessForInstrumentLineShape);
        InstrumentLineshapeEstimationFromDoas sut{ pixelToWavelengthMapping, initialGuessForInstrumentLineShape };

        // Act
        const auto output = sut.EstimateInstrumentLineShape(fraunhoferSpectrumGenerator, *measuredSpectrum, settings);

        // Assert that the result is nearly gaussian with nearly the correct fwhm.
        // Since the initial fwhm guess is rather far off the tolerance can be quite large.
        REQUIRE(fabs(actualFwhm - output.result.lineShape.Fwhm()) < 0.02 * actualFwhm); // 2% margin
        REQUIRE(fabs(2.0 - output.result.lineShape.k) < 0.04);
        REQUIRE(fabs(output.result.shift) < 0.1); // in pixels
    }
}

TEST_CASE("EstimateInstrumentLineShape with measured spectrum being Super Gaussian of 0.5nm fwhm and k = 2.5 - returns correct fitted function", "[InstrumentLineshapeEstimationFromDoas]")
{
    std::vector<double> pixelToWavelengthMapping = GeneratePixelToWavelengthMapping(330.0, 350.0, 0.05);

    InstrumentLineshapeEstimationFromDoas::LineShapeEstimationSettings settings;
    settings.startPixel = pixelToWavelengthMapping.size() / 20; // 5% margin in either end
    settings.endPixel = pixelToWavelengthMapping.size() - pixelToWavelengthMapping.size() / 20;

    std::vector<std::pair<std::string, double>> noCrossSections;
    FraunhoferSpectrumGeneration fraunhoferSpectrumGenerator{ GetSolarAtlasFileName(), noCrossSections };

    const double actualFwhm = 0.5;
    SuperGaussianLineShape actualParameterizedLineShape;
    actualParameterizedLineShape.k = 2.5;
    actualParameterizedLineShape.w = 0.5 * actualFwhm / std::pow(0.69314718056, 1.0 / actualParameterizedLineShape.k);
    REQUIRE(std::abs(actualFwhm - actualParameterizedLineShape.Fwhm()) < 0.01); // verify that the calculation above is correct

    CCrossSectionData actualnstrumentLineShape = SampleInstrumentLineShape(actualParameterizedLineShape);

    // Create the measured spectrum using the actual Fwhm
    auto measuredSpectrum = fraunhoferSpectrumGenerator.GetFraunhoferSpectrum(pixelToWavelengthMapping, actualnstrumentLineShape);

    SECTION("Correct Estimation of Instrument Line Shape with actual line shape as input")
    {
        // Create the SUT with the correct initial values, should return the same solution back...
        InstrumentLineshapeEstimationFromDoas sut{ pixelToWavelengthMapping, actualnstrumentLineShape };

        // Act
        const auto output = sut.EstimateInstrumentLineShape(fraunhoferSpectrumGenerator, *measuredSpectrum, settings);

        // Assert
        REQUIRE(fabs(actualFwhm - output.result.lineShape.Fwhm()) < 0.01 * actualFwhm); // 1% margin
        REQUIRE(fabs(actualParameterizedLineShape.k - output.result.lineShape.k) < 0.03); // The power must be nearly correct.
        REQUIRE(fabs(output.result.shift) < 0.1); // in pixels
    }

    SECTION("Correct Estimation of Instrument Line Shape with initial guess being Gaussian of 0.3nm fwhm")
    {
        const double initialGuessForFwhm = 0.3;
        CCrossSectionData initialGuessForInstrumentLineShape;
        CreateGaussian(GaussianFwhmToSigma(initialGuessForFwhm), 0.05, initialGuessForInstrumentLineShape);
        InstrumentLineshapeEstimationFromDoas sut{ pixelToWavelengthMapping, initialGuessForInstrumentLineShape };

        // Act
        const auto output = sut.EstimateInstrumentLineShape(fraunhoferSpectrumGenerator, *measuredSpectrum, settings);

        // Assert that the result has nearly the correct power and nearly the correct fwhm.
        // Since the initial fwhm guess is rather far off the tolerance can be quite large.
        REQUIRE(fabs(actualFwhm - output.result.lineShape.Fwhm()) < 0.02 * actualFwhm); // 2% margin
        REQUIRE(fabs(actualParameterizedLineShape.k - output.result.lineShape.k) < 0.06);
        REQUIRE(fabs(output.result.shift) < 0.1); // in pixels
    }

    SECTION("Correct Estimation of Instrument Line Shape with initial guess being Gaussian of 0.7nm fwhm")
    {
        const double initialGuessForFwhm = 0.7;
        CCrossSectionData initialGuessForInstrumentLineShape;
        CreateGaussian(GaussianFwhmToSigma(initialGuessForFwhm), 0.05, initialGuessForInstrumentLineShape);
        InstrumentLineshapeEstimationFromDoas sut{ pixelToWavelengthMapping, initialGuessForInstrumentLineShape };

        // Act
        const auto output = sut.EstimateInstrumentLineShape(fraunhoferSpectrumGenerator, *measuredSpectrum, settings);

        // Assert that the result has nearly the correct power and nearly the correct fwhm.
        // Since the initial fwhm guess is rather far off the tolerance can be quite large.
        REQUIRE(fabs(actualFwhm - output.result.lineShape.Fwhm()) < 0.02 * actualFwhm); // 2% margin
        REQUIRE(fabs(actualParameterizedLineShape.k - output.result.lineShape.k) < 0.04);
        REQUIRE(fabs(output.result.shift) < 0.1); // in pixels
    }
}