#include "catch.hpp"
#include <string.h>
#include <SpectralEvaluation/Calibration/InstrumentLineShapeEstimation.h>
#include <SpectralEvaluation/Calibration/FraunhoferSpectrumGeneration.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>

namespace novac
{

    std::vector<double> GeneratePixelToWavelengthMapping(double lambdaMin, double lambdaMax, double lambdaStep)
    {
        std::vector<double> result;
        for (double lambda = lambdaMin; lambda < lambdaMax; lambda += lambdaStep)
        {
            result.push_back(lambda);
        }
        return result;
    }

    double GaussianFwhmToSigma(double fwhm)
    {
        return fwhm / (2.0 * std::sqrt(2.0 * std::log(2.0)));
    }

    std::string GetSolarAtlasFileName()
    {
#ifdef _MSC_VER
        return std::string("../TestData/SOLARFL_330-350nm.xs");
#else
        return std::string("TestData/SOLARFL_330-350nm.xs");
#endif // _MSC_VER
    }

    TEST_CASE("GetFwhm returns correct value", "[InstrumentLineShapeEstimationFromKeypointDistance]")
    {
        SECTION("Correct Fwhm of ApproximateGaussian with actual fwhm of 0.5 nm")
        {
            const double actualFwhm = 0.5;
            const double gaussianSigma = GaussianFwhmToSigma(actualFwhm);
            CCrossSectionData actualnstrumentLineShape;
            CreateGaussian(gaussianSigma, 0.1 * actualFwhm, actualnstrumentLineShape);

            // Act
            const double estimatedFwhm = GetFwhm(actualnstrumentLineShape);

            // Assert
            REQUIRE(fabs(actualFwhm - estimatedFwhm) < 0.1 * actualFwhm); // 1% margin
        }

        SECTION("Correct Fwhm of ApproximateGaussian with actual fwhm of 0.01 nm")
        {
            const double actualFwhm = 0.01;
            const double gaussianSigma = GaussianFwhmToSigma(actualFwhm);
            CCrossSectionData actualnstrumentLineShape;
            CreateGaussian(gaussianSigma, 0.1 * actualFwhm, actualnstrumentLineShape);

            // Act
            const double estimatedFwhm = GetFwhm(actualnstrumentLineShape);

            // Assert
            REQUIRE(fabs(actualFwhm - estimatedFwhm) < 0.1 * actualFwhm); // 1% margin
        }

        SECTION("Correct Fwhm of badly sampled ApproximateGaussian with actual fwhm of 0.5 nm")
        {
            const double actualFwhm = 0.5;
            const double gaussianSigma = GaussianFwhmToSigma(actualFwhm);
            CCrossSectionData actualnstrumentLineShape;
            CreateGaussian(gaussianSigma, actualFwhm, actualnstrumentLineShape);

            // Act
            const double estimatedFwhm = GetFwhm(actualnstrumentLineShape);

            // Assert
            REQUIRE(estimatedFwhm == Approx(actualFwhm).margin(0.20 * actualFwhm)); // 20% margin
        }
    }

    TEST_CASE("InstrumentLineShapeEstimationFromKeypointDistance Correct result when actual ILS is gaussian with fwhm 0,5 nm",
        "[InstrumentLineShapeEstimationFromKeypointDistance][UnitTest]")
    {
        std::vector<double> pixelToWavelengthMapping = GeneratePixelToWavelengthMapping(330.0, 350.0, 0.05);
        InstrumentLineShapeEstimationFromKeypointDistance sut{ pixelToWavelengthMapping };

        std::vector<std::pair<std::string, double>> noCrossSections;
        FraunhoferSpectrumGeneration fraunhoferSpectrumGenerator{ GetSolarAtlasFileName(), noCrossSections };

        const double actualFwhm = 0.5;
        const double gaussianSigma = GaussianFwhmToSigma(actualFwhm);
        CCrossSectionData actualnstrumentLineShape;
        CreateGaussian(gaussianSigma, 0.05, actualnstrumentLineShape);
        auto measuredSpectrum = fraunhoferSpectrumGenerator.GetFraunhoferSpectrum(pixelToWavelengthMapping, actualnstrumentLineShape);

        // Act
        CCrossSectionData estimatedLineShape;
        double estimatedFwhm = 0.0;
        sut.EstimateInstrumentLineShape(fraunhoferSpectrumGenerator, *measuredSpectrum, estimatedLineShape, estimatedFwhm);

        // Assert
        REQUIRE(estimatedFwhm == Approx(actualFwhm).margin(0.20 * actualFwhm)); // 20% margin
    }

    TEST_CASE("InstrumentLineShapeEstimationFromKeypointDistance Correct result when actual instrument line shape is gaussian with fwhm 0,4 nm",
        "[InstrumentLineShapeEstimationFromKeypointDistance]")
    {
        std::vector<double> pixelToWavelengthMapping = GeneratePixelToWavelengthMapping(330.0, 350.0, 0.05);
        InstrumentLineShapeEstimationFromKeypointDistance sut{ pixelToWavelengthMapping };

        std::vector<std::pair<std::string, double>> noCrossSections;
        FraunhoferSpectrumGeneration fraunhoferSpectrumGenerator{ GetSolarAtlasFileName(), noCrossSections };

        const double actualFwhm = 0.4;
        const double gaussianSigma = GaussianFwhmToSigma(actualFwhm);
        CCrossSectionData actualnstrumentLineShape;
        CreateGaussian(gaussianSigma, 0.05, actualnstrumentLineShape);
        auto measuredSpectrum = fraunhoferSpectrumGenerator.GetFraunhoferSpectrum(pixelToWavelengthMapping, actualnstrumentLineShape);

        // Act
        CCrossSectionData estimatedLineShape;
        double estimatedFwhm = 0.0;
        sut.EstimateInstrumentLineShape(fraunhoferSpectrumGenerator, *measuredSpectrum, estimatedLineShape, estimatedFwhm);

        // Assert
        REQUIRE(estimatedFwhm == Approx(actualFwhm).margin(0.20 * actualFwhm)); // 20% margin
    }

    TEST_CASE("InstrumentLineShapeEstimationFromKeypointDistance Retrieves correct result when actual instrument line shape is gaussian with fwhm 0,8 nm",
        "[InstrumentLineShapeEstimationFromKeypointDistance]")
    {
        std::vector<double> pixelToWavelengthMapping = GeneratePixelToWavelengthMapping(330.0, 350.0, 0.05);
        InstrumentLineShapeEstimationFromKeypointDistance sut{ pixelToWavelengthMapping };

        std::vector<std::pair<std::string, double>> noCrossSections;
        FraunhoferSpectrumGeneration fraunhoferSpectrumGenerator{ GetSolarAtlasFileName(), noCrossSections };

        const double actualFwhm = 0.8;
        const double gaussianSigma = GaussianFwhmToSigma(actualFwhm);
        CCrossSectionData actualnstrumentLineShape;
        CreateGaussian(gaussianSigma, 0.1, actualnstrumentLineShape);
        auto measuredSpectrum = fraunhoferSpectrumGenerator.GetFraunhoferSpectrum(pixelToWavelengthMapping, actualnstrumentLineShape);

        // Act
        CCrossSectionData estimatedLineShape;
        double estimatedFwhm = 0.0;
        sut.EstimateInstrumentLineShape(fraunhoferSpectrumGenerator, *measuredSpectrum, estimatedLineShape, estimatedFwhm);

        // Assert
        REQUIRE(estimatedFwhm == Approx(actualFwhm).margin(0.20 * actualFwhm)); // 20% margin
    }
}
