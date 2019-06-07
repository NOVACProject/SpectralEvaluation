#include "catch.hpp"
#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>

struct PlumeMeasurement
{
    std::vector<double> scanAngles;
    std::vector<double> phi;
    std::vector<double> columns;
    std::vector<double> columnErrors;
    std::vector<bool>   badEvaluation;
};

PlumeMeasurement GenerateGaussianPlume(double amplitude, double plumeCenter, double fwhmInScanAngles)
{
    PlumeMeasurement plume;

    const int length = 51;
    const double sigma2 = (fwhmInScanAngles / 2.3548) * (fwhmInScanAngles / 2.3548); // sigma-squared

    // A scan always has 51 measurements (novac standard)
    plume.scanAngles.resize(length);
    plume.phi.resize(length);
    plume.columns.resize(length);
    plume.columnErrors.resize(length);
    plume.badEvaluation.resize(length);

    const double scanStepSize = 180.0 / (double)(length - 1);
    for (int stepIdx = 0; stepIdx < length; ++stepIdx)
    {
        const double alpha = -90.0 + scanStepSize * stepIdx;
        plume.scanAngles[stepIdx] = alpha;
        plume.phi[stepIdx] = 0.0; // don't simulate the Heidelberg instrument

        plume.columns[stepIdx] = amplitude * exp(-(alpha - plumeCenter) * (alpha - plumeCenter) / (2.0 * sigma2));
        plume.columns[stepIdx] = 0.05 * plume.columns[stepIdx];

        plume.badEvaluation[stepIdx] = false;
    }

    return plume;
}

double CalculateGaussianOneOverEFromFwhm(double gaussianFwhm)
{
    const double sigma = (gaussianFwhm / 2.3548);

    // 1/e occurrs when x^2 = 2*sigma^2
    const double x = sigma * std::sqrt(2.0);

    return 2 * x; // full width.
}

TEST_CASE("FindPlume", "[PlumeProperties]")
{
    CPlumeInScanProperty plumeProperties;
    const double angleTolerance = 3.6; // one step on the scanner.

    SECTION("No plume - returns false.")
    {
        PlumeMeasurement plume = GenerateGaussianPlume(0.0, 1.0, 20.0);
        
        REQUIRE(false == FindPlume(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, (long)plume.scanAngles.size(), plumeProperties));
    }

    SECTION("Plume overhead - returns true and sets correct plume-center and edges.")
    {
        const double plumeCenter = 0.0;
        PlumeMeasurement plume = GenerateGaussianPlume(200.0, plumeCenter, 20.0);
        const double expectedPlumeWidth = CalculateGaussianOneOverEFromFwhm(20.0);

        REQUIRE(true == FindPlume(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, (long)plume.scanAngles.size(), plumeProperties));
        REQUIRE(std::abs(plumeProperties.plumeCenter) < angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeCenter2) < 1.0);
        REQUIRE(std::abs(plumeProperties.plumeEdgeHigh - plumeCenter - expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeEdgeLow - plumeCenter + expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeHalfHigh - plumeCenter - 20.0 * 0.5) < 1.5 * angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeHalfLow - plumeCenter + 20.0 * 0.5) < 1.5 * angleTolerance);
    }

    SECTION("Plume at -45degrees - returns true and sets correct plume-center.")
    {
        const double plumeCenter = -45.0;
        PlumeMeasurement plume = GenerateGaussianPlume(200.0, plumeCenter, 20.0);
        const double expectedPlumeWidth = CalculateGaussianOneOverEFromFwhm(20.0);

        REQUIRE(true == FindPlume(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, (long)plume.scanAngles.size(), plumeProperties));
        REQUIRE(std::abs(plumeProperties.plumeCenter - plumeCenter) < 4.0);
        REQUIRE(std::abs(plumeProperties.plumeCenter2) < 1.0);
        REQUIRE(std::abs(plumeProperties.plumeEdgeHigh - plumeCenter - expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeEdgeLow - plumeCenter + expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeHalfHigh - plumeCenter - 20.0 * 0.5) < 1.5 * angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeHalfLow - plumeCenter + 20.0 * 0.5) < 1.5 * angleTolerance);
    }

    SECTION("Plume at +45degrees - returns true and sets correct plume-center.")
    {
        const double plumeCenter = 45.0;
        PlumeMeasurement plume = GenerateGaussianPlume(200.0, plumeCenter, 20.0);
        const double expectedPlumeWidth = CalculateGaussianOneOverEFromFwhm(20.0);

        REQUIRE(true == FindPlume(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, (long)plume.scanAngles.size(), plumeProperties));
        REQUIRE(std::abs(plumeProperties.plumeCenter - plumeCenter) < 4.0);
        REQUIRE(std::abs(plumeProperties.plumeCenter2) < 1.0);
        REQUIRE(std::abs(plumeProperties.plumeEdgeHigh - plumeCenter - expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeEdgeLow - plumeCenter + expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeHalfHigh - plumeCenter - 20.0 * 0.5) < 1.5 * angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeHalfLow - plumeCenter + 20.0 * 0.5) < 1.5 * angleTolerance);
    }

    SECTION("Plume at -75degrees - returns true and sets correct plume-center.")
    {
        PlumeMeasurement plume = GenerateGaussianPlume(200.0, -75.0, 20.0);

        REQUIRE(true == FindPlume(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, (long)plume.scanAngles.size(), plumeProperties));
        REQUIRE(std::abs(plumeProperties.plumeCenter + 75.0) < 4.0);
        REQUIRE(std::abs(plumeProperties.plumeCenter2) < 1.0);
    }

    SECTION("Plume at +75degrees - returns true and sets correct plume-center.")
    {
        PlumeMeasurement plume = GenerateGaussianPlume(200.0, +75.0, 20.0);

        REQUIRE(true == FindPlume(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, (long)plume.scanAngles.size(), plumeProperties));
        REQUIRE(std::abs(plumeProperties.plumeCenter - 75.0) < 4.0);
        REQUIRE(std::abs(plumeProperties.plumeCenter2) < 1.0);
    }


    // TODO: Continue here...

}