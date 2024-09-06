#include "catch.hpp"
#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>

using namespace novac;

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

    SECTION("Actual measurement #1")
    {
        const int numPoints = 44;
        std::vector<double> columns = { 17.1010, -0.0151, -18.752, 0.0077, -1.376, -6.768, 1.322, 5.878, 16.303, 13.216, 18.263, 18.570, 18.566, 14.322, 10.781, 14.349, 10.949, 8.482, 3.636, 0.778, 1.369, 1.440, -1.00, -4.751, -10.700, -12.756, -13.213, -8.257, -18.936, -15.519, -19.413, -20.048, -26.563, -24.695, -29.760, -30.331, -29.405, -34.421, -40.300, -39.843, -48.137, -52.413, -52.963, -47.123 };
        std::vector<double> columnErrors = { 7.86, 4.57, 4.43, 4.30, 4.55, 4.25, 4.01, 4.11, 4.25, 3.83, 4.24, 4.29,4.45, 4.46, 3.83, 3.92, 3.76, 4.09, 3.76, 3.46, 3.46, 3.69, 3.82, 3.72, 3.75, 3.78, 4.05, 3.78, 3.80, 4.01, 4.28, 4.27, 4.26, 4.00, 4.55, 4.21, 4.58, 4.03, 4.53, 4.53, 4.61, 4.94, 4.73, 4.05 };

        std::vector<double> phi(numPoints);
        std::fill(begin(phi), end(phi), 0.0);

        std::vector<bool> badEvaluation(numPoints);
        std::fill(begin(badEvaluation), end(badEvaluation), false);

        std::vector<double> scanAngles(numPoints);
        for (int idx = 0; idx < numPoints; ++idx)
        {
            scanAngles[idx] = -90.0 + idx * 3.0;
        }

        REQUIRE(true == FindPlume(scanAngles, phi, columns, columnErrors, badEvaluation, numPoints, plumeProperties));
        REQUIRE(std::abs(plumeProperties.plumeCenter - 75.0) < 4.0);
        REQUIRE(std::abs(plumeProperties.plumeCenter2) < 1.0);
    }

    // TODO: Continue here...
}