#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>
#include <SpectralEvaluation/File/ScanEvaluationLogFileHandler.h>
#include "catch.hpp"
#include "TestData.h"

#include <SpectralEvaluation/VectorUtils.h>

using namespace novac;

struct PlumeMeasurement
{
    std::vector<double> scanAngles;
    std::vector<double> phi;
    std::vector<double> columns;
    std::vector<double> columnErrors;
    std::vector<bool>   badEvaluation;
};

PlumeMeasurement GenerateGaussianPlume(double amplitude, double plumeCenter, double fwhmInScanAngles, double offset = 0.0)
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
        plume.columns[stepIdx] = plume.columns[stepIdx] + offset;

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

TEST_CASE("FindPlume with Gaussian plume", "[PlumeProperties]")
{
    CPlumeInScanProperty plumeProperties;
    const double angleTolerance = 3.6; // one step on the scanner.

    SECTION("No plume - returns false.")
    {
        const double columnOffsetInMeasurement = 0.0;
        PlumeMeasurement plume = GenerateGaussianPlume(0.0, 1.0, 20.0);

        REQUIRE(false == FindPlume(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, (long)plume.scanAngles.size(), columnOffsetInMeasurement, plumeProperties));
    }

    SECTION("Plume overhead - returns true and sets correct plume properties.")
    {
        const double acutalPlumeCenter = 0.0;
        const double actualPlumeFwhm = 20.0;
        const double columnOffsetInMeasurement = 0.0;
        PlumeMeasurement plume = GenerateGaussianPlume(200.0, acutalPlumeCenter, actualPlumeFwhm);
        const double expectedPlumeWidth = CalculateGaussianOneOverEFromFwhm(actualPlumeFwhm);

        REQUIRE(true == FindPlume(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, (long)plume.scanAngles.size(), columnOffsetInMeasurement, plumeProperties));

        // Assert, the calculated values should be as expected
        REQUIRE(std::abs(plumeProperties.plumeCenter) < angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeCenter2) < 1.0);
        REQUIRE(std::abs(plumeProperties.plumeEdgeHigh - acutalPlumeCenter - expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeEdgeLow - acutalPlumeCenter + expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeHalfHigh - acutalPlumeCenter - actualPlumeFwhm * 0.5) < 1.5 * angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeHalfLow - acutalPlumeCenter + actualPlumeFwhm * 0.5) < 1.5 * angleTolerance);
    }

    SECTION("Plume overhead - large offset - returns true and sets correct plume properties.")
    {
        const double acutalPlumeCenter = 0.0;
        const double actualPlumeFwhm = 20.0;
        const double columnOffsetInMeasurement = -200.0;
        PlumeMeasurement plume = GenerateGaussianPlume(200.0, acutalPlumeCenter, actualPlumeFwhm, columnOffsetInMeasurement);
        const double expectedPlumeWidth = CalculateGaussianOneOverEFromFwhm(actualPlumeFwhm);

        REQUIRE(true == FindPlume(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, (long)plume.scanAngles.size(), columnOffsetInMeasurement, plumeProperties));

        // Assert, the calculated values should be as expected
        REQUIRE(std::abs(plumeProperties.plumeCenter) < angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeCenter2) < 1.0);
        REQUIRE(std::abs(plumeProperties.plumeEdgeHigh - acutalPlumeCenter - expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeEdgeLow - acutalPlumeCenter + expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeHalfHigh - acutalPlumeCenter - actualPlumeFwhm * 0.5) < 1.5 * angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeHalfLow - acutalPlumeCenter + actualPlumeFwhm * 0.5) < 1.5 * angleTolerance);
    }

    SECTION("Plume at -45degrees - returns true and sets correct plume properties.")
    {
        const double plumeCenter = -45.0;
        const double columnOffsetInMeasurement = 0.0;
        PlumeMeasurement plume = GenerateGaussianPlume(200.0, plumeCenter, 20.0);
        const double expectedPlumeWidth = CalculateGaussianOneOverEFromFwhm(20.0);

        REQUIRE(true == FindPlume(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, (long)plume.scanAngles.size(), columnOffsetInMeasurement, plumeProperties));

        // Assert, the calculated values should be as expected
        REQUIRE(std::abs(plumeProperties.plumeCenter - plumeCenter) < 4.0);
        REQUIRE(std::abs(plumeProperties.plumeCenter2) < 1.0);
        REQUIRE(std::abs(plumeProperties.plumeEdgeHigh - plumeCenter - expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeEdgeLow - plumeCenter + expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeHalfHigh - plumeCenter - 20.0 * 0.5) < 1.5 * angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeHalfLow - plumeCenter + 20.0 * 0.5) < 1.5 * angleTolerance);
    }

    SECTION("Plume at -45degrees - large offset  - returns true and sets correct plume properties.")
    {
        const double acutalPlumeCenter = -45.0;
        const double actualPlumeFwhm = 20.0;
        const double columnOffsetInMeasurement = -400.0;
        PlumeMeasurement plume = GenerateGaussianPlume(200.0, acutalPlumeCenter, actualPlumeFwhm, columnOffsetInMeasurement);
        const double expectedPlumeWidth = CalculateGaussianOneOverEFromFwhm(20.0);

        REQUIRE(true == FindPlume(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, (long)plume.scanAngles.size(), columnOffsetInMeasurement, plumeProperties));

        // Assert, the calculated values should be as expected
        REQUIRE(std::abs(plumeProperties.plumeCenter - acutalPlumeCenter) < 4.0);
        REQUIRE(std::abs(plumeProperties.plumeCenter2) < 1.0);
        REQUIRE(std::abs(plumeProperties.plumeEdgeHigh - acutalPlumeCenter - expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeEdgeLow - acutalPlumeCenter + expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeHalfHigh - acutalPlumeCenter - 20.0 * 0.5) < 1.5 * angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeHalfLow - acutalPlumeCenter + 20.0 * 0.5) < 1.5 * angleTolerance);
    }

    SECTION("Plume at +45degrees - returns true and sets correct plume properties.")
    {
        const double plumeCenter = 45.0;
        const double columnOffsetInMeasurement = 0.0;
        PlumeMeasurement plume = GenerateGaussianPlume(200.0, plumeCenter, 20.0);
        const double expectedPlumeWidth = CalculateGaussianOneOverEFromFwhm(20.0);

        REQUIRE(true == FindPlume(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, (long)plume.scanAngles.size(), columnOffsetInMeasurement, plumeProperties));

        // Assert, the calculated values should be as expected
        REQUIRE(std::abs(plumeProperties.plumeCenter - plumeCenter) < 4.0);
        REQUIRE(std::abs(plumeProperties.plumeCenter2) < 1.0);
        REQUIRE(std::abs(plumeProperties.plumeEdgeHigh - plumeCenter - expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeEdgeLow - plumeCenter + expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeHalfHigh - plumeCenter - 20.0 * 0.5) < 1.5 * angleTolerance);
        REQUIRE(std::abs(plumeProperties.plumeHalfLow - plumeCenter + 20.0 * 0.5) < 1.5 * angleTolerance);
    }

    SECTION("Plume at -75degrees - returns true and sets correct plume properties.")
    {
        const double acutalPlumeCenter = -75.0;
        const double actualPlumeFwhm = 20.0;
        const double columnOffsetInMeasurement = 0.0;
        PlumeMeasurement plume = GenerateGaussianPlume(200.0, acutalPlumeCenter, actualPlumeFwhm, columnOffsetInMeasurement);

        REQUIRE(true == FindPlume(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, (long)plume.scanAngles.size(), columnOffsetInMeasurement, plumeProperties));

        // Assert, the calculated values should be as expected
        REQUIRE(std::abs(plumeProperties.plumeCenter + 75.0) < 4.0);
        REQUIRE(std::abs(plumeProperties.plumeCenter2) < 1.0);
    }

    SECTION("Plume at +75degrees - returns true and sets correct plume properties.")
    {
        const double columnOffsetInMeasurement = 0.0;
        PlumeMeasurement plume = GenerateGaussianPlume(200.0, +75.0, 20.0);

        REQUIRE(true == FindPlume(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, (long)plume.scanAngles.size(), columnOffsetInMeasurement, plumeProperties));

        // Assert, the calculated values should be as expected
        REQUIRE(std::abs(plumeProperties.plumeCenter - 75.0) < 4.0);
        REQUIRE(std::abs(plumeProperties.plumeCenter2) < 1.0);
    }

    SECTION("Plume at -45degrees - bad values edges of scan - returns true and sets plume properties within scan range.")
    {
        const double acutalPlumeCenter = -45.0;
        const double actualPlumeFwhm = 90.0;
        const double columnOffsetInMeasurement = -100.0;
        PlumeMeasurement plume = GenerateGaussianPlume(200.0, acutalPlumeCenter, actualPlumeFwhm, columnOffsetInMeasurement);
        // add a couple of 'bad' measurements in the beginning..
        plume.columns[0] = -5000.0;
        plume.columns[1] = -5000.0;

        REQUIRE(true == FindPlume(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, (long)plume.scanAngles.size(), columnOffsetInMeasurement, plumeProperties));

        // Assert, the calculated values should be as expected
        REQUIRE(plumeProperties.plumeCenter < plume.scanAngles.back());
        REQUIRE(plumeProperties.plumeCenter > plume.scanAngles.front());
    }
}

// Test labelled as integration test since we do need to read data from file.
TEST_CASE("FindPlume with measured plume (BroSo2 ratio measurement)", "[PlumeProperties][IntegrationTest]")
{
    novac::CScanEvaluationLogFileHandler evaluationFileHandler;
    const bool evaluationFileIsOk = evaluationFileHandler.ReadEvaluationLog(TestData::GetBrORatioEvaluationFile1());
    REQUIRE(evaluationFileIsOk); // check assumption on the setup
    REQUIRE(evaluationFileHandler.m_scan.size() == 1); // check assumption on the setup

    CPlumeInScanProperty plumeProperties;
    const double offset = CalculatePlumeOffset(evaluationFileHandler.m_scan[0], 0, plumeProperties);

    // Act
    REQUIRE(true == FindPlume(evaluationFileHandler.m_scan[0], 0, offset, plumeProperties));

    SECTION("Correct plume center.")
    {
        REQUIRE(plumeProperties.plumeCenter == Approx(-23.12).margin(0.01));
    }

    SECTION("Correct plume edges.")
    {
        // Notice that this scan has some dark spectra in the beginning, hence the plume is not entirely complete.
        REQUIRE(plumeProperties.plumeEdgeLow == Approx(-50.0).margin(0.01));
        REQUIRE(plumeProperties.plumeEdgeHigh == Approx(0.0).margin(0.01));

        REQUIRE(plumeProperties.plumeHalfLow == Approx(-50.0).margin(0.01));
        REQUIRE(plumeProperties.plumeHalfHigh == Approx(0.0).margin(0.01));
    }

    SECTION("Does not set completeness nor offset.")
    {
        REQUIRE(plumeProperties.completeness == Approx(0.0));
        REQUIRE(plumeProperties.offset == Approx(offset));
    }
}


TEST_CASE("CalculatePlumeCompleteness with Gaussian plume", "[PlumeProperties]")
{
    CPlumeInScanProperty plumeProperties;

    SECTION("No plume - returns false and sets completeness.")
    {
        plumeProperties.completeness = -1e19;
        PlumeMeasurement plume = GenerateGaussianPlume(0.0, 1.0, 20.0);

        REQUIRE(false == CalculatePlumeCompleteness(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, 0.0, (long)plume.scanAngles.size(), plumeProperties));
        REQUIRE(plumeProperties.completeness == Approx(0.0));
    }

    SECTION("Plume overhead - returns true and sets completeness to 1.0.")
    {
        const double plumeCenter = 0.0;
        PlumeMeasurement plume = GenerateGaussianPlume(200.0, plumeCenter, 20.0);

        REQUIRE(true == CalculatePlumeCompleteness(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, 0.0, (long)plume.scanAngles.size(), plumeProperties));

        REQUIRE(plumeProperties.completeness == Approx(1.0));
    }

    SECTION("Narrow plume at -45degrees - returns true and sets completeness to 1.0.")
    {
        const double plumeCenter = -45.0;
        PlumeMeasurement plume = GenerateGaussianPlume(200.0, plumeCenter, 20.0);

        REQUIRE(true == CalculatePlumeCompleteness(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, 0.0, (long)plume.scanAngles.size(), plumeProperties));

        // Assert, the plume is relatively narrow and does not extend down to horizon. Hence everything is visible.
        REQUIRE(plumeProperties.completeness == Approx(1.0).margin(0.1));
    }

    SECTION("Narrow plume at -45degree - large negative offset - returns true and sets completeness to 1.0.")
    {
        const double plumeCenter = -45.0;
        const double plumeAmplitude = 200.0;
        const double columnOffset = -400.0;
        PlumeMeasurement plume = GenerateGaussianPlume(plumeAmplitude, plumeCenter, 20.0, columnOffset);

        REQUIRE(true == CalculatePlumeCompleteness(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, 0.0, (long)plume.scanAngles.size(), plumeProperties));

        // Assert, the plume is relatively narrow and does not extend down to horizon. Hence everything is visible.
        REQUIRE(plumeProperties.completeness == Approx(1.0).margin(0.1));
    }

    SECTION("Wide plume at -45degrees - returns true and sets completeness to 1.0.")
    {
        const double plumeCenter = -45.0;
        PlumeMeasurement plume = GenerateGaussianPlume(200.0, plumeCenter, 40.0);

        REQUIRE(true == CalculatePlumeCompleteness(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, 0.0, (long)plume.scanAngles.size(), plumeProperties));

        // Assert, the plume is relatively wide but still does not extend down to horizon. Hence everything is visible.
        REQUIRE(plumeProperties.completeness == Approx(1.0).margin(0.1));
    }

    SECTION("Plume at +45degrees - returns true and sets correct completeness.")
    {
        const double plumeCenter = 45.0;
        PlumeMeasurement plume = GenerateGaussianPlume(200.0, plumeCenter, 20.0);

        REQUIRE(true == CalculatePlumeCompleteness(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, 0.0, (long)plume.scanAngles.size(), plumeProperties));

        // Assert, the plume is relatively narrow and does not extend down to horizon. Hence everything is visible.
        REQUIRE(plumeProperties.completeness == Approx(1.0).margin(0.1));
    }

    SECTION("Wide plume at +45degrees - returns true and sets correct completeness.")
    {
        const double plumeCenter = 45.0;
        const double plumeFwhm = 90.0;
        PlumeMeasurement plume = GenerateGaussianPlume(200.0, plumeCenter, plumeFwhm);

        REQUIRE(true == CalculatePlumeCompleteness(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, 0.0, (long)plume.scanAngles.size(), plumeProperties));

        // Assert, the plume is relatively wide and does extend down to horizon.
        REQUIRE(plumeProperties.completeness == Approx(0.7).margin(0.1));
    }


    SECTION("Plume at -75degrees - returns true and sets correct completeness.")
    {
        PlumeMeasurement plume = GenerateGaussianPlume(200.0, -75.0, 20.0);

        REQUIRE(true == CalculatePlumeCompleteness(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, 0.0, (long)plume.scanAngles.size(), plumeProperties));

        // Assert, the plume is relatively narrow but also at a low angle and does hence extend down to horizon.
        REQUIRE(plumeProperties.completeness == Approx(0.7).margin(0.1));
    }

    SECTION("Plume at +75degrees - returns true and sets correct completeness.")
    {
        PlumeMeasurement plume = GenerateGaussianPlume(200.0, +75.0, 20.0);

        REQUIRE(true == CalculatePlumeCompleteness(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, 0.0, (long)plume.scanAngles.size(), plumeProperties));

        // Assert, the plume is relatively narrow but also at a low angle and does hence extend down to horizon.
        REQUIRE(plumeProperties.completeness == Approx(0.7).margin(0.1));
    }
}

// Test labelled as integration test since we do need to read data from file.
TEST_CASE("CalculatePlumeCompleteness with measured plume (BroSo2 ratio measurement)", "[PlumeProperties][IntegrationTest]")
{
    novac::CScanEvaluationLogFileHandler evaluationFileHandler;
    const bool evaluationFileIsOk = evaluationFileHandler.ReadEvaluationLog(TestData::GetBrORatioEvaluationFile1());
    REQUIRE(evaluationFileIsOk); // check assumption on the setup
    REQUIRE(evaluationFileHandler.m_scan.size() == 1); // check assumption on the setup
    CPlumeInScanProperty plumeProperties;
    CalculatePlumeOffset(evaluationFileHandler.m_scan[0], 0, plumeProperties); // The offset needs to be filled in first, as a preparation

    // Act
    REQUIRE(true == CalculatePlumeCompleteness(evaluationFileHandler.m_scan[0], 0, plumeProperties));

    SECTION("Sets correct completeness.")
    {
        REQUIRE(plumeProperties.completeness == Approx(0.7).margin(0.01));
    }

    SECTION("Sets correct plume center.")
    {
        REQUIRE(plumeProperties.plumeCenter == Approx(-23.12).margin(0.01));
    }

    SECTION("Sets correct plume edges.")
    {
        // Notice that this scan has some dark spectra in the beginning, hence the plume is not entirely complete.
        REQUIRE(plumeProperties.plumeEdgeLow == Approx(-50.0).margin(0.01));
        REQUIRE(plumeProperties.plumeEdgeHigh == Approx(0.0).margin(0.01));

        REQUIRE(plumeProperties.plumeHalfLow == Approx(-50.0).margin(0.01));
        REQUIRE(plumeProperties.plumeHalfHigh == Approx(0.0).margin(0.01));
    }
}


// Test labelled as integration test since we do need to read data from file.
TEST_CASE("CalculatePlumeOffset with measured plume (BroSo2 ratio measurement)", "[PlumeProperties][IntegrationTest]")
{
    novac::CScanEvaluationLogFileHandler evaluationFileHandler;
    const bool evaluationFileIsOk = evaluationFileHandler.ReadEvaluationLog(TestData::GetBrORatioEvaluationFile1());
    REQUIRE(evaluationFileIsOk); // check assumption on the setup
    REQUIRE(evaluationFileHandler.m_scan.size() == 1); // check assumption on the setup
    CPlumeInScanProperty plumeProperties;
    CalculatePlumeOffset(evaluationFileHandler.m_scan[0], 0, plumeProperties); // The offset needs to be filled in first, as a preparation

    // Act
    const double calculatedOffset = CalculatePlumeOffset(evaluationFileHandler.m_scan[0], 0, plumeProperties);

    // Assert
    REQUIRE(calculatedOffset == Approx(-1.2e18).margin(1e17));
    REQUIRE(plumeProperties.offset == Approx(-1.2e18).margin(1e17));
}