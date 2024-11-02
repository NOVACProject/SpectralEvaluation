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

static PlumeMeasurement GenerateGaussianPlume(double amplitude, double plumeCenter, double fwhmInScanAngles, double offset = 0.0)
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

static BasicScanEvaluationResult GenerateGaussianScanResult(double amplitude, double plumeCenter, double fwhmInScanAngles, double offset = 0.0)
{
    BasicScanEvaluationResult plume;

    // A scan always has 51 measurements (novac standard)
    const int length = 51;
    const double sigma2 = (fwhmInScanAngles / 2.3548) * (fwhmInScanAngles / 2.3548); // sigma-squared

    const double scanStepSize = 180.0 / (double)(length - 1);
    for (int stepIdx = 0; stepIdx < length; ++stepIdx)
    {
        const double alpha = -90.0 + scanStepSize * stepIdx;
        
        CSpectrumInfo specInfo;
        specInfo.m_scanAngle = static_cast<float>(alpha);
        specInfo.m_scanAngle2 = 0.0; // don't simulate the Heidelberg instrument

        CEvaluationResult evaluation;
        const double column = offset + amplitude * exp(-(alpha - plumeCenter) * (alpha - plumeCenter) / (2.0 * sigma2));

        novac::CReferenceFitResult result;
        result.m_specieName = "SO2";
        result.m_column = column;
        result.m_columnError = 0.0;
        evaluation.m_referenceResult.push_back(result);

        plume.AppendResult(evaluation, specInfo);
    }

    return plume;
}

static double CalculateGaussianOneOverEFromFwhm(double gaussianFwhm)
{
    const double sigma = (gaussianFwhm / 2.3548);

    // 1/e occurrs when x^2 = 2*sigma^2
    const double x = sigma * std::sqrt(2.0);

    return 2 * x; // full width.
}

TEST_CASE("CalculatePlumeProperties with Gaussian plume", "[PlumeProperties]")
{
    const double angleTolerance = 3.6; // one step on the scanner.
    std::string message;

    SECTION("No plume - returns null.")
    {
        auto plume = GenerateGaussianScanResult(0.0, 1.0, 20.0);

        // Act
        auto plumeProperties = CalculatePlumeProperties(plume, StandardMolecule::SO2, message);

        // Assert
        REQUIRE(plumeProperties == nullptr);
    }

    SECTION("Plume overhead - returns correct plume properties.")
    {
        const double acutalPlumeCenter = 0.0;
        const double actualPlumeFwhm = 20.0;
        auto plume = GenerateGaussianScanResult(200.0, acutalPlumeCenter, actualPlumeFwhm);
        REQUIRE(plume.NumberOfEvaluatedSpectra() == 51);
        REQUIRE(0 == plume.GetSpecieIndex(StandardMolecule::SO2));
        const double expectedPlumeWidth = CalculateGaussianOneOverEFromFwhm(actualPlumeFwhm);

        // Act
        auto plumeProperties = CalculatePlumeProperties(plume, StandardMolecule::SO2, message);

        // Assert, the calculated values should be as expected
        REQUIRE(nullptr != plumeProperties);
        REQUIRE(std::abs(plumeProperties->plumeCenter.Value()) < angleTolerance);
        REQUIRE(std::abs(plumeProperties->plumeCenter2.Value()) < angleTolerance); //plumecenter2 is set to zero, all the scan2 angles are zero
        REQUIRE(std::abs(plumeProperties->plumeEdgeHigh.Value() - acutalPlumeCenter - expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties->plumeEdgeLow.Value() - acutalPlumeCenter + expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties->plumeHalfHigh.Value() - acutalPlumeCenter - actualPlumeFwhm * 0.5) < 1.5 * angleTolerance);
        REQUIRE(std::abs(plumeProperties->plumeHalfLow.Value() - acutalPlumeCenter + actualPlumeFwhm * 0.5) < 1.5 * angleTolerance);
    }

    SECTION("Plume overhead - large offset - returns correct plume properties.")
    {
        const double acutalPlumeCenter = 0.0;
        const double actualPlumeFwhm = 20.0;
        const double columnOffsetInMeasurement = -200.0;
        auto plume = GenerateGaussianScanResult(200.0, acutalPlumeCenter, actualPlumeFwhm, columnOffsetInMeasurement);
        const double expectedPlumeWidth = CalculateGaussianOneOverEFromFwhm(actualPlumeFwhm);

        // Act
        auto plumeProperties = CalculatePlumeProperties(plume, StandardMolecule::SO2, message);

        // Assert, the calculated values should be as expected
        REQUIRE(nullptr != plumeProperties);
        REQUIRE(plumeProperties->plumeCenter.HasValue());
        REQUIRE(std::abs(plumeProperties->plumeCenter.Value()) < angleTolerance);
        REQUIRE(std::abs(plumeProperties->plumeCenter2.Value()) < angleTolerance); //plumecenter2 is set to zero, all the scan2 angles are zero
        REQUIRE(std::abs(plumeProperties->plumeEdgeHigh.Value() - acutalPlumeCenter - expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties->plumeEdgeLow.Value() - acutalPlumeCenter + expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties->plumeHalfHigh.Value() - acutalPlumeCenter - actualPlumeFwhm * 0.5) < 1.5 * angleTolerance);
        REQUIRE(std::abs(plumeProperties->plumeHalfLow.Value() - acutalPlumeCenter + actualPlumeFwhm * 0.5) < 1.5 * angleTolerance);
    }

    SECTION("Plume at -45degrees - returns correct plume properties.")
    {
        const double plumeCenter = -45.0;
        auto plume = GenerateGaussianScanResult(200.0, plumeCenter, 20.0);
        const double expectedPlumeWidth = CalculateGaussianOneOverEFromFwhm(20.0);

        // Act
        auto plumeProperties = CalculatePlumeProperties(plume, StandardMolecule::SO2, message);

        // Assert, the calculated values should be as expected
        REQUIRE(nullptr != plumeProperties);
        REQUIRE(plumeProperties->plumeCenter.HasValue());
        REQUIRE(std::abs(plumeProperties->plumeCenter.Value() - plumeCenter) < 4.0);
        REQUIRE(std::abs(plumeProperties->plumeCenter2.Value()) < 1.0);
        REQUIRE(std::abs(plumeProperties->plumeEdgeHigh.Value() - plumeCenter - expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties->plumeEdgeLow.Value() - plumeCenter + expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties->plumeHalfHigh.Value() - plumeCenter - 20.0 * 0.5) < 1.5 * angleTolerance);
        REQUIRE(std::abs(plumeProperties->plumeHalfLow.Value() - plumeCenter + 20.0 * 0.5) < 1.5 * angleTolerance);
    }

    SECTION("Plume at -45degrees - large offset  - correct plume properties.")
    {
        const double acutalPlumeCenter = -45.0;
        const double actualPlumeFwhm = 20.0;
        const double columnOffsetInMeasurement = -400.0;
        auto plume = GenerateGaussianScanResult(200.0, acutalPlumeCenter, actualPlumeFwhm, columnOffsetInMeasurement);
        const double expectedPlumeWidth = CalculateGaussianOneOverEFromFwhm(20.0);

        // Act
        auto plumeProperties = CalculatePlumeProperties(plume, StandardMolecule::SO2, message);

        // Assert, the calculated values should be as expected
        REQUIRE(nullptr != plumeProperties);
        REQUIRE(plumeProperties->plumeCenter.HasValue());
        REQUIRE(std::abs(plumeProperties->plumeCenter.Value() - acutalPlumeCenter) < 4.0);
        REQUIRE(std::abs(plumeProperties->plumeCenter2.Value()) < 1.0);
        REQUIRE(std::abs(plumeProperties->plumeEdgeHigh.Value() - acutalPlumeCenter - expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties->plumeEdgeLow.Value() - acutalPlumeCenter + expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties->plumeHalfHigh.Value() - acutalPlumeCenter - 20.0 * 0.5) < 1.5 * angleTolerance);
        REQUIRE(std::abs(plumeProperties->plumeHalfLow.Value() - acutalPlumeCenter + 20.0 * 0.5) < 1.5 * angleTolerance);
    }

    SECTION("Plume at +45degrees - returns correct plume properties.")
    {
        const double plumeCenter = 45.0;
        auto plume = GenerateGaussianScanResult(200.0, plumeCenter, 20.0);
        const double expectedPlumeWidth = CalculateGaussianOneOverEFromFwhm(20.0);

        // Act
        auto plumeProperties = CalculatePlumeProperties(plume, StandardMolecule::SO2, message);

        // Assert, the calculated values should be as expected
        REQUIRE(nullptr != plumeProperties);
        REQUIRE(plumeProperties->plumeCenter.HasValue());
        REQUIRE(std::abs(plumeProperties->plumeCenter.Value() - plumeCenter) < 4.0);
        REQUIRE(std::abs(plumeProperties->plumeCenter2.Value()) < 1.0);
        REQUIRE(std::abs(plumeProperties->plumeEdgeHigh.Value() - plumeCenter - expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties->plumeEdgeLow.Value() - plumeCenter + expectedPlumeWidth * 0.5) < angleTolerance);
        REQUIRE(std::abs(plumeProperties->plumeHalfHigh.Value() - plumeCenter - 20.0 * 0.5) < 1.5 * angleTolerance);
        REQUIRE(std::abs(plumeProperties->plumeHalfLow.Value() - plumeCenter + 20.0 * 0.5) < 1.5 * angleTolerance);
    }

    SECTION("Plume at -75degrees - returns correct plume properties.")
    {
        const double acutalPlumeCenter = -75.0;
        const double actualPlumeFwhm = 20.0;
        const double columnOffsetInMeasurement = 0.0;
        auto plume = GenerateGaussianScanResult(200.0, acutalPlumeCenter, actualPlumeFwhm, columnOffsetInMeasurement);

        // Act
        auto plumeProperties = CalculatePlumeProperties(plume, StandardMolecule::SO2, message);

        // Assert, the calculated values should be as expected
        REQUIRE(nullptr != plumeProperties);
        REQUIRE(plumeProperties->plumeCenter.HasValue());
        REQUIRE(std::abs(plumeProperties->plumeCenter.Value() + 75.0) < 4.0);
        REQUIRE(std::abs(plumeProperties->plumeCenter2.Value()) < 1.0);
    }

    SECTION("Plume at +75degrees - returns correct plume properties.")
    {
        auto plume = GenerateGaussianScanResult(200.0, +75.0, 20.0);

        // Act
        auto plumeProperties = CalculatePlumeProperties(plume, StandardMolecule::SO2, message);

        // Assert, the calculated values should be as expected
        REQUIRE(nullptr != plumeProperties);
        REQUIRE(plumeProperties->plumeCenter.HasValue());
        REQUIRE(std::abs(plumeProperties->plumeCenter.Value() - 75.0) < 4.0);
        REQUIRE(std::abs(plumeProperties->plumeCenter2.Value()) < 1.0);
    }

    SECTION("Plume at -45degrees - bad values edges of scan - returns true and sets plume properties within scan range.")
    {
        const double acutalPlumeCenter = -45.0;
        const double actualPlumeFwhm = 90.0;
        const double columnOffsetInMeasurement = -100.0;
        auto plume = GenerateGaussianScanResult(200.0, acutalPlumeCenter, actualPlumeFwhm, columnOffsetInMeasurement);
        // add a couple of 'bad' measurements in the beginning...
        plume.m_spec[0].m_referenceResult[0].m_column = -5000.0;
        plume.m_spec[0].m_referenceResult[1].m_column = -5000.0;

        // Act
        auto plumeProperties = CalculatePlumeProperties(plume, StandardMolecule::SO2, message);

        // Assert, the calculated values should be as expected
        REQUIRE(nullptr != plumeProperties);
        REQUIRE(plumeProperties->plumeCenter.HasValue());
        REQUIRE(plumeProperties->plumeCenter.Value() < 90.0);
        REQUIRE(plumeProperties->plumeCenter.Value() > -90.0);
    }
}

// Test labelled as integration test since we do need to read data from file.
TEST_CASE("CalculatePlumeProperties with measured plume (BroSo2 ratio measurement)", "[PlumeProperties][IntegrationTest]")
{
    std::string message;
    novac::CScanEvaluationLogFileHandler evaluationFileHandler;
    const bool evaluationFileIsOk = evaluationFileHandler.ReadEvaluationLog(TestData::GetBrORatioEvaluationFile1());
    REQUIRE(evaluationFileIsOk); // check assumption on the setup
    REQUIRE(evaluationFileHandler.m_scan.size() == 1); // check assumption on the setup

    // Act
    auto plumeProperties = CalculatePlumeProperties(evaluationFileHandler.m_scan[0], StandardMolecule::SO2, message);

    SECTION("Correct plume center.")
    {
        REQUIRE(plumeProperties->plumeCenter.Value() == Approx(-23.12).margin(0.01));
    }

    SECTION("Correct plume edges.")
    {
        // Notice that this scan has some dark spectra in the beginning, hence the plume is not entirely complete.
        REQUIRE(plumeProperties->plumeEdgeLow.Value() == Approx(-50.0).margin(0.01));
        REQUIRE(plumeProperties->plumeEdgeHigh.Value() == Approx(0.0).margin(0.01));

        REQUIRE(plumeProperties->plumeHalfLow.Value() == Approx(-50.0).margin(0.01));
        REQUIRE(plumeProperties->plumeHalfHigh.Value() == Approx(0.0).margin(0.01));
    }

    SECTION("Sets completeness and offset.")
    {
        REQUIRE(plumeProperties->completeness.Value() == Approx(0.7).margin(0.01));
        REQUIRE(plumeProperties->offset.Value() == Approx(-1.215e18));
    }
}

TEST_CASE("CalculatePlumeCompleteness with Gaussian plume", "[PlumeProperties]")
{
    CPlumeInScanProperty plumeProperties;

    SECTION("No plume - returns false.")
    {
        // Arrange
        auto plume = GenerateGaussianPlume(0.0, 1.0, 20.0);

        // Act
        REQUIRE(false == CalculatePlumeCompleteness(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, 0.0, (long)plume.scanAngles.size(), plumeProperties));

        // Assert
        REQUIRE(plumeProperties.completeness.HasValue() == false);
    }

    SECTION("Plume overhead - returns true and sets completeness to 1.0.")
    {
        // Arrange
        const double plumeCenter = 0.0;
        auto plume = GenerateGaussianPlume(200.0, plumeCenter, 20.0);

        // Act
        REQUIRE(true == CalculatePlumeCompleteness(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, 0.0, (long)plume.scanAngles.size(), plumeProperties));

        // Assert
        REQUIRE(plumeProperties.completeness.Value() == Approx(1.0));
    }

    SECTION("Narrow plume at -45degrees - returns true and sets completeness to 1.0.")
    {
        // Arrange
        const double plumeCenter = -45.0;
        auto plume = GenerateGaussianPlume(200.0, plumeCenter, 20.0);

        // Act
        REQUIRE(true == CalculatePlumeCompleteness(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, 0.0, (long)plume.scanAngles.size(), plumeProperties));

        // Assert, the plume is relatively narrow and does not extend down to horizon. Hence everything is visible.
        REQUIRE(plumeProperties.completeness.Value() == Approx(1.0).margin(0.1));
    }

    SECTION("Narrow plume at -45degree - large negative offset - returns true and sets completeness to 1.0.")
    {
        // Arrange
        const double plumeCenter = -45.0;
        const double plumeAmplitude = 200.0;
        const double columnOffset = -400.0;
        auto plume = GenerateGaussianPlume(plumeAmplitude, plumeCenter, 20.0, columnOffset);

        // Act
        REQUIRE(true == CalculatePlumeCompleteness(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, 0.0, (long)plume.scanAngles.size(), plumeProperties));

        // Assert, the plume is relatively narrow and does not extend down to horizon. Hence everything is visible.
        REQUIRE(plumeProperties.completeness.Value() == Approx(1.0).margin(0.1));
    }

    SECTION("Wide plume at -45degrees - returns true and sets completeness to 1.0.")
    {
        // Arrange
        const double plumeCenter = -45.0;
        auto plume = GenerateGaussianPlume(200.0, plumeCenter, 40.0);

        // Act
        REQUIRE(true == CalculatePlumeCompleteness(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, 0.0, (long)plume.scanAngles.size(), plumeProperties));

        // Assert, the plume is relatively wide but still does not extend down to horizon. Hence everything is visible.
        REQUIRE(plumeProperties.completeness.Value() == Approx(1.0).margin(0.1));
    }

    SECTION("Plume at +45degrees - returns correct completeness.")
    {
        // Arrange
        const double plumeCenter = 45.0;
        auto plume = GenerateGaussianPlume(200.0, plumeCenter, 20.0);

        // Act
        REQUIRE(true == CalculatePlumeCompleteness(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, 0.0, (long)plume.scanAngles.size(), plumeProperties));

        // Assert, the plume is relatively narrow and does not extend down to horizon. Hence everything is visible.
        REQUIRE(plumeProperties.completeness.Value() == Approx(1.0).margin(0.1));
    }

    SECTION("Wide plume at +45degrees - returns correct completeness.")
    {
        // Arrange
        const double plumeCenter = 45.0;
        const double plumeFwhm = 90.0;
        auto plume = GenerateGaussianPlume(200.0, plumeCenter, plumeFwhm);

        // Act
        REQUIRE(true == CalculatePlumeCompleteness(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, 0.0, (long)plume.scanAngles.size(), plumeProperties));

        // Assert, the plume is relatively wide and does extend down to horizon.
        REQUIRE(plumeProperties.completeness.Value() == Approx(0.7).margin(0.1));
    }


    SECTION("Plume at -75degrees - returns correct completeness.")
    {
        // Arrange
        auto plume = GenerateGaussianPlume(200.0, -75.0, 20.0);

        // Act
        REQUIRE(true == CalculatePlumeCompleteness(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, 0.0, (long)plume.scanAngles.size(), plumeProperties));

        // Assert, the plume is relatively narrow but also at a low angle and does hence extend down to horizon.
        REQUIRE(plumeProperties.completeness.Value() == Approx(0.7).margin(0.1));
    }

    SECTION("Plume at +75degrees - returns correct completeness.")
    {
        // Arrange
        auto plume = GenerateGaussianPlume(200.0, +75.0, 20.0);

        // Act
        REQUIRE(true == CalculatePlumeCompleteness(plume.scanAngles, plume.phi, plume.columns, plume.columnErrors, plume.badEvaluation, 0.0, (long)plume.scanAngles.size(), plumeProperties));

        // Assert, the plume is relatively narrow but also at a low angle and does hence extend down to horizon.
        REQUIRE(plumeProperties.completeness.Value() == Approx(0.7).margin(0.1));
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
    plumeProperties.offset = CalculatePlumeOffset(evaluationFileHandler.m_scan[0], 0); // The offset needs to be filled in first, as a preparation

    // Act
    REQUIRE(true == CalculatePlumeCompleteness(evaluationFileHandler.m_scan[0], 0, plumeProperties));

    SECTION("Sets correct completeness.")
    {
        REQUIRE(plumeProperties.completeness.Value() == Approx(0.7).margin(0.01));
    }

    SECTION("Sets correct plume center.")
    {
        REQUIRE(plumeProperties.plumeCenter.Value() == Approx(-23.12).margin(0.01));
    }

    SECTION("Sets correct plume edges.")
    {
        // Notice that this scan has some dark spectra in the beginning, hence the plume is not entirely complete.
        REQUIRE(plumeProperties.plumeEdgeLow.Value() == Approx(-50.0).margin(0.01));
        REQUIRE(plumeProperties.plumeEdgeHigh.Value() == Approx(0.0).margin(0.01));

        REQUIRE(plumeProperties.plumeHalfLow.Value() == Approx(-50.0).margin(0.01));
        REQUIRE(plumeProperties.plumeHalfHigh.Value() == Approx(0.0).margin(0.01));
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

    // Act
    const auto calculatedOffset = CalculatePlumeOffset(evaluationFileHandler.m_scan[0], 0);

    // Assert
    REQUIRE(calculatedOffset.HasValue());
    REQUIRE(calculatedOffset.Value() == Approx(-1.2e18).margin(1e17));
}