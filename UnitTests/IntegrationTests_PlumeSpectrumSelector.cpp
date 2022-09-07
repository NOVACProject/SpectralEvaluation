#include <iostream>

#include <SpectralEvaluation/Configuration/RatioEvaluationSettings.h>
#include <SpectralEvaluation/Evaluation/PlumeSpectrumSelector.h>
#include <SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/File/ScanEvaluationLogFileHandler.h>
#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>
#include "catch.hpp"
#include "TestData.h"

namespace novac
{
    TEST_CASE("PlumeSpectrumSelector returns nullptr if input not readable", "[PlumeSpectrumSelector][IntegrationTest][Ratios]")
    {
        novac::CScanFileHandler fileHandler;
        novac::BasicScanEvaluationResult evaluationResult;
        novac::CPlumeInScanProperty plumeInScanProperties;
        Configuration::RatioEvaluationSettings settings;
        PlumeSpectrumSelector sut;

        SECTION("Invalid scan file")
        {
            const auto result = sut.CreatePlumeSpectra(fileHandler, evaluationResult, plumeInScanProperties, settings);

            REQUIRE(nullptr == result);
        }

        SECTION("Invalid evaluation result")
        {
            REQUIRE(true == fileHandler.CheckScanFile(TestData::GetBrORatioScanFile1()));

            const auto result = sut.CreatePlumeSpectra(fileHandler, evaluationResult, plumeInScanProperties, settings);

            REQUIRE(nullptr == result);
        }
    }

    // happy case where there is a good plume and we should be able to extract in plume and out of plume spectra
    TEST_CASE("PlumeSpectrumSelector returns expected in and out of plume for good scan - scan file 1", "[PlumeSpectrumSelector][IntegrationTest][Ratios]")
    {
        novac::CScanFileHandler fileHandler;
        novac::CScanEvaluationLogFileHandler evaluationFileHandler;
        novac::BasicScanEvaluationResult evaluationResult;
        Configuration::RatioEvaluationSettings settings;
        PlumeSpectrumSelector sut;

        // Prepare the test by reading in the .pak-file and the evaluation result and calculate the plume-properties from the result.
        const bool scanFileIsOk = fileHandler.CheckScanFile(TestData::GetBrORatioScanFile1());
        REQUIRE(scanFileIsOk); // check assumption on the setup

        const bool evaluationFileIsOk = evaluationFileHandler.ReadEvaluationLog(TestData::GetBrORatioEvaluationFile1());
        REQUIRE(evaluationFileIsOk); // check assumption on the setup
        REQUIRE(evaluationFileHandler.m_scan.size() == 1); // check assumption on the setup

        novac::CPlumeInScanProperty plumeInScanProperties;
        novac::CalculatePlumeOffset(evaluationFileHandler.m_scan[0], 0, plumeInScanProperties);
        REQUIRE(true == novac::CalculatePlumeCompleteness(evaluationFileHandler.m_scan[0], 0, plumeInScanProperties));
        const int so2Index = 1; // here, O3 is actually the first specie

        // Act
        const auto result = sut.CreatePlumeSpectra(fileHandler, evaluationFileHandler.m_scan[0], plumeInScanProperties, settings, so2Index);

        // Assert
        REQUIRE(result != nullptr);
        REQUIRE(result->inPlumeSpectrum != nullptr);
        REQUIRE(result->referenceSpectrum != nullptr);
        REQUIRE(result->darkSpectrum != nullptr);

        // There should be 10 spectra selected for in-plume (each having 15 readouts)
        REQUIRE(result->inPlumeSpectrum->m_info.m_numSpec == 150);
        REQUIRE(result->inPlumeSpectrumIndices.size() == 10);
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 15));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 16));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 17));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 18));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 19));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 20));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 21));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 22));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 23));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 24));

        // There should be 10 spectra selected for out-of-plume (each having 15 readouts)
        REQUIRE(result->referenceSpectrum->m_info.m_numSpec == 150);
        REQUIRE(result->referenceSpectrumIndices.size() == 10);
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 12));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 13));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 14));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 25));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 26));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 27));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 28));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 29));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 30));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 32));
    }

    TEST_CASE("PlumeSpectrumSelector returns expected in and out of plume for good scan - scan file 2", "[PlumeSpectrumSelector][IntegrationTest][Ratios]")
    {
        novac::CScanFileHandler fileHandler;
        novac::CScanEvaluationLogFileHandler evaluationFileHandler;
        novac::BasicScanEvaluationResult evaluationResult;
        Configuration::RatioEvaluationSettings settings;
        PlumeSpectrumSelector sut;

        // Prepare the test by reading in the .pak-file and the evaluation result and calculate the plume-properties from the result.
        const bool scanFileIsOk = fileHandler.CheckScanFile(TestData::GetBrORatioScanFile2());
        REQUIRE(scanFileIsOk); // check assumption on the setup

        const bool evaluationFileIsOk = evaluationFileHandler.ReadEvaluationLog(TestData::GetBrORatioEvaluationFile2());
        REQUIRE(evaluationFileIsOk); // check assumption on the setup
        REQUIRE(evaluationFileHandler.m_scan.size() == 1); // check assumption on the setup

        novac::CPlumeInScanProperty plumeInScanProperties;
        novac::CalculatePlumeOffset(evaluationFileHandler.m_scan[0], 0, plumeInScanProperties);
        REQUIRE(true == novac::CalculatePlumeCompleteness(evaluationFileHandler.m_scan[0], 1, plumeInScanProperties));
        const int so2Index = 1; // here, O3 is actually the first specie

        // Act
        std::string errorMessage;
        const auto result = sut.CreatePlumeSpectra(fileHandler, evaluationFileHandler.m_scan[0], plumeInScanProperties, settings, so2Index, &errorMessage);

        // Assert
        REQUIRE(result != nullptr);
        REQUIRE(result->inPlumeSpectrum != nullptr);
        REQUIRE(result->referenceSpectrum != nullptr);
        REQUIRE(result->darkSpectrum != nullptr);

        // There should be 10 spectra selected for in-plume (each having 15 readouts)
        REQUIRE(result->inPlumeSpectrum->m_info.m_numSpec == 150);
        REQUIRE(result->inPlumeSpectrumIndices.size() == 10);
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 21));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 22));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 23));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 24));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 25));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 26));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 27));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 28));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 29));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 30));

        // There should be 10 spectra selected for out-of-plume (each having 15 readouts)
        REQUIRE(result->referenceSpectrum->m_info.m_numSpec == 150);
        REQUIRE(result->referenceSpectrumIndices.size() == 10);
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 41));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 42));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 43));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 44));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 45));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 46));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 47));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 48));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 49));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 50));
    }

    // another happy case where there is a good plume and we should be able to extract in plume and out of plume spectra (but trickier than the one above)
    TEST_CASE("PlumeSpectrumSelector returns expected in and out of plume for good scan - scan file 3", "[PlumeSpectrumSelector][IntegrationTest][Ratios]")
    {
        novac::CScanFileHandler fileHandler;
        novac::CScanEvaluationLogFileHandler evaluationFileHandler;
        novac::BasicScanEvaluationResult evaluationResult;
        Configuration::RatioEvaluationSettings settings;
        PlumeSpectrumSelector sut;

        // Prepare the test by reading in the .pak-file and the evaluation result and calculate the plume-properties from the result.
        const bool scanFileIsOk = fileHandler.CheckScanFile(TestData::GetBrORatioScanFile3());
        REQUIRE(scanFileIsOk); // check assumption on the setup

        const bool evaluationFileIsOk = evaluationFileHandler.ReadEvaluationLog(TestData::GetBrORatioEvaluationFile3());
        REQUIRE(evaluationFileIsOk); // check assumption on the setup
        REQUIRE(evaluationFileHandler.m_scan.size() == 1); // check assumption on the setup

        novac::CPlumeInScanProperty plumeInScanProperties;
        novac::CalculatePlumeOffset(evaluationFileHandler.m_scan[0], 0, plumeInScanProperties);
        REQUIRE(true == novac::CalculatePlumeCompleteness(evaluationFileHandler.m_scan[0], 1, plumeInScanProperties));
        const int so2Index = 1; // here, O3 is actually the first specie

        // Act
        std::string errorMessage;
        const auto result = sut.CreatePlumeSpectra(fileHandler, evaluationFileHandler.m_scan[0], plumeInScanProperties, settings, so2Index, &errorMessage);

        // Assert
        REQUIRE(result != nullptr);
        REQUIRE(result->inPlumeSpectrum != nullptr);
        REQUIRE(result->referenceSpectrum != nullptr);
        REQUIRE(result->darkSpectrum != nullptr);

        // There should be 10 spectra selected for in-plume (each having 15 readouts)
        REQUIRE(result->inPlumeSpectrum->m_info.m_numSpec == 150);
        REQUIRE(result->inPlumeSpectrumIndices.size() == 10);
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 14));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 15));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 16));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 17));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 18));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 19));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 20));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 21));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 22));
        REQUIRE(result->inPlumeSpectrumIndices.end() != std::find(begin(result->inPlumeSpectrumIndices), end(result->inPlumeSpectrumIndices), 23));

        // There should be 10 spectra selected for out-of-plume (each having 15 readouts)
        REQUIRE(result->referenceSpectrum->m_info.m_numSpec == 150);
        REQUIRE(result->referenceSpectrumIndices.size() == 10);
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 11));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 12));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 13));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 24));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 25));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 26));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 27));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 28));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 31));
        REQUIRE(result->referenceSpectrumIndices.end() != std::find(begin(result->referenceSpectrumIndices), end(result->referenceSpectrumIndices), 32));
    }
}
