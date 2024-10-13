#include <SpectralEvaluation/Evaluation/DarkSpectrum.h> // the file to test
#include <SpectralEvaluation/File/ScanFileHandler.h> // The implementation of the IScanSpectrumSource
#include <SpectralEvaluation/Configuration/DarkSettings.h>
#include "catch.hpp"
#include "TestData.h"

namespace novac
{
TEST_CASE("GetDark Can retrieve the the measured dark spectrum", "[DarkSpectrum][GetDark]")
{
    // Prepare the test by reading in a .pak-file which we know contains a dark spectrum.
    novac::ConsoleLog log;
    novac::LogContext context;
    CScanFileHandler fileHandler(log);
    const bool scanFileIsOk = fileHandler.CheckScanFile(context, TestData::GetBrORatioScanFile1());
    REQUIRE(scanFileIsOk); // check assumption on the setup
    CSpectrum measuredSpectrum;
    const bool scanFileContainsSky = 0 == fileHandler.GetSky(measuredSpectrum);
    REQUIRE(scanFileContainsSky);
    std::string errorMessage;
    CSpectrum resultingSpectrum;
    Configuration::CDarkSettings settings;
    settings.m_darkSpecOption = Configuration::DARK_SPEC_OPTION::MEASURED_IN_SCAN;

    // Act
    const bool result = GetDark(fileHandler, measuredSpectrum, settings, resultingSpectrum, errorMessage);

    // Assert
    REQUIRE(result == true);
    REQUIRE(errorMessage.size() == 0);
    REQUIRE(resultingSpectrum.m_length == 2048);
    REQUIRE(resultingSpectrum.m_info.m_name == "dark");
    REQUIRE(resultingSpectrum.m_info.m_numSpec == 15);
}
}