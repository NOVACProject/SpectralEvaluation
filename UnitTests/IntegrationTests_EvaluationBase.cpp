#include <SpectralEvaluation/File/SpectrumIO.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Evaluation/EvaluationBase.h>
#include "catch.hpp"
#include "TestData.h"

namespace novac
{

//Region Helper methods
static CFitWindow PrepareFitWindow()
{
    const auto references = TestData::GetReferences_2009175M1();
    REQUIRE(3 == references.size()); // Assumption here

    CFitWindow window;
    window.fitLow = 475;
    window.fitHigh = 643;
    window.fitType = novac::FIT_TYPE::FIT_HP_DIV;
    window.nRef = (int)references.size();
    int refIdx = 0;
    for (auto& reference : references)
    {
        window.ref[refIdx].m_path = reference;
        window.ref[refIdx].m_columnOption = novac::SHIFT_TYPE::SHIFT_FREE;
        window.ref[refIdx].m_shiftOption = novac::SHIFT_TYPE::SHIFT_FIX;
        window.ref[refIdx].m_shiftValue = 0.0;
        window.ref[refIdx].m_squeezeOption = novac::SHIFT_TYPE::SHIFT_FIX;
        window.ref[refIdx].m_squeezeValue = 1.0;

        window.ref[refIdx].ReadCrossSectionDataFromFile();

        ++refIdx;
    }

    return window;
}

static CSpectrum ReadSpectrumNumber(const std::string& scanFile, int number)
{
    CSpectrumIO reader;

    CSpectrum spectrum;
    bool success = reader.ReadSpectrum(scanFile, number, spectrum);
    REQUIRE(success);
    spectrum.Div(spectrum.NumSpectra());

    return spectrum;
}

static CSpectrum ReadSkySpectrum(const std::string& scanFile)
{
    return ReadSpectrumNumber(scanFile, 0);
}

static CSpectrum ReadDarkSpectrum(const std::string& scanFile)
{
    return ReadSpectrumNumber(scanFile, 1);
}

//endregion


TEST_CASE("Evaluate Avaspec spectrum number eight in scan", "[Evaluate][EvaluationBase][2009175M1_211214_1817_0]")
{
    const auto scanFile = TestData::GetMeasuredSpectrumName_2009175M1();

    novac::ConsoleLog log;
    novac::LogContext context;

    CSpectrumIO reader;
    CFitWindow window = PrepareFitWindow();

    CEvaluationBase sut(log);
    sut.SetFitWindow(window);

    // Read the spectra (and divide them by the number of readouts already)
    CSpectrum skySpectrum = ReadSkySpectrum(scanFile);
    const CSpectrum darkSpectrum = ReadDarkSpectrum(scanFile);
    CSpectrum spectrumToEvaluate = ReadSpectrumNumber(scanFile, 8);
    REQUIRE(spectrumToEvaluate.ScanAngle() == Approx(-68.0)); // Verification that we did indeed read the correct spectrum

    // Prepare the spectra
    skySpectrum.Sub(darkSpectrum);
    spectrumToEvaluate.Sub(darkSpectrum);

    sut.SetSkySpectrum(skySpectrum);

    // Act
    const int returnCode = sut.Evaluate(spectrumToEvaluate);

    // Assert
    REQUIRE(returnCode == 0);

    REQUIRE(window.nRef == sut.m_result.m_referenceResult.size());

    REQUIRE(sut.m_result.m_delta == Approx(0.0752).margin(0.001));
    REQUIRE(sut.m_result.m_chiSquare == Approx(0.0119).margin(0.001));

    // Verify that the results are as expected
    REQUIRE(sut.m_result.m_referenceResult[0].m_column == Approx(10.5).margin(0.1));
    REQUIRE(sut.m_result.m_referenceResult[0].m_columnError == Approx(12.8).margin(0.1));
    REQUIRE(sut.m_result.m_referenceResult[0].m_shift == 0.0);
    REQUIRE(sut.m_result.m_referenceResult[0].m_squeeze == 1.0);

    REQUIRE(sut.m_result.m_referenceResult[1].m_column == Approx(-0.0029).margin(0.01));
    REQUIRE(sut.m_result.m_referenceResult[1].m_columnError == Approx(0.00376).margin(0.01));
    REQUIRE(sut.m_result.m_referenceResult[1].m_shift == 0.0);
    REQUIRE(sut.m_result.m_referenceResult[1].m_squeeze == 1.0);

    REQUIRE(sut.m_result.m_referenceResult[2].m_column == Approx(1050.0).margin(10.0));
    REQUIRE(sut.m_result.m_referenceResult[2].m_columnError == Approx(97.0).margin(1.0));
    REQUIRE(sut.m_result.m_referenceResult[2].m_shift == 0.0);
    REQUIRE(sut.m_result.m_referenceResult[2].m_squeeze == 1.0);
}

TEST_CASE("Evaluate Avaspec spectrum number 21 in scan", "[Evaluate][EvaluationBase][2009175M1_211214_1817_0]")
{
    const auto scanFile = TestData::GetMeasuredSpectrumName_2009175M1();

    novac::ConsoleLog log;
    novac::LogContext context;

    CSpectrumIO reader;
    CFitWindow window = PrepareFitWindow();

    CEvaluationBase sut(log);
    sut.SetFitWindow(window);

    // Read the spectra (and divide them by the number of readouts already)
    CSpectrum skySpectrum = ReadSkySpectrum(scanFile);
    const CSpectrum darkSpectrum = ReadDarkSpectrum(scanFile);
    CSpectrum spectrumToEvaluate = ReadSpectrumNumber(scanFile, 21);
    REQUIRE(spectrumToEvaluate.ScanAngle() == Approx(-21.0)); // Verification that we did indeed read the correct spectrum

    // Prepare the spectra
    skySpectrum.Sub(darkSpectrum);
    spectrumToEvaluate.Sub(darkSpectrum);

    sut.SetSkySpectrum(skySpectrum);

    // Act
    const int returnCode = sut.Evaluate(spectrumToEvaluate);

    // Assert
    REQUIRE(returnCode == 0);

    REQUIRE(window.nRef == sut.m_result.m_referenceResult.size());

    REQUIRE(sut.m_result.m_delta == Approx(0.0429).margin(0.001));
    REQUIRE(sut.m_result.m_chiSquare == Approx(0.0066).margin(0.001));

    // Verify that the results are as expected
    REQUIRE(sut.m_result.m_referenceResult[0].m_column == Approx(-15.5).margin(0.1));
    REQUIRE(sut.m_result.m_referenceResult[0].m_columnError == Approx(9.54).margin(0.1));
    REQUIRE(sut.m_result.m_referenceResult[0].m_shift == 0.0);
    REQUIRE(sut.m_result.m_referenceResult[0].m_squeeze == 1.0);

    REQUIRE(sut.m_result.m_referenceResult[1].m_column == Approx(-0.00359).margin(0.001));
    REQUIRE(sut.m_result.m_referenceResult[1].m_columnError == Approx(0.00280).margin(0.001));
    REQUIRE(sut.m_result.m_referenceResult[1].m_shift == 0.0);
    REQUIRE(sut.m_result.m_referenceResult[1].m_squeeze == 1.0);

    REQUIRE(sut.m_result.m_referenceResult[2].m_column == Approx(-169.0).margin(1.0));
    REQUIRE(sut.m_result.m_referenceResult[2].m_columnError == Approx(72.4).margin(1.0));
    REQUIRE(sut.m_result.m_referenceResult[2].m_shift == 0.0);
    REQUIRE(sut.m_result.m_referenceResult[2].m_squeeze == 1.0);
}

TEST_CASE("EvaluateShift Avaspec spectrum number 28 in scan", "[Evaluate][EvaluationBase][2009175M1_211214_1817_0]")
{
    const auto scanFile = TestData::GetMeasuredSpectrumName_2009175M1();

    CSpectrumIO reader;

    novac::ConsoleLog log;
    novac::LogContext context;

    CFitWindow window = PrepareFitWindow();
    window.UV = 0; // Avaspec and the UV option are not great together
    window.fraunhoferRef.m_path = TestData::GetSyntheticFraunhoferSpectrumName_2009175M1();
    window.fraunhoferRef.ReadCrossSectionDataFromFile();
    novac::HighPassFilter_Ring(*window.fraunhoferRef.m_data); // filter the fraunhofer reference, to match the other references.

    CEvaluationBase sut(log);
    sut.SetFitWindow(window);

    // Read the spectra (and divide them by the number of readouts already)
    CSpectrum skySpectrum = ReadSkySpectrum(scanFile);
    const CSpectrum darkSpectrum = ReadDarkSpectrum(scanFile);
    CSpectrum spectrumToEvaluate = ReadSpectrumNumber(scanFile, 28);
    REQUIRE(spectrumToEvaluate.ScanAngle() == Approx(3.0)); // Verification that we did indeed read the correct spectrum

    // Prepare the spectra
    skySpectrum.Sub(darkSpectrum);
    spectrumToEvaluate.Sub(darkSpectrum);

    novac::ShiftEvaluationResult result;

    // Act
    const int returnCode = sut.EvaluateShift(context, spectrumToEvaluate, result);

    // Assert
    REQUIRE(returnCode == 0);

    REQUIRE(result.shift == Approx(-1.64).margin(0.01));
    REQUIRE(result.shiftError == Approx(0.061).margin(0.01));
    REQUIRE(result.squeeze == Approx(1.0));
    REQUIRE(result.squeezeError == Approx(0.0));
    REQUIRE(result.chi2 < 0.24);
}

}