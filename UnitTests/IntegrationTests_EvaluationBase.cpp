#include <SpectralEvaluation/File/SpectrumIO.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Evaluation/EvaluationBase.h>
#include "catch.hpp"
#include "TestData.h"

namespace novac
{

//Region Helper methods
CFitWindow PrepareFitWindow()
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

        int retCode = window.ref[refIdx].ReadCrossSectionDataFromFile();
        REQUIRE(retCode == 0);

        ++refIdx;
    }

    return window;
}

CSpectrum ReadSpectrumNumber(const std::string& scanFile, int number)
{
    CSpectrumIO reader;

    CSpectrum spectrum;
    bool success = reader.ReadSpectrum(scanFile, number, spectrum);
    REQUIRE(success);
    spectrum.Div(spectrum.NumSpectra());

    return spectrum;
}

CSpectrum ReadSkySpectrum(const std::string& scanFile)
{
    return ReadSpectrumNumber(scanFile, 0);
}

CSpectrum ReadDarkSpectrum(const std::string& scanFile)
{
    return ReadSpectrumNumber(scanFile, 1);
}

//endregion


TEST_CASE("Evaluate Avaspec spectrum number eight in scan", "[Evaluate][EvaluationBase]")
{
    const auto scanFile = TestData::GetMeasuredSpectrumName_2009175M1();

    CSpectrumIO reader;

    CFitWindow window = PrepareFitWindow();

    // Read the spectra
    CEvaluationBase sut;
    sut.SetFitWindow(window);

    CSpectrum skySpectrum = ReadSkySpectrum(scanFile);
    CSpectrum darkSpectrum = ReadDarkSpectrum(scanFile);
    CSpectrum spectrumToEvaluate = ReadSpectrumNumber(scanFile, 8);
    REQUIRE(spectrumToEvaluate.ScanAngle() == Approx(-68.0)); // Verification that we did indeed read the correct spectrum

    // Prepare the spectra
    skySpectrum.Sub(darkSpectrum);
    spectrumToEvaluate.Sub(darkSpectrum);

    sut.SetSkySpectrum(skySpectrum);

    // Act
    int returnCode = sut.Evaluate(spectrumToEvaluate);

    // Assert
    REQUIRE(returnCode == 0);

    REQUIRE(window.nRef == int(sut.m_result.m_referenceResult.size()));

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

TEST_CASE("Evaluate Avaspec spectrum number 21 in scan", "[Evaluate][EvaluationBase]")
{
    const auto scanFile = TestData::GetMeasuredSpectrumName_2009175M1();

    CSpectrumIO reader;

    CFitWindow window = PrepareFitWindow();

    // Read the spectra
    CEvaluationBase sut;
    sut.SetFitWindow(window);

    CSpectrum skySpectrum = ReadSkySpectrum(scanFile);
    CSpectrum darkSpectrum = ReadDarkSpectrum(scanFile);
    CSpectrum spectrumToEvaluate = ReadSpectrumNumber(scanFile, 21);
    REQUIRE(spectrumToEvaluate.ScanAngle() == Approx(-21.0)); // Verification that we did indeed read the correct spectrum

    // Prepare the spectra
    skySpectrum.Sub(darkSpectrum);
    spectrumToEvaluate.Sub(darkSpectrum);

    sut.SetSkySpectrum(skySpectrum);

    // Act
    int returnCode = sut.Evaluate(spectrumToEvaluate);

    // Assert
    REQUIRE(returnCode == 0);

    REQUIRE(window.nRef == int(sut.m_result.m_referenceResult.size()));

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


TEST_CASE("EvaluateShift Avaspec spectrum number 28 in scan", "[Evaluate][EvaluationBase]")
{
    const auto scanFile = TestData::GetMeasuredSpectrumName_2009175M1();

    CSpectrumIO reader;

    CFitWindow window = PrepareFitWindow();
    window.fraunhoferRef.m_path = TestData::GetSyntheticFraunhoferSpectrumName_2009175M1();
    int retCode = window.fraunhoferRef.ReadCrossSectionDataFromFile();
    REQUIRE(retCode == 0);

    // Read the spectra
    CEvaluationBase sut;
    sut.SetFitWindow(window);

    CSpectrum skySpectrum = ReadSkySpectrum(scanFile);
    CSpectrum darkSpectrum = ReadDarkSpectrum(scanFile);
    CSpectrum spectrumToEvaluate = ReadSpectrumNumber(scanFile, 28);
    REQUIRE(spectrumToEvaluate.ScanAngle() == Approx(3.0)); // Verification that we did indeed read the correct spectrum

    // Prepare the spectra
    skySpectrum.Sub(darkSpectrum);
    spectrumToEvaluate.Sub(darkSpectrum);

    sut.SetSkySpectrum(skySpectrum);

    double resultingShift = 0.0;
    double resultingShiftError = 0.0;
    double resultingSqueeze = 0.0;
    double resultingSqueezeError = 0.0;

    // Act
    int returnCode = sut.EvaluateShift(spectrumToEvaluate, resultingShift, resultingShiftError, resultingSqueeze, resultingSqueezeError);

    // Assert
    REQUIRE(returnCode == 0);

    REQUIRE(resultingShift == Approx(0.219).margin(0.01));
    REQUIRE(resultingShiftError == Approx(0.091).margin(0.01));
    REQUIRE(resultingSqueeze == Approx(1.0));
    REQUIRE(resultingSqueezeError == Approx(0.0));
}

}