#pragma once

// ------------ This file contains a basic helper class for organizing our test data ------------

namespace novac
{
class TestData
{
public:
    static std::string GetTestDataDirectory()
    {
#ifdef _MSC_VER
        return std::string("../TestData/");
#else
        return std::string("TestData/");
#endif // _MSC_VER 
    }

    // region Instrument Calibration

    static std::string GetMeasuredSpectrumName_I2J8549()
    {
        return GetTestDataDirectory() + std::string("I2J8549/I2J8549_170216_1230_0.pak");
    }

    static std::string GetInitialPixelToWavelengthCalibration_I2J8549()
    {
        return GetTestDataDirectory() + std::string("I2J8549/I2J8549.clb");
    }

    static std::string GetInitialInstrumentLineShapefile_I2J8549()
    {
        return GetTestDataDirectory() + std::string("I2J8549/I2J8549_302nm.slf");
    }

    static std::string GetMeasuredSpectrumName_2009175M1()
    {
        return GetTestDataDirectory() + std::string("2009175M1/2009175M1_211214_1817_0.pak");
    }

    static std::string GetSyntheticFraunhoferSpectrumName_2009175M1()
    {
        return GetTestDataDirectory() + std::string("2009175M1/2002127M1_Fraunhofer.txt");
    }

    static std::string GetInitialPixelToWavelengthCalibration_2009175M1()
    {
        return GetTestDataDirectory() + std::string("2009175M1/DD2J3040_MASTER_SO2_HP500_PPMM.txt");
    }

    static std::vector<std::string> GetReferences_2009175M1()
    {
        return std::vector <std::string>{
            GetTestDataDirectory() + std::string("2009175M1/2002127M1_SO2_Bogumil_293K.txt"),
                GetTestDataDirectory() + std::string("2009175M1/2002127M1_Ring_HP500_PPMM.txt"),
                GetTestDataDirectory() + std::string("2009175M1/2002127M1_O3_Voigt_223K.txt")
        };
    }

    static std::string GetInitialPixelToWavelengthCalibration_D2J2200()
    {
        return GetTestDataDirectory() + std::string("D2J2200/D2J2200_Master.clb");
    }

    static std::string GetMeasuredMercurySpectrum_D2J2200()
    {
        return GetTestDataDirectory() + std::string("D2J2200/Hg_D2J2200_all.Master.Sample.txt");
    }

    static std::string GetSkySpectrumName_MAYP11440()
    {
        return GetTestDataDirectory() + std::string("MAYP11440/sky_0.STD");
    }

    static std::string GetDarkSpectrumName_MAYP11440()
    {
        return GetTestDataDirectory() + std::string("MAYP11440/dark_0.STD");
    }

    static std::string GetInitialPixelToWavelengthCalibration_MAYP11440()
    {
        return GetTestDataDirectory() + std::string("MAYP11440/MAYP11440_SO2_293K_Bogumil_334nm.txt");
    }

    static std::string GetMeasuredSpectrumName_FLMS14634()
    {
        return GetTestDataDirectory() + std::string("FLMS14634/00007_0.STD");
    }

    static std::string GetDarkSpectrumName_FLMS14634()
    {
        return GetTestDataDirectory() + std::string("FLMS14634/dark_0.STD");
    }

    static std::string GetInitialPixelToWavelengthCalibration_FLMS14634()
    {
        return GetTestDataDirectory() + std::string("FLMS14634/FLMS14634.clb");
    }

    static std::string GetInitialInstrumentLineShapefile_FLMS14634()
    {
        return GetTestDataDirectory() + std::string("FLMS14634/FLMS14634_302nm.slf");
    }

    static std::string GetSyntheticFraunhoferSpectrumName_FLMS14634()
    {
        return GetTestDataDirectory() + std::string("FLMS14634/FLMS14634_Fraunhofer.txt");
    }

    static std::string GetMeasuredSpectrumName_I2J0093()
    {
        return GetTestDataDirectory() + std::string("I2P0093/00000_0.STD");
    }

    static std::string GetDarkSpectrumName_I2P0093()
    {
        return GetTestDataDirectory() + std::string("I2P0093/dark_0.STD");
    }

    static std::string GetInitialPixelToWavelengthCalibration_I2P0093()
    {
        return GetTestDataDirectory() + std::string("I2P0093/I2P0093_Master.clb");
    }

    static std::string GetInitialInstrumentLineShape_I2P0093()
    {
        return GetTestDataDirectory() + std::string("I2P0093/I2P0093_302nm_Master.slf");
    }

    static std::string GetSolarAtlasFile()
    {
        return GetTestDataDirectory() + std::string("SOLARFL_296-440nm.xs");
    }

    static std::string GetMercurySpectrumWithoutWavelengthCalibration()
    {
        return GetTestDataDirectory() + std::string("MercurySpectra/hglampnov152021.std");
    }

    static std::string GetMercurySpectrumWithoutWavelengthCalibration_DarkSpectrumFile()
    {
        return GetTestDataDirectory() + std::string("MercurySpectra/hglampnov152021_dark.std");
    }

    static std::string GetWavelengthCalibratedStdFileName()
    {
        return GetTestDataDirectory() + std::string("WavelengthCalibratedStdFile.std");
    }

    static std::string GetHighResolutionSO2CrossSectionFile()
    {
        return GetTestDataDirectory() + std::string("SO2_Bogumil(2003)_293K_239-395nm.txt");
    }

    static std::string GetQDoasConvolvedSO2CrossSectionFile()
    {
        return GetTestDataDirectory() + std::string("SO2_QDOAS.xs");
    }

    static std::string GetTemporaryInstrumentCalibrationStdFileName()
    {
        return GetTestDataDirectory() + std::string("Temporary_InstrumentCalibration.std");
    }

    static std::string GetTemporaryFitWindowFileName()
    {
        return GetTestDataDirectory() + std::string("Temporary_FitWindowFile.nfw");
    }

    static std::string GetTemporaryConfigurationFileName()
    {
        return GetTestDataDirectory() + std::string("Temporary_Congfiguration.config");
    }

    // endregion

    // region Evaluation log file formats

    static std::string GetEvaluationLogfile1()
    {
        return GetTestDataDirectory() + std::string("EvaluationLogs/") + std::string("EvaluationLog_1904156M1_2020.10.07.txt");
    }

    static std::string GetEvaluationLogfile2()
    {
        return GetTestDataDirectory() + std::string("EvaluationLogs/") + std::string("ReEvaluationLog_SO2.txt");
    }

    // endregion

    //  region BrO to SO2 ratio calculations

    static std::string GetBrORatioScanFile1()
    {
        return GetTestDataDirectory() + std::string("BrORatio/") + std::string("D2J2124_160331_1608_0.pak");
    }

    static std::string GetBrORatioEvaluationFile1()
    {
        return GetTestDataDirectory() + std::string("BrORatio/") + std::string("D2J2124_160331_1608_0.txt");
    }

    static std::string GetBrORatioFitWindowFileSO2()
    {
        return GetTestDataDirectory() + std::string("BrORatio/") + std::string("SO2-314-326nm.nfw");
    }

    static std::string GetBrORatioFitWindowFileBrO()
    {
        return GetTestDataDirectory() + std::string("BrORatio/") + std::string("BrO-330-352nm.nfw");
    }

    // Notice that this is a file where we don't expect to get a ratio due to too wide plume
    static std::string GetBrORatioScanFile2()
    {
        return GetTestDataDirectory() + std::string("BrORatio/") + std::string("D2J2124_160331_2049_0.pak");
    }

    // Notice that this is a file where we don't expect to get a ratio due to too wide plume
    static std::string GetBrORatioEvaluationFile2()
    {
        return GetTestDataDirectory() + std::string("BrORatio/") + std::string("D2J2124_160331_2049_0.txt");
    }

    static std::string GetBrORatioScanFile3()
    {
        return GetTestDataDirectory() + std::string("BrORatio/") + std::string("D2J2124_160331_1510_0.pak");
    }

    static std::string GetBrORatioEvaluationFile3()
    {
        return GetTestDataDirectory() + std::string("BrORatio/") + std::string("D2J2124_160331_1510_0.txt");
    }

    // endregion

};
}
