#include <SpectralEvaluation/Evaluation/ScanEvaluationBase.h>
#include <SpectralEvaluation/Evaluation/DarkSpectrum.h>
#include <SpectralEvaluation/Evaluation/FitWindow.h>
#include <SpectralEvaluation/Evaluation/EvaluationBase.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Spectra/SpectrometerModel.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/File/SpectrumIO.h>
#include <SpectralEvaluation/File/STDFile.h>
#include <SpectralEvaluation/File/TXTFile.h>
#include <SpectralEvaluation/StringUtils.h>
#include <SpectralEvaluation/Configuration/DarkSettings.h>
#include <SpectralEvaluation/Configuration/SkySettings.h>

#include <sstream>

namespace novac
{
bool ReadSpectrumFromFile(const std::string& fullFilename, CSpectrum& spec)
{
    if (fullFilename.size() < 3)
    {
        return false;
    }
    if (CSTDFile::ReadSpectrum(spec, fullFilename))
    {
        return true;
    }
    if (CTXTFile::ReadSpectrum(spec, fullFilename))
    {
        return true;
    }
    return false;
}

ScanEvaluationBase::ScanEvaluationBase(novac::ILogger& log)
    : m_log(log)
{
}

ScanEvaluationBase::~ScanEvaluationBase()
{
}

bool ScanEvaluationBase::GetDark(CScanFileHandler& scan, const CSpectrum& spec, CSpectrum& dark, const Configuration::CDarkSettings* darkSettings)
{
    // Using DarkSpectrum.h
    if (darkSettings == nullptr)
    {
        Configuration::CDarkSettings defaultDarkSettings;
        return ::novac::GetDark(scan, spec, defaultDarkSettings, dark, m_lastErrorMessage);
    }
    else
    {
        return ::novac::GetDark(scan, spec, *darkSettings, dark, m_lastErrorMessage);
    }
}

bool ScanEvaluationBase::GetSky(CScanFileHandler& scan, const Configuration::CSkySettings& settings, CSpectrum& sky)
{
    novac::LogContext context; // TODO: Get from input

    // If the sky spectrum is the first spectrum in the scan
    if (settings.skyOption == Configuration::SKY_OPTION::MEASURED_IN_SCAN)
    {
        if (0 != scan.GetSky(sky))
        {
            m_lastErrorMessage = "Could not get sky spectrum from scan, no spectrum is labelled as sky.";
            return false;
        }

        if (sky.m_info.m_interlaceStep > 1)
        {
            sky.InterpolateSpectrum();
        }

        return true;
    }

    // If the sky spectrum is the average of all credible spectra
    if (settings.skyOption == Configuration::SKY_OPTION::AVERAGE_OF_GOOD_SPECTRA_IN_SCAN)
    {
        if (!GetSkySpectrumFromAverageOfGoodSpectra(scan, sky))
        {
            m_lastErrorMessage = "Failed to get sky spectrum from scan, no spectra are good enough.";
        }
        return true;
    }

    // If the user wants to use another spectrum than 'sky' as reference-spectrum...
    if (settings.skyOption == Configuration::SKY_OPTION::SPECTRUM_INDEX_IN_SCAN)
    {
        if (0 == scan.GetSpectrum(context, sky, settings.indexInScan))
        {
            m_lastErrorMessage = "Failed to get sky spectrum from index in scan, could not find given spectrum.";
            return false;
        }

        if (sky.m_info.m_interlaceStep > 1)
        {
            sky.InterpolateSpectrum();
        }

        return true;
    }

    // If the user has supplied a special sky-spectrum to use
    if (settings.skyOption == Configuration::SKY_OPTION::USER_SUPPLIED)
    {
        return GetSkySpectrumFromFile(settings.skySpectrumFile, sky);
    }

    return false;
}

bool ScanEvaluationBase::GetSkySpectrumFromAverageOfGoodSpectra(CScanFileHandler& scan, CSpectrum& sky) const
{
    const int interlaceSteps = scan.GetInterlaceSteps();
    const int startChannel = scan.GetStartChannel();
    const int fitLow = m_fitLow / interlaceSteps - startChannel;
    const int fitHigh = m_fitHigh / interlaceSteps - startChannel;
    int nofSpectraAveraged = 0;

    novac::LogContext context; // TODO: Get from input

    CSpectrum tmp;
    scan.GetSky(tmp);
    scan.ResetCounter();

    // Get the maximum intensity of this spectrometer model (with a little bit of margin)
    const double spectrometerDynamicRange = (CSpectrometerDatabase::GetInstance().GetModel(tmp.m_info.m_specModelName).maximumIntensityForSingleReadout - 20);

    double spectrumIntensityInFitRegion = tmp.MaxValue(fitLow, fitHigh);
    if (spectrumIntensityInFitRegion < spectrometerDynamicRange * tmp.NumSpectra() && !tmp.IsDark())
    {
        sky = tmp;
        nofSpectraAveraged = 1;
    }
    else
    {
        sky.Clear();
    }

    while (scan.GetNextSpectrum(context, tmp))
    {
        spectrumIntensityInFitRegion = tmp.MaxValue(fitLow, fitHigh);
        if (spectrumIntensityInFitRegion < spectrometerDynamicRange * tmp.NumSpectra() && !tmp.IsDark())
        {
            sky.Add(tmp);
            ++nofSpectraAveraged;
        }
    }
    scan.ResetCounter();

    if (sky.m_info.m_interlaceStep > 1)
    {
        sky.InterpolateSpectrum();
    }

    return (nofSpectraAveraged > 0);
}

bool ScanEvaluationBase::GetSkySpectrumFromFile(const std::string& filename, CSpectrum& sky)
{
    if (filename.size() < 5)
    {
        m_lastErrorMessage = "Cannot get the sky-spectrum, the filename was empty";
        return false;
    }
    const std::string extension = Right(filename, 4);

    if (EqualsIgnoringCase(extension, ".pak"))
    {
        CSpectrumIO reader;
        return reader.ReadSpectrum(filename, 0, sky);
    }

    if (EqualsIgnoringCase(extension, ".std"))
    {
        return CSTDFile::ReadSpectrum(sky, filename);
    }

    // If we don't recognize the sky-spectrum format
    m_lastErrorMessage = "Unknown format for sky spectrum. Please use .pak or .std";
    return false;
}

double GetSpectrumMaximumIntensity(const CSpectrum& spectrum)
{
    if (EqualsIgnoringCase(spectrum.m_info.m_specModelName, "S2000"))
    {
        // Default spectrometer model. May or may not be correct
        const auto& model = CSpectrometerDatabase::GetInstance().GuessModelFromSerial(spectrum.m_info.m_device);
        return model.maximumIntensityForSingleReadout;
    }
    else
    {
        // The model is not the default, assume this is correct.
        return CSpectrometerDatabase::GetInstance().GetModel(spectrum.m_info.m_specModelName).maximumIntensityForSingleReadout;
    }
}

int ScanEvaluationBase::GetIndexOfSpectrumWithBestIntensity(novac::LogContext context, const CFitWindow& fitWindow, const novac::SpectrometerModel& spectrometerModel, CScanFileHandler& scan)
{
    const int INDEX_OF_SKYSPECTRUM = -1;
    const int NO_SPECTRUM_INDEX = -2;
    CSpectrum sky;

    // Find the spectrum for which we should determine shift & squeeze
    //      This spectrum should have high enough intensity in the fit-region
    //      without being saturated.
    int indexOfMostSuitableSpectrum = NO_SPECTRUM_INDEX;
    scan.GetSky(sky);
    double bestSaturation = -1.0;
    const double skyFitIntensity = sky.MaxValue(fitWindow.fitLow, fitWindow.fitHigh);

    const double skyFitSaturation = (sky.NumSpectra() > 0) ?
        (skyFitIntensity / (spectrometerModel.maximumIntensityForSingleReadout * sky.NumSpectra())) :
        skyFitIntensity / spectrometerModel.maximumIntensityForSingleReadout;

    if (skyFitSaturation < 0.9 && skyFitSaturation > 0.1)
    {
        indexOfMostSuitableSpectrum = INDEX_OF_SKYSPECTRUM;
        bestSaturation = skyFitSaturation;
    }

    scan.ResetCounter(); // start from the beginning

    CSpectrum spectrum;
    int curIndex = 0;
    while (scan.GetSpectrum(context, spectrum, curIndex))
    {
        const double fitIntensity = spectrum.MaxValue(fitWindow.fitLow, fitWindow.fitHigh);
        const double numSpectra = std::max(1.0, (double)spectrum.NumSpectra());
        const double fitSaturation = fitIntensity / (numSpectra * spectrometerModel.maximumIntensityForSingleReadout);

        // Check if this spectrum is good...
        if (fitSaturation < 0.9 && fitSaturation > 0.1 && fitSaturation > bestSaturation)
        {
            indexOfMostSuitableSpectrum = curIndex;
            bestSaturation = fitSaturation;
        }

        // Go to the next spectrum
        ++curIndex;
    }

    return indexOfMostSuitableSpectrum;
}

CEvaluationBase* ScanEvaluationBase::FindOptimumShiftAndSqueezeFromFraunhoferReference(
    novac::LogContext context,
    const CFitWindow& fitWindow,
    const novac::SpectrometerModel& spectrometerModel,
    const Configuration::CDarkSettings& darkSettings,
    CScanFileHandler& scan)
{
    // Check that the Fraunhofer reference has been read in
    if (fitWindow.fraunhoferRef.m_data == nullptr || fitWindow.fraunhoferRef.m_data->m_crossSection.size() == 0)
    {
        m_log.Information(context, "Cannot determine shift and squeeze from Fraunhofer reference. Reference has no values.");
        return nullptr;
    }

    const int INDEX_OF_SKYSPECTRUM = -1;
    const int NO_SPECTRUM_INDEX = -2;

    // 1. Find the spectrum for which we should determine shift & squeeze
    //      This spectrum should have high enough intensity in the fit-region
    //      without being saturated.
    const int indexOfMostSuitableSpectrum = GetIndexOfSpectrumWithBestIntensity(context, fitWindow, spectrometerModel, scan);

    // 2. Get the spectrum we should evaluate...
    CSpectrum spectrum;
    if (indexOfMostSuitableSpectrum == NO_SPECTRUM_INDEX)
    {
        m_log.Information(context, "Could not find any suitable spectrum to determine shift from.");
        return nullptr; // we could not find any good spectrum to use...
    }
    else if (indexOfMostSuitableSpectrum == INDEX_OF_SKYSPECTRUM)
    {
        scan.GetSky(spectrum);
        m_log.Information(context, "Determining shift and squeeze from given Fraunhofer Reference spectrum and sky spectrum");
    }
    else
    {
        scan.GetSpectrum(context, spectrum, indexOfMostSuitableSpectrum);
        std::stringstream msg;
        msg << "Determining shift and squeeze from given Fraunhofer Reference spectrum and spectrum " << indexOfMostSuitableSpectrum;
        m_log.Information(context, msg.str());
    }

    if (spectrum.NumSpectra() > 0 && !m_averagedSpectra)
    {
        spectrum.Div(spectrum.NumSpectra());
    }

    // verification here, the measured spectrum should now have an intensity range which matches the range of the spectrometer model.
    if (spectrum.MaxValue(0, spectrum.m_length) > spectrometerModel.maximumIntensityForSingleReadout + 1)
    {
        std::stringstream msg;
        msg << "Unexpected spectrum range. Max value: " << spectrum.MaxValue(0, spectrum.m_length);
        msg << " is larger than expected maximum intensity of: " << spectrometerModel.maximumIntensityForSingleReadout << " for a device of type " << spectrometerModel.modelName;
        m_log.Error(context, msg.str());
    }

    CSpectrum dark;
    if (!GetDark(scan, spectrum, dark, &darkSettings))
    {
        m_log.Information(context, "Failed to get dark spectrum, determination of shift-and-squeeze from Fraunhofer lines failed.");
        return nullptr; // fail
    }
    if (dark.NumSpectra() > 0 && !m_averagedSpectra)
    {
        dark.Div(dark.NumSpectra());
    }
    spectrum.Sub(dark);

    // 3. Do the evaluation.
    CFitWindow copyOfFitWindow = fitWindow;
    CEvaluationBase shiftEvaluator(copyOfFitWindow, m_log);

    novac::ShiftEvaluationResult shiftResult;
    if (shiftEvaluator.EvaluateShift(context, spectrum, shiftResult))
    {
        // We failed to make the fit, what shall we do now??
        std::stringstream msg;
        msg << "Failed to determine shift and squeeze from Fraunhofer reference in scan " << scan.GetFileName() << ". Will proceed with default shift and squeeze";
        m_lastErrorMessage = msg.str();
        m_log.Information(context, msg.str());
        return nullptr;
    }

    if (std::abs(shiftResult.shiftError) < 1 &&
        std::abs(shiftResult.squeezeError) < 0.01 &&
        std::abs(shiftResult.chi2) < 2.0)
    {
        CFitWindow improvedFitWindow = fitWindow;

        // The fit is good enough to use the values
        for (size_t it = 0; it < improvedFitWindow.reference.size(); ++it)
        {
            improvedFitWindow.reference[it].m_shiftOption = SHIFT_TYPE::SHIFT_FIX;
            improvedFitWindow.reference[it].m_squeezeOption = SHIFT_TYPE::SHIFT_FIX;
            improvedFitWindow.reference[it].m_shiftValue = shiftResult.shift;
            improvedFitWindow.reference[it].m_squeezeValue = shiftResult.squeeze;
        }

        std::stringstream msg;
        msg << "Determined shift: " << shiftResult.shift << " +- " << shiftResult.shiftError << "; Squeeze: " << shiftResult.squeeze << " +- " << shiftResult.squeezeError << ". Doas chi2: " << shiftResult.chi2;
        m_log.Information(context, msg.str());
        return new CEvaluationBase(improvedFitWindow, m_log);
    }
    else
    {
        std::stringstream msg;
        msg << "Inadequate result from deciding shift. Output was, shift: " << shiftResult.shift << " +- " << shiftResult.shiftError << "; Squeeze: " << shiftResult.squeeze << " +- " << shiftResult.squeezeError << ". Doas chi2: " << shiftResult.chi2;
        m_log.Information(context, msg.str());
    }

    m_log.Information(context, "Fit not good enough. Will proceed with default shift and squeeze.");
    return nullptr;
}

}