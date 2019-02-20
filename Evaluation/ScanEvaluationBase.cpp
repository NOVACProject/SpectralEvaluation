#include "ScanEvaluationBase.h"
#include "../Spectra/Spectrum.h"
#include "../File/ScanFileHandler.h"
#include "../File/STDFile.h"
#include "../File/TXTFile.h"

namespace Evaluation
{
    ScanEvaluationBase::ScanEvaluationBase()
    {

    }

    ScanEvaluationBase::~ScanEvaluationBase()
    {

    }

    bool ScanEvaluationBase::ReadSpectrumFromFile(const std::string& fullFilename, CSpectrum& spec)
    {
        if (fullFilename.size() < 3)
        {
            return false;
        }
        if (SpectrumIO::CSTDFile::ReadSpectrum(spec, fullFilename))
        {
            return true;
        }
        if (SpectrumIO::CTXTFile::ReadSpectrum(spec, fullFilename))
        {
            return true;
        }
        return false;
    }

    bool ScanEvaluationBase::GetOffsetSpectrum(FileHandler::CScanFileHandler &scan, const Configuration::CDarkSettings& darkSettings, CSpectrum &offsetSpectrum)
    {
        if (darkSettings.m_offsetOption == Configuration::DARK_MODEL_OPTION::USER_SUPPLIED)
        {
            return ReadSpectrumFromFile(darkSettings.m_offsetSpec, offsetSpectrum);
        }
        else
        {
            return 0 != scan.GetOffset(offsetSpectrum);
        }
    }

    bool ScanEvaluationBase::GetDarkCurrentSpectrum(FileHandler::CScanFileHandler &scan, const Configuration::CDarkSettings& darkSettings, CSpectrum &darkCurrent, bool& needsOffsetCorrection)
    {
        needsOffsetCorrection = true;
        if (darkSettings.m_darkCurrentOption == Configuration::DARK_MODEL_OPTION::USER_SUPPLIED)
        {
            if (ReadSpectrumFromFile(darkSettings.m_darkCurrentSpec, darkCurrent))
            {
                needsOffsetCorrection = false;
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            return 0 != scan.GetDarkCurrent(darkCurrent);
        }
    }

    bool ScanEvaluationBase::ModelDarkSpectrum(FileHandler::CScanFileHandler &scan, const Configuration::CDarkSettings& darkSettings, const CSpectrum &spec, CSpectrum &dark)
    {
        bool offsetCorrectDC = true;

        CSpectrum offset;
        if (!GetOffsetSpectrum(scan, darkSettings, offset))
        {
            m_lastErrorMessage = "Could not retrieve an offset spectrum from the scan.";
            return false;
        }

        CSpectrum darkCurrent;
        if (!GetDarkCurrentSpectrum(scan, darkSettings, darkCurrent, offsetCorrectDC))
        {
            m_lastErrorMessage = "Could not retrieve a dark-current spectrum from the scan.";
            return false;
        }

        if (offset.m_length == darkCurrent.m_length && offset.m_length > 0)
        {
            // 3c-1 Scale the offset spectrum to the measured
            offset.Mult(spec.NumSpectra() / (double)offset.NumSpectra());
            offset.m_info.m_numSpec = spec.NumSpectra();

            // 3c-2 Remove offset from the dark-current spectrum
            if (offsetCorrectDC)
            {
                CSpectrum offset_dc = offset;
                offset_dc.Mult(darkCurrent.NumSpectra() / (double)offset_dc.NumSpectra());
                darkCurrent.Sub(offset_dc);
            }

            // 3c-3 Scale the dark-current spectrum to the measured
            darkCurrent.Mult((spec.NumSpectra() * spec.ExposureTime()) / (double)(darkCurrent.NumSpectra() * darkCurrent.ExposureTime()));
            darkCurrent.m_info.m_numSpec = spec.NumSpectra();

            // 3d. Make the dark-spectrum
            dark.Clear();
            dark.m_length = offset.m_length;
            dark.m_info.m_interlaceStep = offset.m_info.m_interlaceStep;
            dark.m_info.m_channel = offset.m_info.m_channel;
            dark.Add(offset);
            dark.Add(darkCurrent);

            // If the dark-spectrum is read out in an interlaced way then interpolate it back to it's original state
            if (dark.m_info.m_interlaceStep > 1)
            {
                dark.InterpolateSpectrum();
            }

            return true;
        }
        else
        {
            m_lastErrorMessage = "Invalid offset/dark-current spectra the lengths does not agree.";
            return false;
        }
    }

    bool ScanEvaluationBase::GetDark(FileHandler::CScanFileHandler& scan, const CSpectrum &spec, CSpectrum &dark, const Configuration::CDarkSettings* darkSettings)
    {
        CSpectrum offset, darkCurrent, offset_dc;

        // 1. The user wants to take the dark spectrum directly from the measurement
        //		as the second spectrum in the scan.
        if (darkSettings == nullptr ||
            darkSettings->m_darkSpecOption == Configuration::DARK_SPEC_OPTION::MEASURE ||
            darkSettings->m_darkSpecOption == Configuration::DARK_SPEC_OPTION::MODEL_SOMETIMES)
        {
            if (0 != scan.GetDark(dark))
            {
                m_lastErrorMessage = "Could not read dark-spectrum from scan " + scan.GetFileName();
                return false;
            }

            // if there is no dark spectrum but one offset and one dark-current spectrum,
            //	then read those instead and model the dark spectrum
            if (dark.m_length == 0)
            {
                if (ModelDarkSpectrum(scan, *darkSettings, spec, dark))
                {
                    m_lastErrorMessage = "Warning: Incorrect settings: check settings for dark current correction";
                    return false;
                }
                else
                {
                    m_lastErrorMessage = "WARNING: NO DARK SPECTRUM FOUND IN SCAN. INCORRECT DARK CURRENT CORRECTION";
                    return false;
                }
            }

            // If the dark-spectrum is read out in an interlaced way then interpolate it back to it's original state
            if (dark.m_info.m_interlaceStep > 1) {
                dark.InterpolateSpectrum();
            }

            // Check so that the exposure-time of the dark-spectrum is same as the 
            //	exposure time of the measured spectrum
            if (dark.ExposureTime() != spec.ExposureTime()) {
                m_lastErrorMessage = "WARNING: EXPOSURE-TIME OF DARK-SPECTRUM IS NOT SAME AS FOR MEASURED SPECTRUM. INCORRECT DARK-CORRECTION!! " + scan.GetFileName();
            }

            // Make sure that there are the same number of exposures in the
            //	dark-spectrum as in the measured spectrum
            if (dark.NumSpectra() != spec.NumSpectra()) {
                dark.Mult(spec.NumSpectra() / (double)dark.NumSpectra());
            }

            return true;
        }

        // 3. The user wants to model the dark spectrum
        if (darkSettings->m_darkSpecOption == Configuration::DARK_SPEC_OPTION::MODEL_ALWAYS)
        {
            if (ModelDarkSpectrum(scan, *darkSettings, spec, dark))
                return true;
            else
                return false;
        }

        // 4. The user has his own favourite dark-spectrum that he wants to use
        if (darkSettings->m_darkSpecOption == Configuration::DARK_SPEC_OPTION::DARK_USER_SUPPLIED)
        {
            if (!ReadSpectrumFromFile(darkSettings->m_offsetSpec, dark))
                return false;

            // If the dark-spectrum is read out in an interlaced way then interpolate it back to it's original state
            if (dark.m_info.m_interlaceStep > 1) {
                dark.InterpolateSpectrum();
            }

            return true;
        }

        // something is not implemented
        return false;
    }

}