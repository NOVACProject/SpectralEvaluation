#include <SpectralEvaluation/Evaluation/DoasFitPreparation.h>
#include <SpectralEvaluation/Evaluation/BasicMath.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/VectorUtils.h>

using namespace novac;

std::vector<double> DoasFitPreparation::PrepareSkySpectrum(const CSpectrum& skySpectrum, FIT_TYPE doasFitType)
{
    if (doasFitType == FIT_TYPE::FIT_HP_DIV)
    {
        throw std::invalid_argument("Cannot prepare the sky spectrum for a HP_DIV type of doas fit, the sky spectrum should not be included for this type of fit.");
    }

    CBasicMath math;

    const int spectrumLength = static_cast<int>(skySpectrum.m_length);

    std::vector<double> filteredSkySpectrum{ skySpectrum.m_data, skySpectrum.m_data + spectrumLength };

    // Remove the offset. TODO: Make the range an input parameter
    RemoveOffset(filteredSkySpectrum, 50, 200);

    if (doasFitType == FIT_TYPE::FIT_HP_SUB)
    {
        math.HighPassBinomial(filteredSkySpectrum.data(), spectrumLength, 500);
    }

    math.Log(filteredSkySpectrum.data(), spectrumLength);

    return filteredSkySpectrum;
}

std::vector<double> DoasFitPreparation::PrepareMeasuredSpectrum(const CSpectrum& measuredSpectrum, const CSpectrum& skySpectrum, FIT_TYPE doasFitType)
{
    if (measuredSpectrum.m_length != skySpectrum.m_length)
    {
        throw std::invalid_argument("Cannot prepare the measured spectrum for a DOAS fit if the measured and the sky spectra does not have equal length.");
    }

    CBasicMath math;

    const int spectrumLength = static_cast<int>(measuredSpectrum.m_length);

    std::vector<double> filteredMeasSpectrum{ measuredSpectrum.m_data, measuredSpectrum.m_data + spectrumLength };

    // Always start by removing the offset. TODO: Make the range an input parameter
    RemoveOffset(filteredMeasSpectrum, 50, 200);

    if (doasFitType == FIT_TYPE::FIT_HP_DIV)
    {
        // Divide the measured spectrum with the sky spectrum
        math.Div(filteredMeasSpectrum.data(), skySpectrum.m_data, spectrumLength, 0.0);

        // high pass filter
        math.HighPassBinomial(filteredMeasSpectrum.data(), spectrumLength, 500);
    }
    else if (doasFitType == FIT_TYPE::FIT_HP_SUB)
    {
        // high pass filter
        math.HighPassBinomial(filteredMeasSpectrum.data(), spectrumLength, 500);
    }
    else if (doasFitType == FIT_TYPE::FIT_POLY)
    {
        // ... nothing more needs to be done here ...
    }
    else
    {
        throw std::invalid_argument("Unknown type of Doas fit passed to DoasFitPreparation::PrepareMeasuredSpectrum");
    }

    // Always end with taking the log of the spectrum, such that we end up in Optical Density space.
    math.Log(filteredMeasSpectrum.data(), spectrumLength);

    return filteredMeasSpectrum;
}

void DoasFitPreparation::RemoveOffset(std::vector<double>& spectrum, int startIndex, int endIndex)
{
    if (startIndex == endIndex)
    {
        return;
    }

    const double skySpectrumOffset = Average(begin(spectrum) + startIndex, begin(spectrum) + endIndex);

    CBasicMath math;
    math.Sub(spectrum.data(), static_cast<int>(spectrum.size()), skySpectrumOffset);
}