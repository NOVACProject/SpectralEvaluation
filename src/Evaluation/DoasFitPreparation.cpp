#include <SpectralEvaluation/Evaluation/DoasFitPreparation.h>
#include <SpectralEvaluation/Evaluation/BasicMath.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Spectra/Scattering.h>
#include <SpectralEvaluation/VectorUtils.h>
#include <limits>

using namespace novac;

std::vector<double> DoasFitPreparation::PrepareSkySpectrum(const CSpectrum& skySpectrum, FIT_TYPE doasFitType)
{
    const IndexRange defaultRange{ 50, 200 };
    return PrepareSkySpectrum(skySpectrum, doasFitType, defaultRange);
}

std::vector<double> DoasFitPreparation::PrepareSkySpectrum(const CSpectrum& skySpectrum, FIT_TYPE doasFitType, const IndexRange& offsetRemovalRange)
{
    if (doasFitType == FIT_TYPE::FIT_HP_DIV)
    {
        throw std::invalid_argument("Cannot prepare the sky spectrum for a HP_DIV type of doas fit, the sky spectrum should not be included for this type of fit.");
    }

    CBasicMath math;

    const int spectrumLength = static_cast<int>(skySpectrum.m_length);

    std::vector<double> filteredSkySpectrum{ skySpectrum.m_data, skySpectrum.m_data + spectrumLength };

    // Remove the offset.
    RemoveOffset(
        filteredSkySpectrum,
        static_cast<int>(offsetRemovalRange.from),
        static_cast<int>(offsetRemovalRange.to));

    if (doasFitType == FIT_TYPE::FIT_HP_SUB)
    {
        math.HighPassBinomial(filteredSkySpectrum.data(), spectrumLength, 500);
    }

    math.Log(filteredSkySpectrum.data(), spectrumLength);

    return filteredSkySpectrum;
}

std::vector<double> DoasFitPreparation::PrepareMeasuredSpectrum(const CSpectrum& measuredSpectrum, const CSpectrum& skySpectrum, FIT_TYPE doasFitType)
{
    const IndexRange defaultRange{ 50, 200 };
    return PrepareMeasuredSpectrum(measuredSpectrum, skySpectrum, doasFitType, defaultRange);
}

std::vector<double> DoasFitPreparation::PrepareMeasuredSpectrum(const CSpectrum& measuredSpectrum, const CSpectrum& skySpectrum, FIT_TYPE doasFitType, const IndexRange& offsetRemovalRange)
{
    if (measuredSpectrum.m_length != skySpectrum.m_length)
    {
        throw std::invalid_argument("Cannot prepare the measured spectrum for a DOAS fit if the measured and the sky spectra does not have equal length.");
    }

    CBasicMath math;

    const int spectrumLength = static_cast<int>(measuredSpectrum.m_length);

    std::vector<double> filteredMeasSpectrum{ measuredSpectrum.m_data, measuredSpectrum.m_data + spectrumLength };

    // Always start by removing the offset.
    RemoveOffset(
        filteredMeasSpectrum,
        static_cast<int>(offsetRemovalRange.from),
        static_cast<int>(offsetRemovalRange.to));

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

    const double spectrumOffset = Average(begin(spectrum) + startIndex, begin(spectrum) + endIndex);

    CBasicMath math;
    math.Sub(spectrum.data(), static_cast<int>(spectrum.size()), spectrumOffset);
}

void DoasFitPreparation::RemoveOffset(CSpectrum& spectrum, int startIndex, int endIndex)
{
    if (startIndex == endIndex)
    {
        return;
    }
    std::vector<double> values(spectrum.m_data + startIndex, spectrum.m_data + endIndex);

    const double spectrumOffset = Average(values);

    CBasicMath math;
    math.Sub(spectrum.m_data, static_cast<int>(spectrum.m_length), spectrumOffset);
}

std::vector<double> DoasFitPreparation::PrepareRingSpectrum(const CSpectrum& skySpectrum, FIT_TYPE doasFitType)
{
    CSpectrum copyOfSky = skySpectrum;
    auto ringSpectrum = Doasis::Scattering::CalcRingSpectrum(copyOfSky);

    const int spectrumLength = static_cast<int>(skySpectrum.m_length);

    if (doasFitType == FIT_TYPE::FIT_HP_DIV || doasFitType == FIT_TYPE::FIT_HP_SUB)
    {
        CBasicMath math;
        math.HighPassBinomial(ringSpectrum.m_data, spectrumLength, 500);
    }
    else if (doasFitType == FIT_TYPE::FIT_POLY)
    {
        // Do nothing here.
    }
    else
    {
        throw std::invalid_argument("Unknown type of Doas fit passed to DoasFitPreparation::PrepareRingSpectrum");
    }

    std::vector<double> filteredMeasSpectrum{ ringSpectrum.m_data, ringSpectrum.m_data + spectrumLength };
    return filteredMeasSpectrum;
}

std::vector<double> DoasFitPreparation::PrepareIntensitySpacePolynomial(const CSpectrum& skySpectrum)
{
    const size_t spectrumLength = static_cast<int>(skySpectrum.m_length);
    std::vector<double> result;
    result.resize(spectrumLength);

    for (size_t ii = 0; ii < spectrumLength; ++ii)
    {
        result[ii] = std::abs(skySpectrum.m_data[ii]) < std::numeric_limits<double>::epsilon() ? 0.0 : 1.0 / skySpectrum.m_data[ii];
    }

    return result;

}