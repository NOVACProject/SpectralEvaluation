#include <SpectralEvaluation/Evaluation/EvaluationResult.h>
#include <SpectralEvaluation/Evaluation/DoasFit.h>
#include <SpectralEvaluation/Spectra/SpectrumInfo.h>
#include <SpectralEvaluation/Spectra/SpectrometerModel.h>
#include <cstring>
#include <cmath>

namespace novac
{

CEvaluationResult::CEvaluationResult()
    : m_delta(0.0), m_chiSquare(0.0), m_stepNum(0)
{
    memset(m_polynomial, 0, 6 * sizeof(double));
    m_evaluationStatus = 0;

    m_referenceResult.reserve(4);
}

CEvaluationResult::CEvaluationResult(const CEvaluationResult& b)
    : m_chiSquare(b.m_chiSquare),
    m_delta(b.m_delta),
    m_stepNum(b.m_stepNum),
    m_evaluationStatus(b.m_evaluationStatus),
    m_referenceResult(b.m_referenceResult)
{
    memcpy(this->m_polynomial, b.m_polynomial, 6 * sizeof(double));
}

CEvaluationResult& CEvaluationResult::operator =(const CEvaluationResult& b)
{
    m_referenceResult = std::vector<CReferenceFitResult>(b.m_referenceResult);
    memcpy(this->m_polynomial, b.m_polynomial, 6 * sizeof(double));
    m_chiSquare = b.m_chiSquare;
    m_delta = b.m_delta;
    m_stepNum = b.m_stepNum;

    m_evaluationStatus = b.m_evaluationStatus;

    return *this;
}

CEvaluationResult::CEvaluationResult(const DoasResult& doasResult)
    : m_chiSquare(doasResult.chiSquare),
    m_delta(doasResult.delta),
    m_stepNum(doasResult.iterations)
{
    const int numberOfCoefficientsToCopy = std::min(6, static_cast<int>(doasResult.polynomialCoefficients.size()));
    for (int ii = 0; ii < numberOfCoefficientsToCopy; ++ii)
    {
        m_polynomial[ii] = doasResult.polynomialCoefficients[ii];
    }
    for (int ii = numberOfCoefficientsToCopy; ii < 6; ++ii)
    {
        m_polynomial[ii] = 0.0;
    }

    for (const auto& referenceResult : doasResult.referenceResult)
    {
        CReferenceFitResult copiedReference(referenceResult.column, referenceResult.columnError, referenceResult.shift, referenceResult.shiftError, referenceResult.squeeze, referenceResult.squeezeError);
        copiedReference.m_specieName = referenceResult.name;
        m_referenceResult.push_back(copiedReference);
    }

    m_evaluationStatus = 0; // no info about this in DoasResult..
}

size_t CEvaluationResult::InsertSpecie(const std::string& name)
{
    CReferenceFitResult ref;
    ref.m_specieName = std::string(name);
    m_referenceResult.push_back(ref);

    return m_referenceResult.size();
}

// TODO: this should be in a separate class.
bool CEvaluationResult::CheckGoodnessOfFit(const CSpectrumInfo& info, const SpectrometerModel* spectrometer, float chi2Limit, float upperLimit, float lowerLimit)
{
    // assume that this is an ok evaluation
    m_evaluationStatus &= ~MARK_BAD_EVALUATION;

    // The maximum intensity for one spectrum (# bits in the ADC)
    const double maxInt = (spectrometer != nullptr) ? (spectrometer->maximumIntensityForSingleReadout) : CSpectrometerDatabase::GetInstance().GetModel(info.m_specModelName).maximumIntensityForSingleReadout;

    // The maximum saturation-level in the fit-region
    double fitSaturation = 0.0;
    if (info.m_fitIntensity <= 1.0)
    {
        fitSaturation = info.m_fitIntensity;
    }
    else if (info.m_average)
    {
        fitSaturation = info.m_fitIntensity / maxInt;
    }
    else
    {
        if (info.m_numSpec > 0)
            fitSaturation = info.m_fitIntensity / (maxInt * info.m_numSpec);
        else
        {
            int numSpec = (int)floor(info.m_peakIntensity / maxInt); // a guess for the number of co-adds
            fitSaturation = info.m_fitIntensity / (maxInt * numSpec);
        }
    }

    // The offset of the spectrum
    double offset = 0.0;
    if (info.m_average)
    {
        offset = info.m_offset / maxInt;
    }
    else if (info.m_numSpec > 0)
    {
        offset = info.m_offset / (maxInt * info.m_numSpec);
    }
    else
    {
        int numSpec = (int)floor(info.m_peakIntensity / maxInt); // a guess for the number of co-adds
        offset = info.m_offset / (maxInt * numSpec);
    }

    // first check the intensity of the spectrum in the fit region
    if (upperLimit > -1)
    {
        if (fitSaturation > upperLimit)
            m_evaluationStatus |= MARK_BAD_EVALUATION;
    }
    else
    {
        if (fitSaturation > 0.99)
            m_evaluationStatus |= MARK_BAD_EVALUATION;
    }

    // first check the intensity of the spectrum in the fit region
    if (lowerLimit > -1)
    {
        if (fitSaturation < lowerLimit)
            m_evaluationStatus |= MARK_BAD_EVALUATION;
    }
    else
    {
        if ((fitSaturation - offset) < 0.025)
            m_evaluationStatus |= MARK_BAD_EVALUATION;
    }

    // then check the chi2 of the fit
    if (chi2Limit > -1)
    {
        if (m_chiSquare > chi2Limit)
            m_evaluationStatus |= MARK_BAD_EVALUATION;
    }
    else
    {
        if (m_chiSquare > 0.9)
            m_evaluationStatus |= MARK_BAD_EVALUATION;
    }

    return (m_evaluationStatus & MARK_BAD_EVALUATION);
}

bool CEvaluationResult::MarkAs(int MARK_FLAG)
{
    // check the flag
    switch (MARK_FLAG)
    {
    case MARK_BAD_EVALUATION: break;
    case MARK_DELETED: break;
    default: return false;
    }

    // set the corresponding bit
    m_evaluationStatus |= MARK_FLAG;

    return true;
}

bool CEvaluationResult::RemoveMark(int MARK_FLAG)
{
    // check the flag
    switch (MARK_FLAG)
    {
    case MARK_BAD_EVALUATION: break;
    case MARK_DELETED: break;
    default: return false;
    }

    // remove the corresponding bit
    m_evaluationStatus &= ~MARK_FLAG;

    return true;
}

}
