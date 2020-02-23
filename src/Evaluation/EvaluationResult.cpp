#include <SpectralEvaluation/Evaluation/EvaluationResult.h>
#include <SpectralEvaluation/Spectra/SpectrumInfo.h>
#include <cstring>
#include <cmath>

namespace Evaluation
{
    CEvaluationResult::CEvaluationResult()
        :m_chiSquare(0.0), m_delta(0.0), m_stepNum(0)
    {
        memset(m_polynomial, 0, 5 * sizeof(float));
        m_evaluationStatus = 0;

        m_referenceResult.reserve(4);
    }

    CEvaluationResult::CEvaluationResult(const CEvaluationResult &b)
    {
        const int nRef = (int)b.m_referenceResult.size();
        m_referenceResult.resize(nRef);

        for (long i = 0; i < nRef; ++i) {
            CReferenceFitResult ref;
            ref.m_column = b.m_referenceResult[i].m_column;
            ref.m_columnError = b.m_referenceResult[i].m_columnError;
            ref.m_shift = b.m_referenceResult[i].m_shift;
            ref.m_shiftError = b.m_referenceResult[i].m_shiftError;
            ref.m_squeeze = b.m_referenceResult[i].m_squeeze;
            ref.m_squeezeError = b.m_referenceResult[i].m_squeezeError;
            ref.m_specieName = b.m_referenceResult[i].m_specieName;
            this->m_referenceResult[i] = ref;
        }
        memcpy(this->m_polynomial, b.m_polynomial, 5 * sizeof(float));

        this->m_chiSquare = b.m_chiSquare;
        this->m_delta = b.m_delta;
        this->m_stepNum = b.m_stepNum;
        this->m_evaluationStatus = b.m_evaluationStatus;
    }

    CEvaluationResult::~CEvaluationResult()
    {
    }

    // makes this a copy of 'b'
    CEvaluationResult &CEvaluationResult::operator =(const CEvaluationResult &b) {
        const int nRef = (int)b.m_referenceResult.size();
        m_referenceResult.resize(nRef);
        for (long i = 0; i < nRef; ++i) {
            CReferenceFitResult ref;
            ref.m_column = b.m_referenceResult[i].m_column;
            ref.m_columnError = b.m_referenceResult[i].m_columnError;
            ref.m_shift = b.m_referenceResult[i].m_shift;
            ref.m_shiftError = b.m_referenceResult[i].m_shiftError;
            ref.m_squeeze = b.m_referenceResult[i].m_squeeze;
            ref.m_squeezeError = b.m_referenceResult[i].m_squeezeError;
            ref.m_specieName = b.m_referenceResult[i].m_specieName;
            this->m_referenceResult[i] = ref;
        }
        memcpy(this->m_polynomial, b.m_polynomial, 5 * sizeof(float));

        this->m_chiSquare = b.m_chiSquare;
        this->m_delta = b.m_delta;
        this->m_stepNum = b.m_stepNum;

        m_evaluationStatus = b.m_evaluationStatus;

        return *this;
    }

    size_t CEvaluationResult::InsertSpecie(const std::string &name)
    {
        CReferenceFitResult ref;
        ref.m_specieName = std::string(name);
        m_referenceResult.push_back(ref);

        return m_referenceResult.size();
    }

    bool CEvaluationResult::CheckGoodnessOfFit(const CSpectrumInfo& info, float chi2Limit, float upperLimit, float lowerLimit)
    {
        // assume that this is an ok evaluation
        m_evaluationStatus &= ~MARK_BAD_EVALUATION;

        // The maximum intensity for one spectrum (# bits in the ADC)
        const double maxInt = CSpectrometerDatabase::GetInstance().GetModel(info.m_specModelName).maximumIntensity;

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
            else {
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
        else if (info.m_numSpec > 0) {
            offset = info.m_offset / (maxInt * info.m_numSpec);
        }
        else {
            int numSpec = (int)floor(info.m_peakIntensity / maxInt); // a guess for the number of co-adds
            offset = info.m_offset / (maxInt * numSpec);
        }

        // first check the intensity of the spectrum in the fit region
        if (upperLimit > -1) {
            if (fitSaturation > upperLimit)
                m_evaluationStatus |= MARK_BAD_EVALUATION;
        }
        else {
            if (fitSaturation > 0.99)
                m_evaluationStatus |= MARK_BAD_EVALUATION;
        }

        // first check the intensity of the spectrum in the fit region
        if (lowerLimit > -1) {
            if (fitSaturation < lowerLimit)
                m_evaluationStatus |= MARK_BAD_EVALUATION;
        }
        else {
            if ((fitSaturation - offset) < 0.025)
                m_evaluationStatus |= MARK_BAD_EVALUATION;
        }

        // then check the chi2 of the fit
        if (chi2Limit > -1) {
            if (m_chiSquare > chi2Limit)
                m_evaluationStatus |= MARK_BAD_EVALUATION;
        }
        else {
            if (m_chiSquare > 0.9)
                m_evaluationStatus |= MARK_BAD_EVALUATION;
        }

        return (m_evaluationStatus & MARK_BAD_EVALUATION);
    }

    bool CEvaluationResult::MarkAs(int MARK_FLAG) {
        // check the flag
        switch (MARK_FLAG) {
            case MARK_BAD_EVALUATION: break;
            case MARK_DELETED: break;
            default: return false;
        }

        // set the corresponding bit
        m_evaluationStatus |= MARK_FLAG;

        return true;
    }

    bool CEvaluationResult::RemoveMark(int MARK_FLAG) {
        // check the flag
        switch (MARK_FLAG) {
            case MARK_BAD_EVALUATION: break;
            case MARK_DELETED: break;
            default: return false;
        }

        // remove the corresponding bit
        m_evaluationStatus &= ~MARK_FLAG;

        return true;
    }

}
