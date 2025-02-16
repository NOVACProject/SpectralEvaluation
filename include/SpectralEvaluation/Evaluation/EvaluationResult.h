#pragma once

#include <vector>
#include <SpectralEvaluation/Evaluation/ReferenceFitResult.h>

namespace novac
{

class CSpectrumInfo;
struct DoasResult;
struct SpectrometerModel;

#define MARK_BAD_EVALUATION 0x0001
#define MARK_DELETED 0x0002

/** CEvaluationResult is a storage container for the results after evaluating one single spectrum.
        It contains an array of CReferenceFitResult:s - one for each
        reference spectrum that was included in the fit, and other information about the fit -
        such as Delta, ChiSquare, the number of iterations required and the fitted polynomial. */
class CEvaluationResult
{
public:
    CEvaluationResult();

    CEvaluationResult(const CEvaluationResult& b);
    CEvaluationResult& operator=(const CEvaluationResult& b);

    // Converts a DoasResult into a CEvaluationResult, these two classes are just two different ways of storing the result of a DOAS fit. 
    CEvaluationResult(const DoasResult& b);

    /** The evaluated result for the reference files */
    std::vector<CReferenceFitResult> m_referenceResult;

    /** The result for the polynomial.
        The polynomial is stored so that m_polynomial[i] is the i:th order term */
    double m_polynomial[6];

    /** The delta of the fit (peak-to-peak value of the residual) */
    double m_delta = 0.0;

    /** The Chi-square of the fit */
    double m_chiSquare = 0.0;

    /** The number of steps required to make the fitting */
    long m_stepNum = 0;

    /** The status of the evaluation.
        if(m_evaluationStatus & BAD_EVALUATION) then the spectrum is marked as a bad evaluation
        if(m_evaluationStatus & DELETED) then the spectrum is marked as deleted (used in post-flux calculations) */
    int m_evaluationStatus = 0;

    // --------------- PUBLIC METHODS ---------------------
    size_t NumberOfSpecies() const { return m_referenceResult.size(); }

    /** Increases the list of references evaluated by inserting a
        new CReferenceFitResult into m_referenceResult and giving it the provided name.
        @return the new number of references. */
    size_t InsertSpecie(const std::string& name);

    /** Checks the goodness of fit for this spectrum.
        @param info - the information about this spectrum
        @param chi2Limit - (optional) the limit of chi2 to use
        @param upperLimit - (optional) the upper limit to use on the saturation level
        @param lowerLimit - (optional) the lower limit to use on the saturation level
        @return true if this spectrum result is evaluated to be good, otherwise returns false. */
    bool CheckGoodnessOfFit(const CSpectrumInfo& info, const SpectrometerModel* spectrometer = nullptr, float chi2Limit = -1, float upperLimit = -1, float lowerLimit = -1);

    /** Returns true if this spectrum is judges as being an ok spectrum */
    bool IsOK() const { return !(m_evaluationStatus & MARK_BAD_EVALUATION); }

    /** Returns false if this spectrum is judges as being a bad spectrum */
    bool IsBad() const { return (m_evaluationStatus & MARK_BAD_EVALUATION); }

    /** Returns true if this spectrum is judges as being an ok spectrum */
    bool IsDeleted() const { return (m_evaluationStatus & MARK_DELETED); }

    /** Marks the current spectrum with the supplied mark_flag.
        Mark flag must be MARK_BAD_EVALUATION, or MARK_DELETED
        @return true on success (i.e. the flag to add is valid). */
    bool MarkAs(int MARK_FLAG);

    /** Removes the current mark from the desired spectrum
        Mark flag must be MARK_BAD_EVALUATION, or MARK_DELETED
        @return true on success (i.e. the flag to remove is valid). */
    bool RemoveMark(int MARK_FLAG);

};
}