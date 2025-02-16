#pragma once

#include <SpectralEvaluation/Log.h>
#include <SpectralEvaluation/Fit/ReferenceSpectrumFunction.h>
#include <SpectralEvaluation/Evaluation/BasicMath.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Evaluation/FitWindow.h>
#include <SpectralEvaluation/Evaluation/EvaluationResult.h>

namespace MathFit
{
class CStandardFit;
}

namespace novac
{
class CSpectrum;

struct ShiftEvaluationResult
{
    double shift = 0.0;
    double shiftError = 0.0;
    double squeeze = 0.0;
    double squeezeError = 0.0;
    double chi2 = 0.0;
};

/** The CEvaluationBase is the base class for all evaluation-classes
    in NovacProgram, NovacPPP and MobileDOAS and collects common elements and routines */
class CEvaluationBase : public CBasicMath
{
public:
    CEvaluationBase(novac::ILogger& log);

    // This object is not copyable due to the nature of its members.
    //  This could be implemented in the future if necessary
    CEvaluationBase(const CEvaluationBase&) = delete;
    CEvaluationBase& operator=(const CEvaluationBase&) = delete;

    explicit CEvaluationBase(const CFitWindow& window, novac::ILogger& log);

    virtual ~CEvaluationBase();

    /** Sets the fit window to use. This will initialize the 'ref' member variables. */
    void SetFitWindow(const CFitWindow& window);

    /** @return a reference to the local fit window */
    const CFitWindow& FitWindow() const { return m_window; }

    /** The result of the last performed evaluation.
        Defined only after Evaluate() or EvaluateShift() have been called. */
    CEvaluationResult m_result;

    /** The residual of the last performed evaluation.
        Defined only after Evaluate() or EvaluateShift() have been called. */
    CVector m_residual;

    /** The scaled reference spectra.
        Defined only after Evaluate() or EvaluateShift() have been called.
        The first reference spectrum is the fitted
            polynomial and the following 'MAX_N_REFERENCES + 1' are the scaled references. */
            // TODO: Implement saving this
    CCrossSectionData m_fitResult[MAX_N_REFERENCES + 2];

    /** The measured spectrum, after all processing is done, right before the fit is performed. */
    std::vector<double> m_measuredData;

    /** The last error from calling 'Evaluate' or 'EvaluateShift'.
        This is set by Evaluate and EvaluateShift and is only set if these methods return an error. */
    std::string m_lastError = "";

    /** Removes the offset from the supplied spectrum */
    // TODO: Change the last parameter from a boolean to the pixel-range which should be used!!
    void RemoveOffset(double* spectrum, int sumChn, IndexRange range) const;

    /** Removes the offset from the supplied spectrum.
        The offset is calculated in the [from, to] region but the offset is subtracted from the entire spectrum. */
    void RemoveOffset(std::vector<double>& spectrum, IndexRange range) const;

    /** Sets the sky-spectrum to use. This will be used in the upcoming evaluations.
        The provided spectrum must have been corrected for dark. */
    int SetSkySpectrum(const CSpectrum& spec);
    int SetSkySpectrum(const CCrossSectionData& spec, bool removeOffset = false);

    /** Evaluate using the parameters set in the local parameter 'm_window'
            and using the sky-spectrum which has been set by a previous call to 'SetSkySpectrum'
        The provided spectrum must have been corrected for dark.
        @return 0 if all is ok.
        @return 1 if any error occurred, or if the window is not defined. This will also set m_lastError. */
    int Evaluate(const CSpectrum& measured, int numSteps = 1000);

    /** Evaluate using the parameters set in the local parameter 'm_window'
            and using the sky-spectrum which has been set by a previous call to 'SetSkySpectrum'
        The provided spectrum must have been corrected for dark.
        @return 0 if all is ok.
        @return 1 if any error occurred, or if the window is not defined. This will also set m_lastError. */
    int Evaluate(const double* measured, size_t measuredLength, int measuredStartChannel = 0, int numSteps = 1000);

    /** Evaluate the optimum shift and squeeze to use for evaluating the measured spectrum.
        This is done by reading in the Fraunhofer reference spectrum in m_window.fraunhoferRef.
        If this Fraunhofer reference spectrum is sampled on the same grid as the references in m_window.reference
            then the returned shift and squeeze is the optimum shift and squeeze to use when evaluating the spectra.
        @param measured the spectrum for which to determine the shift & squeeze
                    relative to the solarReference-spectrum found in 'window'
        @return 0 if the fit succeeds and the shift & squeeze could be determined
        @return 1 if any error occured, see m_lastError for the error message. */
    int EvaluateShift(novac::LogContext context, const CSpectrum& measured, ShiftEvaluationResult& result);

    /** Returns the evaluation result for the last spectrum
           @return a reference to a 'CEvaluationResult' - data structure which holds the information from the last evaluation */
    const CEvaluationResult& GetEvaluationResult() const { return m_result; }

    /** Returns the polynomial that was fitted in the last evaluation */
    // TODO: Change to std::vector<double>
    const double* GetPolynomial() const { return m_result.m_polynomial; }

    /** @return the number of references included in the fit.
        This may be higher than the number of referene from the fit window if
        the sky-spectrum is included in the fit or an offset polynomial is included */
    size_t NumberOfReferencesFitted() const { return m_ref.size(); }

    // The CReferenceSpectrumFunctions are used in the evaluation process to model
    // the reference spectra for the different species that are being fitted.
    //  The vector must hold pointers to the references, as these cannot be copied...
    std::vector<CReferenceSpectrumFunction*> m_ref;

    // If the sky-spectrum is to be fitted in the evaluation, then this represents
    //  the reference-spectrum function of the sky spectrum.
    CReferenceSpectrumFunction* m_skyReference = nullptr;

    // If an polynomial is to be included in the evaluation in intensity-space 
    //  then this represents the reference-spectrum function of the polynomial.
    CReferenceSpectrumFunction* m_intensitySpacePolynomial = nullptr;

    // TODO: Combine these loose CReferenceSpectrumFunction into a structure!
    // If a calculated ring spectrum is to be included in the evaluation 
    //  then this represents the reference-spectrum of that ring.
    CReferenceSpectrumFunction* m_ringSpectrum = nullptr;

    // If a calculated ring spectrum scaled by lamda^4 is to be included in the evaluation 
    //  then this represents the reference-spectrum of that ring.
    CReferenceSpectrumFunction* m_ringSpectrumLambda4 = nullptr;

protected:

    novac::ILogger& m_log;

    /** The sky spectrum to use in the evaluations.
        This is set by calling 'SetSkySpectrum' which must be called prior to calling 'Evaluate' */
    CCrossSectionData m_sky;

    /** The fit window, defines the parameters for the fit */
    CFitWindow m_window;

    /** Simple vector for holding the channel number information (element #i in this vector contains the value (i+1) */
    CVector vXData;

    /** Initializes the elements of the CReferenceSpectrumFunction-array 'ref' using the information in m_window
        This must be called once prior to calling 'Evaluate', after the FitWindow has been set and the references read in from disk. */
    int CreateReferenceSpectra();

    /** Clears the elemnts in m_ref. This MUST be called before setting new references as the program will otherwise leak memory */
    void ClearRefereneSpectra();

    /** This generates the member vector 'vxData' which is used to define the x-axis
        values in the DOAS fit. Here this is filled with the channel indices [0, numberOfChannels[ . */
    void CreateXDataVector(int numberOfChannels);

    // Prepares the spectra for evaluation
    [[deprecated]]
    void PrepareSpectra(double* sky, double* meas, const CFitWindow& window);

    // Prepares the spectra for evaluation
    void PrepareSpectra(std::vector<double>& sky, std::vector<double>& meas, const CFitWindow& window) const;

    void PrepareSpectra_HP_Div(std::vector<double>& sky, std::vector<double>& meas, const CFitWindow& window) const;

    void PrepareSpectra_HP_Sub(std::vector<double>& sky, std::vector<double>& meas, const CFitWindow& window) const;

    void PrepareSpectra_Poly(std::vector<double>& sky, std::vector<double>& meas, const CFitWindow& window) const;

    /** Updates the m_residual and m_result.delta */
    void SaveResidual(MathFit::CStandardFit& firstFit);

    // This sets up the member 'm_intensitySpacePolynomial'
    void CreateReferenceForIntensitySpacePolynomial(const std::vector<double>& I0);

    // This sets up the member 'm_skyReference'
    void CreateReferenceForSkySpectrum();

    // This sets up the member 'm_ringSpectrum'
    void CreateReferenceForRingSpectrum(const CSpectrum& ring);

    // This sets up the member 'm_ringSpectrumLambda4'
    void CreateReferenceForRingSpectrumLambda4(const CSpectrum& ring);

    // @return the name of the reference with the given index into m_ref
    std::string GetReferenceName(size_t referenceIndex) const;

private:
    // The following method is deprecated and should be removed eventually.
    void PrepareSpectra_HP_Div(double* sky, double* meas, const CFitWindow& window);

    // The following method is deprecated and should be removed eventually.
    void PrepareSpectra_HP_Sub(double* sky, double* meas, const CFitWindow& window);

    // The following method is deprecated and should be removed eventually.
    void PrepareSpectra_Poly(double* sky, double* meas, const CFitWindow& window);

};
}