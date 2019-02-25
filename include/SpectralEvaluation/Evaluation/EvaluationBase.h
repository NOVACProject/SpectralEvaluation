#pragma once

#include "BasicMath.h"
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Evaluation/FitWindow.h>
#include "EvaluationResult.h"
#include "../Fit/ReferenceSpectrumFunction.h"

class CSpectrum;

namespace MathFit
{
    class CStandardFit;
}

namespace Evaluation
{
    /** The CEvaluationBase is the base class for all evaluation-classes 
        in NovacProgram, NovacPPP and MobileDOAS and collects common elements and routines */
    class CEvaluationBase : public CBasicMath
    {
    public:
        CEvaluationBase();

        // This object is not copyable due to the nature of its members.
        //  This could be implemented in the future if necessary
        CEvaluationBase(const CEvaluationBase&) = delete;
        CEvaluationBase& operator=(const CEvaluationBase&) = delete;

        explicit CEvaluationBase(const CFitWindow &window);

        virtual ~CEvaluationBase();

        /** Sets the fit window to use. This will initialize the 'ref' member variables. */
        void SetFitWindow(const CFitWindow &window);

        /** @return a reference to the local fit window */
        const CFitWindow& FitWindow() const { return m_window;}

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

        /** The last error from calling 'Evaluate' or 'EvaluateShift'.
            This is set by Evaluate and EvaluateShift and is only set if these methods return an error. */
        std::string m_lastError = "";

        /** Removes the offset from the supplied spectrum */
        // TODO: Change the last parameter from begin a boolean to instead begin the pixel-range which should be used!!
        void RemoveOffset(double *spectrum, int sumChn, bool UV = true);

        /** Sets the sky-spectrum to use. This will be used in the upcoming evaluations. 
            The provided spectrum must have been corrected for dark. */
        int SetSkySpectrum(const CSpectrum &spec);

        /** Evaluate using the parameters set in the local parameter 'm_window'
                and using the sky-spectrum which has been set by a previous call to 'SetSkySpectrum'
            The provided spectrum must have been corrected for dark.
            @return 0 if all is ok.
            @return 1 if any error occurred, or if the window is not defined. This will also set m_lastError. */
        int Evaluate(const CSpectrum& measured, int numSteps = 1000);

        /** Evaluate the optimum shift and squeeze to use for evaluating the measured spectrum.
            This is done by reading in the Fraunhofer reference spectrum in m_window.fraunhoferRef.
            If this Fraunhofer reference spectrum is sampled on the same grid as the references in m_window.reference
                then the returned shift and squeeze is the optimum shift and squeeze to use when evaluating the spectra.
            @param measured the spectrum for which to determine the shift & squeeze
                        relative to the solarReference-spectrum found in 'window'
            @return 0 if the fit succeeds and the shift & squeeze could be determined
            @return 1 if any error occured, see m_lastError for the error message. */
        int EvaluateShift(const CSpectrum &measured, double &shift, double &shiftError, double &squeeze, double &squeezeError);

        /** Returns the evaluation result for the last spectrum
               @return a reference to a 'CEvaluationResult' - data structure which holds the information from the last evaluation */
        const CEvaluationResult& GetEvaluationResult() const { return m_result; }

        /** Returns the polynomial that was fitted in the last evaluation */
        const double *GetPolynomial() const { return m_result.m_polynomial; }

        /** @return the number of references included in the fit. 
            This may be higher than the number of referene from the fit window if e.g. the sky-spectrum is included in the fit */
        size_t NumberOfReferencesFitted() const { return m_ref.size(); }

    protected:

        // The CReferenceSpectrumFunctions are used in the evaluation process to model
        // the reference spectra for the different species that are being fitted.
        //  The vector must hold pointers to the references, as these cannot be copied...
        std::vector<CReferenceSpectrumFunction*> m_ref;

        /** The sky spectrum to use in the evaluations. This is set by calling 'SetSkySpectrum' which must be called prior to calling 'Evaluate' */
        std::vector<double> m_sky;

        /** The fit window, defines the parameters for the fit */
        CFitWindow m_window;

        /** Simple vector for holding the channel number information (element #i in this vector contains the value (i+1) */
        CVector vXData;

        /** Initializes the elements of the CReferenceSpectrumFunction-array 'ref' using the information in m_window
            This must be called once prior to calling 'Evaluate', after the FitWindow has been set and the references read in from disk. */
        int CreateReferenceSpectra();

        /** Clears the elemnts in m_ref. This MUST be called before setting new references as the program will otherwise leak memory */
        void ClearRefereneSpectra();

        /** This generates the member 'vxData' and fills it with the values it should contain. */
        void CreateXDataVector(int numberOfChannels);

        // Prepares the spectra for evaluation
        void PrepareSpectra(double *sky, double *meas, const CFitWindow &window);

        // Prepares the spectra for evaluation
        void PrepareSpectra_HP_Div(double *sky, double *meas, const CFitWindow &window);

        // Prepares the spectra for evaluation
        void PrepareSpectra_HP_Sub(double *sky, double *meas, const CFitWindow &window);

        // Prepares the spectra for evaluation
        void PrepareSpectra_Poly(double *sky, double *meas, const CFitWindow &window);

        /** Updates the m_residual and m_result.delta */
        void SaveResidual(MathFit::CStandardFit& firstFit);
    };
}