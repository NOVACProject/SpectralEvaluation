#pragma once

#include "BasicMath.h"
#include "FitWindow.h"
#include "EvaluationResult.h"
#include "../Fit/ReferenceSpectrumFunction.h"

namespace Evaluation
{
    /** The CEvaluationBase is the base class for all evaluation-classes 
        in NovacProgram, NovacPPP and MobileDOAS and collects common elements and routines */
    class CEvaluationBase : public CBasicMath
    {
    public:
        CEvaluationBase();

        virtual ~CEvaluationBase();

        /** The fit window, defines the parameters for the fit */
        CFitWindow m_window;

        /** The result of the last performed evaluation. Defined only after Evaluate() has been called. */
        CEvaluationResult m_result;

        /** Removes the offset from the supplied spectrum */
        // TODO: Change the last parameter from begin a boolean to instead begin the pixel-range which should be used!!
        void RemoveOffset(double *spectrum, int sumChn, bool UV = true);

    protected:

        // The CReferenceSpectrumFunctions are used in the evaluation process to model
        // the reference spectra for the different species that are being fitted.
        CReferenceSpectrumFunction *ref[MAX_N_REFERENCES];
        CReferenceSpectrumFunction *solarSpec;

        /** Simple vector for holding the channel number information (element #i in this vector contains the value (i+1) */
        CVector vXData;

        /** Initializes the elements of the CReferenceSpectrumFunction-array 'ref' using the  information in m_window */
        int CreateReferenceSpectra();

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



    };
}