#include "EvaluationBase.h"
#include "CrossSectionData.h"

namespace Evaluation
{
    CEvaluationBase::CEvaluationBase()
    {

    }

    CEvaluationBase::~CEvaluationBase()
    {

    }

    void CEvaluationBase::CreateXDataVector(int numberOfChannels)
    {
        vXData.SetSize(numberOfChannels);
        for (int i = 0; i < numberOfChannels; ++i)
        {
            vXData.SetAt(i, (TFitData)(1.0f + (double)i));
        }
    }

    int CEvaluationBase::CreateReferenceSpectra()
    {
        CVector yValues;

        for (int i = 0; i < m_window.nRef; i++)
        {
            // reset all reference's parameters
            ref[i]->ResetLinearParameter();
            ref[i]->ResetNonlinearParameter();

            // enable amplitude normalization. This should normally be done in order to avoid numerical
            // problems during fitting.
            ref[i]->SetNormalize(true);

            // set the spectral data of the reference spectrum to the object. This also causes an internal
            // transformation of the spectral data into a B-Spline that will be used to interpolate the 
            // reference spectrum during shift and squeeze operations
            //if(!ref[i]->SetData(vXData.SubVector(0, m_referenceData[i].GetSize()), m_referenceData[i]))
            yValues.SetSize(m_window.ref[i].m_data->GetSize());
            for (unsigned int k = 0; k < m_window.ref[i].m_data->GetSize(); ++k) {
                yValues.SetAt(k, m_window.ref[i].m_data->GetAt(k));
            }

            {
                auto tempXVec = vXData.SubVector(0, m_window.ref[i].m_data->GetSize());
                if (!ref[i]->SetData(tempXVec, yValues))
                {
                    Error0("Error initializing spline object!");
                    return(1);
                }
            }

            // Chech the options for the column value
            switch (m_window.ref[i].m_columnOption) {
            case SHIFT_FIX:   ref[i]->FixParameter(CReferenceSpectrumFunction::CONCENTRATION, m_window.ref[i].m_columnValue * ref[i]->GetAmplitudeScale()); break;
            case SHIFT_LINK:  ref[(int)m_window.ref[i].m_columnValue]->LinkParameter(CReferenceSpectrumFunction::CONCENTRATION, *ref[i], CReferenceSpectrumFunction::CONCENTRATION); break;
            }

            // Check the options for the shift
            switch (m_window.ref[i].m_shiftOption) {
            case SHIFT_FIX:		ref[i]->FixParameter(CReferenceSpectrumFunction::SHIFT, m_window.ref[i].m_shiftValue); break;
            case SHIFT_LINK:	ref[(int)m_window.ref[i].m_shiftValue]->LinkParameter(CReferenceSpectrumFunction::SHIFT, *ref[i], CReferenceSpectrumFunction::SHIFT); break;
            case SHIFT_LIMIT:	ref[i]->SetParameterLimits(CReferenceSpectrumFunction::SHIFT, (TFitData)m_window.ref[i].m_shiftValue, (TFitData)m_window.ref[i].m_shiftMaxValue, 1); break;
            default:			ref[i]->SetParameterLimits(CReferenceSpectrumFunction::SHIFT, (TFitData)-10.0, (TFitData)10.0, (TFitData)1e0); break;
                ref[i]->SetDefaultParameter(CReferenceSpectrumFunction::SHIFT, (TFitData)0.0);
            }

            // Check the options for the squeeze
            switch (m_window.ref[i].m_squeezeOption) {
            case SHIFT_FIX:		ref[i]->FixParameter(CReferenceSpectrumFunction::SQUEEZE, m_window.ref[i].m_squeezeValue); break;
            case SHIFT_LINK:	ref[(int)m_window.ref[i].m_squeezeValue]->LinkParameter(CReferenceSpectrumFunction::SQUEEZE, *ref[i], CReferenceSpectrumFunction::SQUEEZE); break;
            case SHIFT_LIMIT:	ref[i]->SetParameterLimits(CReferenceSpectrumFunction::SQUEEZE, (TFitData)m_window.ref[i].m_squeezeValue, (TFitData)m_window.ref[i].m_squeezeMaxValue, 1e7); break;
            default:			ref[i]->SetDefaultParameter(CReferenceSpectrumFunction::SQUEEZE, (TFitData)1.0);
                ref[i]->SetParameterLimits(CReferenceSpectrumFunction::SQUEEZE, (TFitData)0.98, (TFitData)1.02, (TFitData)1e0); break;
            }
        }

        // THIS FOLLOWING CODE COMES FROM THE NOVACPPP. IMPLEMENT HERE ALSO!!
        /* 
        // If we should also include the sky-spectrum in the fit
        if (m_skySpectrum.m_length > 0 && (m_window.fitType == FIT_HP_SUB || m_window.fitType == FIT_POLY)) {
            // reset all reference's parameters
            ref[i]->ResetLinearParameter();
            ref[i]->ResetNonlinearParameter();

            // enable amplitude normalization. This should normally be done in order to avoid numerical
            // problems during fitting.
            ref[i]->SetNormalize(true);

            // set the spectral data of the reference spectrum to the object. This also causes an internal
            // transformation of the spectral data into a B-Spline that will be used to interpolate the 
            // reference spectrum during shift and squeeze operations
            //if(!ref[i]->SetData(vXData.SubVector(0, m_referenceData[i].GetSize()), m_referenceData[i]))
            yValues.SetSize(m_skySpectrum.m_length);
            for (int k = 0; k < m_skySpectrum.m_length; ++k) {
                yValues.SetAt(k, m_sky[k]);
            }

            {
                auto tempXVec = vXData.SubVector(0, m_skySpectrum.m_length);
                if (!ref[i]->SetData(tempXVec, yValues))
                {
                    Error0("Error initializing spline object!");
                    return(1);
                }
            }

            // Chech the options for the column value
            if (m_window.fitType == FIT_POLY)
                ref[i]->FixParameter(CReferenceSpectrumFunction::CONCENTRATION, -1.0 * ref[i]->GetAmplitudeScale());
            else
                ref[i]->FixParameter(CReferenceSpectrumFunction::CONCENTRATION, 1.0 * ref[i]->GetAmplitudeScale());

            // Check the options for the shift & squeeze
            if (m_window.shiftSky) {
                ref[i]->SetParameterLimits(CReferenceSpectrumFunction::SHIFT, (TFitData)-3.0, (TFitData)3.0, 1);
                ref[i]->SetParameterLimits(CReferenceSpectrumFunction::SQUEEZE, (TFitData)0.95, (TFitData)1.05, 1e7);
            }
            else {
                ref[i]->FixParameter(CReferenceSpectrumFunction::SHIFT, 0.0);
                ref[i]->FixParameter(CReferenceSpectrumFunction::SQUEEZE, 1.0);
            }
        } */

        return 0;
    }

    void CEvaluationBase::RemoveOffset(double *spectrum, int sumChn, bool UV) {
        int offsetFrom = (UV) ? 50 : 2;
        int offsetTo = (UV) ? 200 : 20;

        //  remove any remaining offset in the measured spectrum
        double avg = 0;
        for (int i = offsetFrom; i < offsetTo; i++) {
            avg += spectrum[i];
        }
        avg = avg / (double)(offsetTo - offsetFrom);

        Sub(spectrum, sumChn, avg);

        return;
    }

    void CEvaluationBase::PrepareSpectra(double *sky, double *meas, const CFitWindow &window) {

        if (window.fitType == FIT_HP_DIV)
            return PrepareSpectra_HP_Div(sky, meas, window);
        if (window.fitType == FIT_HP_SUB)
            return PrepareSpectra_HP_Sub(sky, meas, window);
        if (window.fitType == FIT_POLY)
            return PrepareSpectra_Poly(sky, meas, window);
    }

    void CEvaluationBase::PrepareSpectra_HP_Div(double *skyArray, double *measArray, const CFitWindow &window) {

        //  1. Remove any remaining offset
        RemoveOffset(measArray, window.specLength, window.UV);
        RemoveOffset(skyArray, window.specLength, window.UV);

        // 2. Divide the measured spectrum with the sky spectrum
        Div(measArray, skyArray, window.specLength, 0.0);

        // 3. high pass filter
        HighPassBinomial(measArray, window.specLength, 500);

        // 4. log(spec)
        Log(measArray, window.specLength);

        // 5. low pass filter
        //	LowPassBinomial(measArray,window.specLength, 5);
    }

    void CEvaluationBase::PrepareSpectra_HP_Sub(double *skyArray, double *measArray, const CFitWindow &window) {

        // 1. remove any remaining offset in the measured spectrum
        RemoveOffset(measArray, window.specLength, window.UV);

        // 2. high pass filter
        HighPassBinomial(measArray, window.specLength, 500);

        // 3. log(spec)
        Log(measArray, window.specLength);
    }

    void CEvaluationBase::PrepareSpectra_Poly(double *skyArray, double *measArray, const CFitWindow &window) {

        // 1. remove any remaining offset in the measured spectrum
        RemoveOffset(measArray, window.specLength, window.UV);

        // 2. log(spec)
        Log(measArray, window.specLength);

        // 3. Multiply the spectrum with -1 to get the correct sign for everything
        for (int i = 0; i < window.specLength; ++i) {
            measArray[i] *= -1.0;
        }
    }

}