#include "EvaluationBase.h"
#include "CrossSectionData.h"
#include "../Spectra/Spectrum.h"

// include all required fit objects
#include "../Fit/ReferenceSpectrumFunction.h"
#include "../Fit/SimpleDOASFunction.h"
#include "../Fit/StandardMetricFunction.h"
#include "../Fit/StandardFit.h"
#include "../Fit/ExpFunction.h"
#include "../Fit/LnFunction.h"
#include "../Fit/PolynomialFunction.h"
#include "../Fit/NegateFunction.h"
#include "../Fit/MulFunction.h"
#include "../Fit/DivFunction.h"
#include "../Fit/GaussFunction.h"
#include "../Fit/DiscreteFunction.h"
#include "../Fit/DOASVector.h"
#include "../Fit/NonlinearParameterFunction.h"

using namespace MathFit;

namespace Evaluation
{
    CEvaluationBase::CEvaluationBase()
    {
        CreateXDataVector(MAX_SPECTRUM_LENGTH);
    }

    CEvaluationBase::CEvaluationBase(const CFitWindow &window)
    {
        this->m_window = window;
        CreateXDataVector(MAX_SPECTRUM_LENGTH);
        CreateReferenceSpectra();
    }

    CEvaluationBase::~CEvaluationBase()
    {
        ClearRefereneSpectra();
    }

    void CEvaluationBase::SetFitWindow(const CFitWindow &window)
    {
        this->m_window = window;
        CreateReferenceSpectra();
    }

    void CEvaluationBase::CreateXDataVector(int numberOfChannels)
    {
        vXData.SetSize(numberOfChannels);
        for (int i = 0; i < numberOfChannels; ++i)
        {
            vXData.SetAt(i, (TFitData)(1.0f + (double)i));
        }
    }

    void CEvaluationBase::ClearRefereneSpectra()
    {
        for (auto* r : m_ref)
        {
            delete r;
        }
        m_ref.clear();
    }

    int CEvaluationBase::CreateReferenceSpectra()
    {
        CVector yValues;

        ClearRefereneSpectra();

        // 1) Create the references
        for (int i = 0; i < m_window.nRef; i++)
        {
            CReferenceSpectrumFunction* newRef = new CReferenceSpectrumFunction();

            // reset all reference's parameters
            newRef->ResetLinearParameter();
            newRef->ResetNonlinearParameter();

            // enable amplitude normalization. This should normally be done in order to avoid numerical
            // problems during fitting.
            newRef->SetNormalize(true);

            // set the spectral data of the reference spectrum to the object. This also causes an internal
            // transformation of the spectral data into a B-Spline that will be used to interpolate the 
            // reference spectrum during shift and squeeze operations
            //if(!newRef->SetData(vXData.SubVector(0, m_referenceData[i].GetSize()), m_referenceData[i]))
            yValues.SetSize(m_window.ref[i].m_data->GetSize());
            for (unsigned int k = 0; k < m_window.ref[i].m_data->GetSize(); ++k)
            {
                yValues.SetAt(k, m_window.ref[i].m_data->GetAt(k));
            }

            {
                auto tempXVec = vXData.SubVector(0, m_window.ref[i].m_data->GetSize());
                if (!newRef->SetData(tempXVec, yValues))
                {
                    Error0("Error initializing spline object!");
                    return(1);
                }
            }

            // Finally add this reference to the vector
            this->m_ref.push_back(newRef);
        }


        // 2) Couple the references
        for (int i = 0; i < m_window.nRef; i++)
        {
            // Chech the options for the column value
            switch (m_window.ref[i].m_columnOption) {
            case SHIFT_FIX:     m_ref[i]->FixParameter(CReferenceSpectrumFunction::CONCENTRATION, m_window.ref[i].m_columnValue * m_ref[i]->GetAmplitudeScale()); break;
            case SHIFT_LINK:    m_ref[(int)m_window.ref[i].m_columnValue]->LinkParameter(CReferenceSpectrumFunction::CONCENTRATION, *m_ref[i], CReferenceSpectrumFunction::CONCENTRATION); break;
            }

            // Check the options for the shift
            switch (m_window.ref[i].m_shiftOption) {
            case SHIFT_FIX:     m_ref[i]->FixParameter(CReferenceSpectrumFunction::SHIFT, m_window.ref[i].m_shiftValue); break;
            case SHIFT_LINK:    m_ref[(int)m_window.ref[i].m_shiftValue]->LinkParameter(CReferenceSpectrumFunction::SHIFT, *m_ref[i], CReferenceSpectrumFunction::SHIFT); break;
            case SHIFT_LIMIT:   m_ref[i]->SetParameterLimits(CReferenceSpectrumFunction::SHIFT, (TFitData)m_window.ref[i].m_shiftValue, (TFitData)m_window.ref[i].m_shiftMaxValue, 1); break;
            default:            m_ref[i]->SetDefaultParameter(CReferenceSpectrumFunction::SHIFT, (TFitData)0.0);
                                m_ref[i]->SetParameterLimits(CReferenceSpectrumFunction::SHIFT, (TFitData)-10.0, (TFitData)10.0, (TFitData)1e0); break; // TODO: Get these limits as parameters!
            }

            // Check the options for the squeeze
            switch (m_window.ref[i].m_squeezeOption) {
            case SHIFT_FIX:     m_ref[i]->FixParameter(CReferenceSpectrumFunction::SQUEEZE, m_window.ref[i].m_squeezeValue); break;
            case SHIFT_LINK:    m_ref[(int)m_window.ref[i].m_squeezeValue]->LinkParameter(CReferenceSpectrumFunction::SQUEEZE, *m_ref[i], CReferenceSpectrumFunction::SQUEEZE); break;
            case SHIFT_LIMIT:   m_ref[i]->SetParameterLimits(CReferenceSpectrumFunction::SQUEEZE, (TFitData)m_window.ref[i].m_squeezeValue, (TFitData)m_window.ref[i].m_squeezeMaxValue, 1e7); break;
            default:            m_ref[i]->SetDefaultParameter(CReferenceSpectrumFunction::SQUEEZE, (TFitData)1.0);
                                m_ref[i]->SetParameterLimits(CReferenceSpectrumFunction::SQUEEZE, (TFitData)0.98, (TFitData)1.02, (TFitData)1e0); break; // TODO: Get these limits as parameters!
            }
        }

        // THIS FOLLOWING CODE COMES FROM THE NOVACPPP. IMPLEMENT HERE ALSO!!
        // If we should also include the sky-spectrum in the fit
        if (m_sky.size() > 0 && (m_window.fitType == FIT_HP_SUB || m_window.fitType == FIT_POLY))
        {
            CReferenceSpectrumFunction* newRef = new CReferenceSpectrumFunction();

            // reset all reference's parameters
            newRef->ResetLinearParameter();
            newRef->ResetNonlinearParameter();

            // enable amplitude normalization. This should normally be done in order to avoid numerical
            // problems during fitting.
            newRef->SetNormalize(true);

            // set the spectral data of the reference spectrum to the object. This also causes an internal
            // transformation of the spectral data into a B-Spline that will be used to interpolate the
            // reference spectrum during shift and squeeze operations
            //if(!ref[i]->SetData(vXData.SubVector(0, m_referenceData[i].GetSize()), m_referenceData[i]))
            yValues.SetSize((int)m_sky.size());
            for (int k = 0; k < (int)m_sky.size(); ++k) {
                yValues.SetAt(k, m_sky[k]);
            }

            {
                auto tempXVec = vXData.SubVector(0, (int)m_sky.size());
                if (!newRef->SetData(tempXVec, yValues))
                {
                    Error0("Error initializing spline object!");
                    return(1);
                }
            }

            // Chech the options for the column value
            if (m_window.fitType == FIT_POLY)
                newRef->FixParameter(CReferenceSpectrumFunction::CONCENTRATION, -1.0 * newRef->GetAmplitudeScale());
            else
                newRef->FixParameter(CReferenceSpectrumFunction::CONCENTRATION, 1.0 * newRef->GetAmplitudeScale());

            // Check the options for the shift & squeeze
            if (m_window.shiftSky) {
                newRef->SetParameterLimits(CReferenceSpectrumFunction::SHIFT, (TFitData)-3.0, (TFitData)3.0, 1);
                newRef->SetParameterLimits(CReferenceSpectrumFunction::SQUEEZE, (TFitData)0.95, (TFitData)1.05, 1e7);
            }
            else {
                newRef->FixParameter(CReferenceSpectrumFunction::SHIFT, 0.0);
                newRef->FixParameter(CReferenceSpectrumFunction::SQUEEZE, 1.0);
            }

            // Finally add this reference to the vector
            this->m_ref.push_back(newRef);
        }

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

        // 2. Divide the measured spectrum with the sky spectrum
        Div(measArray, skyArray, window.specLength, 0.0);

        // 3. high pass filter
        HighPassBinomial(measArray, window.specLength, 500);

        // 4. log(spec)
        Log(measArray, window.specLength);

        // 5. low pass filter
        //	LowPassBinomial(measArray,window.specLength, 5);
    }

    void CEvaluationBase::PrepareSpectra_HP_Sub(double* /*skyArray*/, double *measArray, const CFitWindow &window) {

        // 1. remove any remaining offset in the measured spectrum
        RemoveOffset(measArray, window.specLength, window.UV);

        // 2. high pass filter
        HighPassBinomial(measArray, window.specLength, 500);

        // 3. log(spec)
        Log(measArray, window.specLength);
    }

    void CEvaluationBase::PrepareSpectra_Poly(double* /*skyArray*/, double *measArray, const CFitWindow &window) {

        // 1. remove any remaining offset in the measured spectrum
        RemoveOffset(measArray, window.specLength, window.UV);

        // 2. log(spec)
        Log(measArray, window.specLength);

        // 3. Multiply the spectrum with -1 to get the correct sign for everything
        for (int i = 0; i < window.specLength; ++i) {
            measArray[i] *= -1.0;
        }
    }

    int CEvaluationBase::SetSkySpectrum(const CSpectrum &spec)
    {
        // make sure that these make sense
        if (m_window.specLength != spec.m_length)
        {
            return 1;
        }

        // set the data
        m_sky = std::vector<double>(spec.m_data, spec.m_data + spec.m_length);

        // --------- prepare the sky spectrum for evaluation ---------
        // Remove any remaining offset of the sky-spectrum
        RemoveOffset(m_sky.data(), (int)m_sky.size(), m_window.UV);

        // High-pass filter the sky-spectrum
        if (m_window.fitType == FIT_HP_SUB)
        {
            HighPassBinomial(m_sky.data(), (int)m_sky.size(), 500);
        }

        // Logaritmate the sky-spectrum
        if (m_window.fitType != FIT_HP_DIV)
        {
            Log(m_sky.data(), (int)m_sky.size());
        }

        if(m_window.fitType != FIT_HP_DIV)
        {
            // Include the sky-spectrum as a reference into the fit
            CreateReferenceSpectra();
        }

        return 0;
    }

    void CEvaluationBase::SaveResidual(CStandardFit& cFirstFit)
    {
        //// get residuum vector and expand it to a DOAS vector object. Do NOT assign the vector data to the new object!
        //// display some statistical stuff about the residual data
        CDOASVector vResiduum;
        vResiduum.Attach(cFirstFit.GetResiduum(), false);
        m_residual.SetSize(vResiduum.GetSize());
        m_residual.Zero();
        m_residual.Add(vResiduum);

        m_result.m_delta = (double)vResiduum.Delta();
    }

    int CEvaluationBase::Evaluate(const CSpectrum &measured, int numSteps)
    {
        assert(vXData.GetSize() >= measured.m_length);

        m_lastError = "";

        int fitLow, fitHigh; // the limits for the DOAS fit

        // Check so that the length of the spectra agree with each other
        if (m_window.specLength != measured.m_length)
        {
            m_lastError = "Failed to evaluate: the length of the measured spectrum does not equal the spectrum length of the fit-window used.";
            return 1;
        }

        // Get the correct limits for the fit
        if (measured.m_info.m_startChannel != 0 || measured.m_length < m_window.specLength)
        {
            // Partial spectra
            fitLow = m_window.fitLow - measured.m_info.m_startChannel;
            fitHigh = m_window.fitHigh - measured.m_info.m_startChannel;
            if (fitLow < 0 || fitHigh > measured.m_length)
            {
                m_lastError = "Invalid fit region, fit-low is negative or fit-high is larger than the length of the measured spectrum.";
                return 1;
            }
        }
        else
        {
            fitLow = m_window.fitLow;
            fitHigh = m_window.fitHigh;
        }

        // Vectors to store the data
        CVector vMeas, vSky;

        // Make a local copy of the data (since we're going to change the contents)
        std::vector<double> measArray(measured.m_data, measured.m_data + measured.m_length );

        //----------------------------------------------------------------
        // --------- prepare the spectrum for evaluation -----------------
        //----------------------------------------------------------------

        PrepareSpectra(m_sky.data(), measArray.data(), m_window);

        //----------------------------------------------------------------

        // Copy the measured spectrum to vMeas
        vMeas.Copy(measArray.data(), m_window.specLength, 1);

        // To perform the fit we need to extract the wavelength (or pixel)
        //	information from the vXData-vector
        CVector vXSec(fitHigh - fitLow);
        vXSec.Copy(vXData.SubVector(fitLow, fitHigh - fitLow));

        ////////////////////////////////////////////////////////////////////////////
        // now we start building the model function needed for fitting.
        //
        // First we create a function object that represents our measured spectrum. Since we do not
        // need any interpolation on the measured data its enough to use a CDiscreteFunction object.
        CDiscreteFunction dataTarget;

        // now set the data of the measured spectrum in regard to the wavelength information
        {
            auto temp = vXData.SubVector(measured.m_info.m_startChannel, m_window.specLength);
            dataTarget.SetData(temp, vMeas);
        }

        // since the DOAS model function consists of the sum of all reference spectra and a polynomial,
        // we first create a summation object
        CSimpleDOASFunction cRefSum;

        // now we add the required CReferenceSpectrumFunction objects that actually represent the 
        // reference spectra used in the DOAS model function
        for (CReferenceSpectrumFunction* reference : this->m_ref)
        {
            cRefSum.AddReference(*reference); // <-- at last add the reference to the summation object
        }

        // create the additional polynomial with the correct order
        //	and add it to the summation object, too
        CPolynomialFunction cPoly(m_window.polyOrder);
        cRefSum.AddReference(cPoly);

        // the last step in the model function will be to define how the difference between the measured data and the modeled
        // data will be determined. In this case we will use the CStandardMetricFunction which actually just calculate the difference
        // between the measured data and the modeled data channel by channel. The fit will try to minimize these differences.
        // So we create the metric object and set the measured spectrum function object and the DOAS model function object as parameters
        CStandardMetricFunction cDiff(dataTarget, cRefSum);

        /////////////////////////////////////////////////////////////////
        // Now its time to create the fit object. The CStandardFit object will 
        // provide a combination of a linear Least Square Fit and a nonlinear Levenberg-Marquardt Fit, which
        // should be sufficient for most needs.
        CStandardFit cFirstFit(cDiff);

        // don't forget to the the already extracted fit range to the fit object!
        // without a valid fit range you'll get an exception.
        cFirstFit.SetFitRange(vXSec);

        // limit the number of fit iteration to 5000. This can still take a long time! More convinient values are
        // between 100 and 1000
        cFirstFit.GetNonlinearMinimizer().SetMaxFitSteps(numSteps);
        cFirstFit.GetNonlinearMinimizer().SetMinChiSquare(0.0001);

        try
        {
            // prepare everything for fitting
            cFirstFit.PrepareMinimize();

            // actually do the fitting
            if (!cFirstFit.Minimize())
            {
                // message.Format("Fit Failed!");
                // ShowMessage(message);
                m_lastError = "Failed to evaluate: fit failed.";
                return 1;
            }

            // finalize the fitting process. This will calculate the error measurements and other statistical stuff
            cFirstFit.FinishMinimize();

            // get the basic fit results
            m_result.m_stepNum   = (long)cFirstFit.GetFitSteps();
            m_result.m_chiSquare = (double)cFirstFit.GetChiSquare();
            m_result.m_referenceResult.resize(this->m_ref.size());

            for (int tmpInt = 0; tmpInt <= m_window.polyOrder; ++tmpInt)
            {
                m_result.m_polynomial[tmpInt] = (double)cPoly.GetCoefficient(tmpInt);
            }

            SaveResidual(cFirstFit);

            // get the fitResult for the polynomial
            CVector tmpVector;
            tmpVector.SetSize(fitHigh - fitLow);
            auto tempXVec = vXData.SubVector(fitLow, fitHigh - fitLow);
            cPoly.GetValues(tempXVec, tmpVector);
            m_fitResult[0].Set(tmpVector, fitHigh - fitLow);

            // finally display the fit results for each reference spectrum including their appropriate error
            for (size_t ii = 0; ii < this->m_ref.size(); ii++)
            {
                m_result.m_referenceResult[ii].m_specieName      = (ii < (size_t)m_window.nRef) ? std::string(m_window.ref[ii].m_specieName) : "Sky";
                m_result.m_referenceResult[ii].m_column          = (double)m_ref[ii]->GetModelParameter(CReferenceSpectrumFunction::CONCENTRATION);
                m_result.m_referenceResult[ii].m_columnError     = (double)m_ref[ii]->GetModelParameterError(CReferenceSpectrumFunction::CONCENTRATION);
                m_result.m_referenceResult[ii].m_shift           = (double)m_ref[ii]->GetModelParameter(CReferenceSpectrumFunction::SHIFT);
                m_result.m_referenceResult[ii].m_shiftError      = (double)m_ref[ii]->GetModelParameterError(CReferenceSpectrumFunction::SHIFT);
                m_result.m_referenceResult[ii].m_squeeze         = (double)m_ref[ii]->GetModelParameter(CReferenceSpectrumFunction::SQUEEZE);
                m_result.m_referenceResult[ii].m_squeezeError    = (double)m_ref[ii]->GetModelParameterError(CReferenceSpectrumFunction::SQUEEZE);

                //// get the final fit result
                CVector tmpVector;
                tmpVector.SetSize(fitHigh - fitLow);
                auto tempXVec = vXData.SubVector(fitLow, fitHigh - fitLow);
                m_ref[ii]->GetValues(tempXVec, tmpVector);
                m_fitResult[ii + 1].Set(tmpVector, fitHigh - fitLow);
            }

            return 0;
        }
        catch (CFitException& e)
        {
            // in case that something went wrong, display the error to the user.
            // normally you will get error in two cases:
            //
            // 1. You forgot to set a valid fit range before you start fitting
            //
            // 2. A matrix inversion failed for some reason inside the fitting loop. Matrix inversions
            //    normally fail when there are linear dependecies in the matrix respectrively you have linear
            //    dependencies in your reference spectrum. Eg. you tried to fit the same reference spectrum twice at once.

            //  e.ReportError();
            //	std::cout << "Failed: " << ++iFalseCount << std::endl;
            //	std::cout << "Steps: " << cFirstFit.GetNonlinearMinimizer().GetFitSteps() << " - Chi: " << cFirstFit.GetNonlinearMinimizer().GetChiSquare() << std::endl;

            // message.Format("A Fit Exception has occurred. Are the reference files OK?");
            // ShowMessage(message);
            m_lastError = "Failed to evaluate: a fit exception occurred.";
            if (nullptr != e.mMessage && strlen(e.mMessage) > 0)
            {
                m_lastError = m_lastError + " Message: '" + std::string{e.mMessage} + "'";
            }

            return (1);
        }
    }

    int CEvaluationBase::EvaluateShift(const CSpectrum &measured, double &shift, double &shiftError, double &squeeze, double &squeezeError)
    {
        assert(vXData.GetSize() >= measured.m_length);

        m_lastError = "";

        int i;
        CVector vMeas;
        CVector yValues;
        int fitLow = m_window.fitLow;
        int fitHigh = m_window.fitHigh;

        // Check so that the length of the spectra agree with each other
        if (m_window.specLength != measured.m_length)
        {
            m_lastError = "Failed to evaluate shift: the length of the measured spectrum does not equal the spectrum length of the fit-window used.";
            return 1;
        }

        // Check that we have a solar-spectrum to check against
        if (m_window.fraunhoferRef.m_path.size() < 6)
        {
            m_lastError = "Failed to evaluate shift: the fraunhofer reference does not exist.";
            return 1;
        }

        if (measured.m_info.m_startChannel != 0 || measured.m_length < m_window.specLength) {
            // Partial spectra
            fitLow = m_window.fitLow - measured.m_info.m_startChannel;
            fitHigh = m_window.fitHigh - measured.m_info.m_startChannel;
            if (fitLow < 0 || fitHigh > measured.m_length)
            {
                m_lastError = "Invalid fit region, fit-low is negative or fit-high is larger than the length of the measured spectrum.";
                return 1;
            }
        }
        else {
            fitLow = m_window.fitLow;
            fitHigh = m_window.fitHigh;
        }

        // initialize the solar-spectrum function
        CReferenceSpectrumFunction *solarSpec = new CReferenceSpectrumFunction();

        // Make a local copy of the data
        double *measArray = (double *)calloc(measured.m_length, sizeof(double));
        memcpy(measArray, measured.m_data, measured.m_length * sizeof(double));

        //----------------------------------------------------------------
        // --------- prepare the spectrum for evaluation -----------------
        //----------------------------------------------------------------

        RemoveOffset(measArray, m_window.specLength, m_window.UV);
        if (m_window.fitType == FIT_HP_DIV || m_window.fitType == FIT_HP_SUB) {
            HighPassBinomial(measArray, m_window.specLength, 500);
        }
        Log(measArray, m_window.specLength);

        if (m_window.fitType == FIT_POLY) {
            for (int j = 0; j < m_window.specLength; ++j)
                measArray[j] *= -1.0;
        }

        // --------- also prepare the solar-spectrum for evaluation -----------------
        CVector localSolarSpectrumData;
        localSolarSpectrumData.SetSize(m_window.fraunhoferRef.m_data->GetSize());
        for (int j = 0; j < m_window.specLength; ++j) {
            localSolarSpectrumData.SetAt(j, m_window.fraunhoferRef.m_data->GetAt(j));
        }

        //----------------------------------------------------------------

        // Copy the measured spectrum to vMeas
        vMeas.Copy(measArray, m_window.specLength, 1);

        // To perform the fit we need to extract the wavelength (or pixel)
        //	information from the vXData-vector
        CVector vXSec(fitHigh - fitLow);
        vXSec.Copy(vXData.SubVector(fitLow, fitHigh - fitLow));

        ////////////////////////////////////////////////////////////////////////////
        // now we start building the model function needed for fitting.
        //
        // First we create a function object that represents our measured spectrum. Since we do not
        // need any interpolation on the measured data its enough to use a CDiscreteFunction object.
        CDiscreteFunction dataTarget;

        // now set the data of the measured spectrum in regard to the wavelength information
        {
            auto temp = vXData.SubVector(measured.m_info.m_startChannel, m_window.specLength);
            dataTarget.SetData(temp, vMeas);
        }

        // since the DOAS model function consists of the sum of all reference spectra and a polynomial,
        // we first create a summation object
        CSimpleDOASFunction cRefSum;

        // reset all reference's parameters
        solarSpec->ResetLinearParameter();
        solarSpec->ResetNonlinearParameter();

        // enable amplitude normalization. This should normally be done in order to avoid numerical
        // problems during fitting.
        solarSpec->SetNormalize(true);

        // set the spectral data of the reference spectrum to the object. This also causes an internal
        // transformation of the spectral data into a B-Spline that will be used to interpolate the 
        // reference spectrum during shift and squeeze operations
        auto tempXVec = vXData.SubVector(0, localSolarSpectrumData.GetSize());
        if (!solarSpec->SetData(tempXVec, localSolarSpectrumData))
        {
            m_lastError = "Failed to evaluate shift: could not initialize spline object for the solar spectrum.";
            free(measArray);
            delete solarSpec;
            return(1);
        }

        // Chech the options for the column value
        if (m_window.fitType == FIT_POLY)
            solarSpec->FixParameter(CReferenceSpectrumFunction::CONCENTRATION, -1.0 * solarSpec->GetAmplitudeScale());
        else
            solarSpec->FixParameter(CReferenceSpectrumFunction::CONCENTRATION, 1.0 * solarSpec->GetAmplitudeScale());

        // Free the shift
        //solarSpec->SetParameterLimits(CReferenceSpectrumFunction::SHIFT,	(TFitData)-6.0, (TFitData)6.0, (TFitData)1e25);
        //solarSpec->FixParameter(CReferenceSpectrumFunction::SHIFT,	(TFitData)1.4);

        // Fix the squeeze
        solarSpec->FixParameter(CReferenceSpectrumFunction::SQUEEZE, (TFitData)1.0);

        // Add the solar-reference to the fit
        cRefSum.AddReference(*solarSpec); // <-- at last add the reference to the summation object

        // Link the shifts of the 'normal' cross sections to the shift of the solar spectrum
        for (i = 0; i < m_window.nRef; ++i) {
            // Link the shift and squeeze to the solar-reference
            solarSpec->LinkParameter(CReferenceSpectrumFunction::SHIFT, *m_ref[i], CReferenceSpectrumFunction::SHIFT);
            solarSpec->LinkParameter(CReferenceSpectrumFunction::SQUEEZE, *m_ref[i], CReferenceSpectrumFunction::SQUEEZE);

            cRefSum.AddReference(*m_ref[i]); // <-- at last add the reference to the summation object
        }

        // create the additional polynomial with the correct order
        //	and add it to the summation object, too
        CPolynomialFunction cPoly(2);
        cRefSum.AddReference(cPoly);

        // the last step in the model function will be to define how the difference between the measured data and the modeled
        // data will be determined. In this case we will use the CStandardMetricFunction which actually just calculate the difference
        // between the measured data and the modeled data channel by channel. The fit will try to minimize these differences.
        // So we create the metric object and set the measured spectrum function object and the DOAS model function object as parameters
        CStandardMetricFunction cDiff(dataTarget, cRefSum);

        /////////////////////////////////////////////////////////////////
        // Now its time to create the fit object. The CStandardFit object will 
        // provide a combination of a linear Least Square Fit and a nonlinear Levenberg-Marquardt Fit, which
        // should be sufficient for most needs.
        CStandardFit cFirstFit(cDiff);

        // don't forget to the the already extracted fit range to the fit object!
        // without a valid fit range you'll get an exception.
        cFirstFit.SetFitRange(vXSec);

        // limit the number of fit iteration to 5000.
        cFirstFit.GetNonlinearMinimizer().SetMaxFitSteps(5000);
        cFirstFit.GetNonlinearMinimizer().SetMinChiSquare(0.0001);

        try
        {
            // prepare everything for fitting
            cFirstFit.PrepareMinimize();

            // actually do the fitting
            if (!cFirstFit.Minimize()) {
                // novac::CString message;
                // message.Format("Fit Failed!");
                // ShowMessage(message);
                m_lastError = "Failed to evaluate shift: fit failed.";
                free(measArray);
                delete solarSpec;
                return 1;
            }

            // finalize the fitting process. This will calculate the error measurements and other statistical stuff
            cFirstFit.FinishMinimize();

            // get the basic fit results
            //long stepNum				= (long)cFirstFit.GetFitSteps();
            // double chiSquare			= (double)cFirstFit.GetChiSquare();
            // unsigned long speciesNum	= (unsigned long)m_window.nRef;

            SaveResidual(cFirstFit);

            // finally get the fit-result
            // double column       = (double)solarSpec->GetModelParameter(CReferenceSpectrumFunction::CONCENTRATION);
            // double columnError  = (double)solarSpec->GetModelParameterError(CReferenceSpectrumFunction::CONCENTRATION);
            shift = (double)solarSpec->GetModelParameter(CReferenceSpectrumFunction::SHIFT);
            shiftError = (double)solarSpec->GetModelParameterError(CReferenceSpectrumFunction::SHIFT);
            squeeze = (double)solarSpec->GetModelParameter(CReferenceSpectrumFunction::SQUEEZE);
            squeezeError = (double)solarSpec->GetModelParameterError(CReferenceSpectrumFunction::SQUEEZE);

            // clean up the evaluation
            free(measArray);
            delete solarSpec;

            return 0;
        }
        catch (CFitException e)
        {
            // in case that something went wrong, display the error to the user.
            // normally you will get error in two cases:
            //
            // 1. You forgot to set a valid fit range before you start fitting
            //
            // 2. A matrix inversion failed for some reason inside the fitting loop. Matrix inversions
            //    normally fail when there are linear dependecies in the matrix respectrively you have linear
            //    dependencies in your reference spectrum. Eg. you tried to fit the same reference spectrum twice at once.

            //  e.ReportError();
            //	std::cout << "Failed: " << ++iFalseCount << std::endl;
            //	std::cout << "Steps: " << cFirstFit.GetNonlinearMinimizer().GetFitSteps() << " - Chi: " << cFirstFit.GetNonlinearMinimizer().GetChiSquare() << std::endl;

            // novac::CString message;
            // message.Format("A Fit Exception has occurred. Are the reference files OK?");
            // ShowMessage(message);

            m_lastError = "Failed to evaluate shift: a fit exception occurred.";

            // clean up the evaluation
            free(measArray);
            delete solarSpec;

            return (1);
        }
    }

}