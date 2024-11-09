#include <SpectralEvaluation/Evaluation/EvaluationBase.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Spectra/Scattering.h>

// include all required fit objects
#include <SpectralEvaluation/Fit/ReferenceSpectrumFunction.h>
#include <SpectralEvaluation/Fit/SimpleDOASFunction.h>
#include <SpectralEvaluation/Fit/StandardMetricFunction.h>
#include <SpectralEvaluation/Fit/StandardFit.h>
#include <SpectralEvaluation/Fit/ExpFunction.h>
#include <SpectralEvaluation/Fit/LnFunction.h>
#include <SpectralEvaluation/Fit/PolynomialFunction.h>
#include <SpectralEvaluation/Fit/NegateFunction.h>
#include <SpectralEvaluation/Fit/MulFunction.h>
#include <SpectralEvaluation/Fit/DivFunction.h>
#include <SpectralEvaluation/Fit/GaussFunction.h>
#include <SpectralEvaluation/Fit/DiscreteFunction.h>
#include <SpectralEvaluation/Fit/DOASVector.h>
#include <SpectralEvaluation/Fit/NonlinearParameterFunction.h>

#include <limits>
#include <sstream>

using namespace MathFit;

namespace novac
{
CEvaluationBase::CEvaluationBase(novac::ILogger& log)
    : m_log(log)
{
    CreateXDataVector(MAX_SPECTRUM_LENGTH);
}

CEvaluationBase::CEvaluationBase(const CFitWindow& window, novac::ILogger& log)
    : m_window(window), m_log(log)
{
    CreateXDataVector(MAX_SPECTRUM_LENGTH);
    CreateReferenceSpectra();
}

CEvaluationBase::~CEvaluationBase()
{
    ClearRefereneSpectra();
}

void CEvaluationBase::SetFitWindow(const CFitWindow& window)
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

CReferenceSpectrumFunction* DefaultReferenceSpectrumFunction();

int CEvaluationBase::CreateReferenceSpectra()
{
    ClearRefereneSpectra();

    // 1) Create the references
    for (size_t i = 0; i < m_window.reference.size(); i++)
    {
        auto newRef = DefaultReferenceSpectrumFunction();

        // set the spectral data of the reference spectrum to the object. This also causes an internal
        // transformation of the spectral data into a B-Spline that will be used to interpolate the 
        // reference spectrum during shift and squeeze operations
        CVector yValues;
        yValues.Copy(m_window.reference[i].m_data->m_crossSection.data(), m_window.reference[i].m_data->GetSize());

        auto tempXVec = vXData.SubVector(0, m_window.reference[i].m_data->GetSize());
        if (!newRef->SetData(tempXVec, yValues))
        {
            Error0("Error initializing spline object!");
            return(1);
        }

        // Finally add this reference to the vector
        this->m_ref.push_back(newRef);
    }

    // 2) Couple the references
    for (size_t i = 0; i < m_window.reference.size(); i++)
    {
        // Chech the options for the column value
        switch (m_window.reference[i].m_columnOption)
        {
        case SHIFT_TYPE::SHIFT_FIX:
            m_ref[i]->FixParameter(CReferenceSpectrumFunction::CONCENTRATION, m_window.reference[i].m_columnValue * m_ref[i]->GetAmplitudeScale());
            break;
        case SHIFT_TYPE::SHIFT_LINK:
            m_ref[(int)m_window.reference[i].m_columnValue]->LinkParameter(CReferenceSpectrumFunction::CONCENTRATION, *m_ref[i], CReferenceSpectrumFunction::CONCENTRATION);
            break;
        case SHIFT_TYPE::SHIFT_FREE:
            m_ref[i]->ReleaseParameter(CReferenceSpectrumFunction::CONCENTRATION);
            break;
        default:
            throw std::invalid_argument("Invalid column option found for reference");
        }

        // Check the options for the shift
        switch (m_window.reference[i].m_shiftOption)
        {
        case SHIFT_TYPE::SHIFT_FIX:
            m_ref[i]->FixParameter(CReferenceSpectrumFunction::SHIFT, m_window.reference[i].m_shiftValue);
            break;
        case SHIFT_TYPE::SHIFT_LINK:
            m_ref[(int)m_window.reference[i].m_shiftValue]->LinkParameter(CReferenceSpectrumFunction::SHIFT, *m_ref[i], CReferenceSpectrumFunction::SHIFT);
            break;
        case SHIFT_TYPE::SHIFT_LIMIT:
            m_ref[i]->SetParameterLimits(CReferenceSpectrumFunction::SHIFT, (TFitData)m_window.reference[i].m_shiftValue, (TFitData)m_window.reference[i].m_shiftMaxValue, 1);
            break;
        default:
            m_ref[i]->SetDefaultParameter(CReferenceSpectrumFunction::SHIFT, (TFitData)0.0);
            m_ref[i]->SetParameterLimits(CReferenceSpectrumFunction::SHIFT, (TFitData)-10.0, (TFitData)10.0, (TFitData)1e0); break;
            // TODO: Get these limits as parameters!
        }

        // Check the options for the squeeze
        switch (m_window.reference[i].m_squeezeOption)
        {
        case SHIFT_TYPE::SHIFT_FIX:
            m_ref[i]->FixParameter(CReferenceSpectrumFunction::SQUEEZE, m_window.reference[i].m_squeezeValue);
            break;
        case SHIFT_TYPE::SHIFT_LINK:
            m_ref[(int)m_window.reference[i].m_squeezeValue]->LinkParameter(CReferenceSpectrumFunction::SQUEEZE, *m_ref[i], CReferenceSpectrumFunction::SQUEEZE);
            break;
        case SHIFT_TYPE::SHIFT_LIMIT:
            m_ref[i]->SetParameterLimits(CReferenceSpectrumFunction::SQUEEZE, (TFitData)m_window.reference[i].m_squeezeValue, (TFitData)m_window.reference[i].m_squeezeMaxValue, 1e7);
            break;
        default:
            m_ref[i]->SetDefaultParameter(CReferenceSpectrumFunction::SQUEEZE, (TFitData)1.0);
            m_ref[i]->SetParameterLimits(CReferenceSpectrumFunction::SQUEEZE, (TFitData)0.98, (TFitData)1.02, (TFitData)1e0);
            break; // TODO: Get these limits as parameters!
        }
    }

    // If we should also include the sky-spectrum in the fit
    if (m_skyReference != nullptr)
    {
        this->m_ref.push_back(m_skyReference);
    }

    // If we should add an additional stray light correction polynomial in intensity space then do this here
    if (m_intensitySpacePolynomial != nullptr)
    {
        this->m_ref.push_back(m_intensitySpacePolynomial);
    }

    // Calculated ring spectra.
    if (m_ringSpectrum != nullptr)
    {
        this->m_ref.push_back(m_ringSpectrum);
    }
    if (m_ringSpectrumLambda4 != nullptr)
    {
        this->m_ref.push_back(m_ringSpectrumLambda4);
    }

    return 0;
}

void CEvaluationBase::RemoveOffset(double* spectrum, int sumChn, bool UV)
{
    // TODO: Read this from the number of optically covered pixels of the SpectrometerModel
    int offsetFrom = (UV) ? 50 : 2;
    int offsetTo = (UV) ? 200 : 20;

    //  remove any remaining offset in the measured spectrum
    double avg = 0;
    for (int i = offsetFrom; i < offsetTo; i++)
    {
        avg += spectrum[i];
    }
    avg = avg / (double)(offsetTo - offsetFrom);

    Sub(spectrum, sumChn, avg);

    return;
}

void CEvaluationBase::RemoveOffset(std::vector<double>& spectrum, size_t from, size_t to)
{
    if (from > to || to > spectrum.size())
    {
        std::stringstream msg;
        msg << "Cannot remove offset, invalid spectrum region [" << from << ", " << to << "] for spectrum of length " << spectrum.size();
        throw std::invalid_argument(msg.str());
    }

    double avg = 0;
    for (size_t i = from; i < to; i++)
    {
        avg += spectrum[i];
    }
    avg = avg / (double)(to - from);

    Sub(spectrum.data(), static_cast<int>(spectrum.size()), avg);

    return;

}

void CEvaluationBase::PrepareSpectra(double* sky, double* meas, const CFitWindow& window)
{

    if (window.fitType == FIT_TYPE::FIT_HP_DIV)
        return PrepareSpectra_HP_Div(sky, meas, window);
    if (window.fitType == FIT_TYPE::FIT_HP_SUB)
        return PrepareSpectra_HP_Sub(sky, meas, window);
    if (window.fitType == FIT_TYPE::FIT_POLY)
        return PrepareSpectra_Poly(sky, meas, window);
}

void CEvaluationBase::PrepareSpectra_HP_Div(double* skyArray, double* measArray, const CFitWindow& window)
{

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

void CEvaluationBase::PrepareSpectra_HP_Sub(double* /*skyArray*/, double* measArray, const CFitWindow& window)
{

    // 1. remove any remaining offset in the measured spectrum
    RemoveOffset(measArray, window.specLength, window.UV);

    // 2. high pass filter
    HighPassBinomial(measArray, window.specLength, 500);

    // 3. log(spec)
    Log(measArray, window.specLength);
}

void CEvaluationBase::PrepareSpectra_Poly(double* /*skyArray*/, double* measArray, const CFitWindow& window)
{
    // 1. remove any remaining offset in the measured spectrum
    RemoveOffset(measArray, window.specLength, window.UV);

    // 2. log(spec)
    Log(measArray, window.specLength);

    // 3. Multiply the spectrum with -1 to get the correct sign for everything
    for (int i = 0; i < window.specLength; ++i)
    {
        measArray[i] *= -1.0;
    }
}

int CEvaluationBase::SetSkySpectrum(const CSpectrum& spec)
{
    CCrossSectionData sky;
    sky.m_crossSection = std::vector<double>(spec.m_data, spec.m_data + spec.m_length);
    if (spec.m_wavelength.size() > 0)
    {
        sky.m_waveLength = spec.m_wavelength;
    }
    return SetSkySpectrum(sky, true);
}

int CEvaluationBase::SetSkySpectrum(const CCrossSectionData& spec, bool removeOffset)
{
    // make sure that these make sense
    if (m_window.specLength != (int)spec.GetSize())
    {
        return 1;
    }

    // set the data
    m_sky = spec;

    // --------- prepare the sky spectrum for evaluation ---------
    // Remove any remaining offset of the sky-spectrum
    if (removeOffset)
    {
        RemoveOffset(m_sky.m_crossSection.data(), m_sky.GetSize(), m_window.UV);
    }

    if (m_window.fitType == FIT_TYPE::FIT_HP_SUB)
    {
        HighPassBinomial(m_sky.m_crossSection.data(), m_sky.GetSize(), 500);
    }

    if (m_window.fitType == FIT_TYPE::FIT_POLY && m_window.includeIntensitySpacePolyominal)
    {
        CreateReferenceForIntensitySpacePolynomial(m_sky.m_crossSection);
    }

    if (m_window.ringCalculation == RING_CALCULATION_OPTION::CALCULATE_RING ||
        m_window.ringCalculation == RING_CALCULATION_OPTION::CALCULATE_RING_X2)
    {
        CSpectrum skySpectrum;
        if (m_sky.m_waveLength.size() == m_sky.m_crossSection.size())
        {
            skySpectrum = CSpectrum(m_sky.m_waveLength, m_sky.m_crossSection);
        }
        else if (m_window.fraunhoferRef.m_data != nullptr)
        {
            skySpectrum = CSpectrum(m_window.fraunhoferRef.m_data->m_waveLength, m_sky.m_crossSection);
        }
        else
        {
            throw std::invalid_argument("Cannot calculate a ring spectrum unless the evaluation has a valid wavelength calibration.");
        }
        auto ringSpectrum = Doasis::Scattering::CalcRingSpectrum(skySpectrum);
        CreateReferenceForRingSpectrum(ringSpectrum);

        if (m_window.ringCalculation == RING_CALCULATION_OPTION::CALCULATE_RING_X2)
        {
            CreateReferenceForRingSpectrumLambda4(ringSpectrum);
        }
    }

    if (m_window.fitType != FIT_TYPE::FIT_HP_DIV)
    {
        // Take the log of the sky-spectrum
        Log(m_sky.m_crossSection.data(), m_sky.GetSize());

        // Include the sky-spectrum in the fit
        CreateReferenceForSkySpectrum();
        CreateReferenceSpectra();
    }

    return 0;
}

void CEvaluationBase::CreateReferenceForIntensitySpacePolynomial(const std::vector<double>& I0)
{
    if (m_intensitySpacePolynomial != nullptr)
    {
        delete m_intensitySpacePolynomial;
    }

    m_intensitySpacePolynomial = DefaultReferenceSpectrumFunction();

    // set the spectral data of the reference spectrum to the object. This also causes an internal
    // transformation of the spectral data into a B-Spline that will be used to interpolate the
    // reference spectrum during shift and squeeze operations
    CVector yValues;
    yValues.SetSize((int)I0.size());
    for (int k = 0; k < (int)I0.size(); ++k)
    {
        yValues.SetAt(k, std::abs(I0[k]) < std::numeric_limits<double>::epsilon() ? 0.0 : 1.0 / I0[k]);
    }

    auto tempXVec = vXData.SubVector(0, (int)I0.size());
    if (!m_intensitySpacePolynomial->SetData(tempXVec, yValues))
    {
        Error0("Error initializing spline object!");
        return;
    }

    // Don't shift / squeeze this reference...
    m_intensitySpacePolynomial->FixParameter(CReferenceSpectrumFunction::SHIFT, 0.0);
    m_intensitySpacePolynomial->FixParameter(CReferenceSpectrumFunction::SQUEEZE, 1.0);
}

void CEvaluationBase::CreateReferenceForSkySpectrum()
{
    if (m_skyReference != nullptr)
    {
        delete m_skyReference;
    }

    m_skyReference = DefaultReferenceSpectrumFunction();

    // set the spectral data of the reference spectrum to the object. This also causes an internal
    // transformation of the spectral data into a B-Spline that will be used to interpolate the
    // reference spectrum during shift and squeeze operations
    CVector yValues;
    yValues.Copy(m_sky.m_crossSection.data(), m_sky.GetSize());

    auto tempXVec = vXData.SubVector(0, m_sky.GetSize());
    if (!m_skyReference->SetData(tempXVec, yValues))
    {
        Error0("Error initializing spline object!");
        return;
    }

    // Check the options for the column value
    const double concentrationMultiplier = (m_window.fitType == FIT_TYPE::FIT_POLY) ? -1.0 : 1.0;
    m_skyReference->FixParameter(CReferenceSpectrumFunction::CONCENTRATION, concentrationMultiplier * m_skyReference->GetAmplitudeScale());

    // Check the options for the shift & squeeze
    if (m_window.shiftSky == 1)
    {
        m_skyReference->SetParameterLimits(CReferenceSpectrumFunction::SHIFT, (TFitData)-3.0, (TFitData)3.0, 1);
        m_skyReference->SetParameterLimits(CReferenceSpectrumFunction::SQUEEZE, (TFitData)0.95, (TFitData)1.05, 1e7);
    }
    else if (m_window.shiftSky == 2)
    {
        m_skyReference->FixParameter(CReferenceSpectrumFunction::SHIFT, (TFitData)m_window.skyShift);
        m_skyReference->FixParameter(CReferenceSpectrumFunction::SQUEEZE, (TFitData)m_window.skySqueeze);
    }
    else
    {
        m_skyReference->FixParameter(CReferenceSpectrumFunction::SHIFT, 0.0);
        m_skyReference->FixParameter(CReferenceSpectrumFunction::SQUEEZE, 1.0);
    }

    if (m_ringSpectrum != nullptr)
    {
        m_ringSpectrum->LinkParameter(CReferenceSpectrumFunction::SHIFT, *m_skyReference, CReferenceSpectrumFunction::SHIFT);
        m_ringSpectrum->LinkParameter(CReferenceSpectrumFunction::SQUEEZE, *m_skyReference, CReferenceSpectrumFunction::SQUEEZE);
    }

    if (m_ringSpectrumLambda4 != nullptr)
    {
        m_ringSpectrumLambda4->LinkParameter(CReferenceSpectrumFunction::SHIFT, *m_skyReference, CReferenceSpectrumFunction::SHIFT);
        m_ringSpectrumLambda4->LinkParameter(CReferenceSpectrumFunction::SQUEEZE, *m_skyReference, CReferenceSpectrumFunction::SQUEEZE);

        // m_ringSpectrumLambda4->FixParameter(CReferenceSpectrumFunction::SHIFT, 0.0);
        // m_ringSpectrumLambda4->FixParameter(CReferenceSpectrumFunction::SQUEEZE, 1.0);
    }
}

void CEvaluationBase::CreateReferenceForRingSpectrum(const CSpectrum& ring)
{
    if (m_ringSpectrum != nullptr)
    {
        delete m_ringSpectrum;
    }

    m_ringSpectrum = DefaultReferenceSpectrumFunction();

    // set the spectral data of the reference spectrum to the object. This also causes an internal
    // transformation of the spectral data into a B-Spline that will be used to interpolate the
    // reference spectrum during shift and squeeze operations
    CVector yValues;
    yValues.Copy((double*)ring.m_data, (int)ring.m_length);

    auto tempXVec = vXData.SubVector(0, (int)ring.m_length);
    if (!m_ringSpectrum->SetData(tempXVec, yValues))
    {
        Error0("Error initializing spline object!");
        return;
    }

    // m_ringSpectrum->FixParameter(CReferenceSpectrumFunction::SHIFT, 0.0);
    // m_ringSpectrum->FixParameter(CReferenceSpectrumFunction::SQUEEZE, 1.0);
}

void CEvaluationBase::CreateReferenceForRingSpectrumLambda4(const CSpectrum& ring)
{
    if (m_ringSpectrumLambda4 != nullptr)
    {
        delete m_ringSpectrumLambda4;
    }

    m_ringSpectrumLambda4 = DefaultReferenceSpectrumFunction();

    // set the spectral data of the reference spectrum to the object. This also causes an internal
    // transformation of the spectral data into a B-Spline that will be used to interpolate the
    // reference spectrum during shift and squeeze operations
    CVector yValues;
    yValues.SetSize((int)ring.m_length);
    for (size_t ii = 0; ii < (size_t)ring.m_length; ++ii)
    {
        double lambda = ring.m_wavelength[ii];
        yValues.SetAt((int)ii, ring.m_data[ii] * std::pow(lambda, 4.0));
    }

    auto tempXVec = vXData.SubVector(0, (int)ring.m_length);
    if (!m_ringSpectrumLambda4->SetData(tempXVec, yValues))
    {
        Error0("Error initializing spline object!");
        return;
    }

    m_ringSpectrumLambda4->FixParameter(CReferenceSpectrumFunction::SHIFT, 0.0);
    m_ringSpectrumLambda4->FixParameter(CReferenceSpectrumFunction::SQUEEZE, 1.0);
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

int CEvaluationBase::Evaluate(const CSpectrum& measured, int numSteps)
{
    return Evaluate(measured.m_data, measured.m_length, measured.m_info.m_startChannel, numSteps);
}

int CEvaluationBase::Evaluate(const double* measured, size_t measuredLength, int measuredStartChannel, int numSteps)
{
    assert(static_cast<size_t>(vXData.GetSize()) >= measuredLength);

    m_lastError = "";

    // Check so that the length of the spectra agree with each other
    if (static_cast<size_t>(m_window.specLength) != measuredLength)
    {
        m_lastError = "Failed to evaluate: the length of the measured spectrum does not equal the spectrum length of the fit-window used.";
        return 1;
    }

    // Get the correct limits for the fit
    int fitLow, fitHigh;
    if (measuredStartChannel != 0 || measuredLength < static_cast<size_t>(m_window.specLength))
    {
        // Partial spectra
        fitLow = m_window.fitLow - measuredStartChannel;
        fitHigh = m_window.fitHigh - measuredStartChannel;
        if (fitLow < 0 || fitHigh > static_cast<int>(measuredLength))
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

    // Make a local copy of the data (since we're going to change the contents)
    std::vector<double> measArray(measured, measured + measuredLength);

    //----------------------------------------------------------------
    // --------- prepare the spectrum for evaluation -----------------
    //----------------------------------------------------------------

    PrepareSpectra(m_sky.m_crossSection.data(), measArray.data(), m_window);
    m_measuredData = std::vector<double>(begin(measArray), end(measArray));

    //----------------------------------------------------------------

    // Copy the measured spectrum to vMeas
    CVector vMeas;
    vMeas.Copy(measArray.data(), m_window.specLength, 1);

    // To perform the fit we need to extract the wavelength (or pixel)
    //  information from the vXData-vector
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
        auto temp = vXData.SubVector(measuredStartChannel, m_window.specLength);
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
        m_result.m_stepNum = (long)cFirstFit.GetFitSteps();
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
            m_result.m_referenceResult[ii].m_specieName = GetReferenceName(ii);
            m_result.m_referenceResult[ii].m_column = (double)m_ref[ii]->GetModelParameter(CReferenceSpectrumFunction::CONCENTRATION);
            m_result.m_referenceResult[ii].m_columnError = (double)m_ref[ii]->GetModelParameterError(CReferenceSpectrumFunction::CONCENTRATION);
            m_result.m_referenceResult[ii].m_shift = (double)m_ref[ii]->GetModelParameter(CReferenceSpectrumFunction::SHIFT);
            m_result.m_referenceResult[ii].m_shiftError = (double)m_ref[ii]->GetModelParameterError(CReferenceSpectrumFunction::SHIFT);
            m_result.m_referenceResult[ii].m_squeeze = (double)m_ref[ii]->GetModelParameter(CReferenceSpectrumFunction::SQUEEZE);
            m_result.m_referenceResult[ii].m_squeezeError = (double)m_ref[ii]->GetModelParameterError(CReferenceSpectrumFunction::SQUEEZE);

            // get the final fit result
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
            m_lastError = m_lastError + " Message: '" + std::string{ e.mMessage } + "'";
        }

        return (1);
    }
}

int CEvaluationBase::EvaluateShift(novac::LogContext context, const CSpectrum& measured, ShiftEvaluationResult& shiftResult)
{
    assert(vXData.GetSize() >= measured.m_length);

    m_lastError = "";

    int fitLow = m_window.fitLow;
    int fitHigh = m_window.fitHigh;

    // Check so that the length of the spectra agree with each other
    if (m_window.specLength != measured.m_length)
    {
        m_log.Error(context, "Failed to evaluate shift: the length of the measured spectrum does not equal the spectrum length of the fit-window used.");
        return 1;
    }

    // Check that we have a solar-spectrum to check against
    if (m_window.fraunhoferRef.m_path.size() < 6)
    {
        m_log.Error(context, "Failed to evaluate shift: the fraunhofer reference does not exist.");
        return 1;
    }

    if (measured.m_info.m_startChannel != 0 || measured.m_length < m_window.specLength)
    {
        // Partial spectra
        fitLow = m_window.fitLow - measured.m_info.m_startChannel;
        fitHigh = m_window.fitHigh - measured.m_info.m_startChannel;
        if (fitLow < 0 || fitHigh > measured.m_length)
        {
            m_log.Error(context, "Invalid fit region, fit-low is negative or fit-high is larger than the length of the measured spectrum.");
            return 1;
        }
    }
    else
    {
        fitLow = m_window.fitLow;
        fitHigh = m_window.fitHigh;
    }

    // initialize the solar-spectrum function
    std::unique_ptr<CReferenceSpectrumFunction> solarSpec = std::make_unique<CReferenceSpectrumFunction>();

    // Make a local copy of the data
    std::vector<double> measArray(measured.m_length, 0.0);
    memcpy(measArray.data(), measured.m_data, measured.m_length * sizeof(double));

    //----------------------------------------------------------------
    // --------- prepare the spectrum for evaluation -----------------
    //----------------------------------------------------------------

    RemoveOffset(measArray.data(), m_window.specLength, m_window.UV);
    if (m_window.fitType == FIT_TYPE::FIT_HP_DIV || m_window.fitType == FIT_TYPE::FIT_HP_SUB)
    {
        HighPassBinomial(measArray.data(), m_window.specLength, 500);
    }
    Log(measArray.data(), m_window.specLength);

    if (m_window.fitType == FIT_TYPE::FIT_POLY)
    {
        for (int j = 0; j < m_window.specLength; ++j)
        {
            measArray[j] *= -1.0;
        }
    }

    // --------- also prepare the solar-spectrum for evaluation -----------------
    CVector localSolarSpectrumData;
    localSolarSpectrumData.SetSize(m_window.fraunhoferRef.m_data->GetSize());
    for (int j = 0; j < m_window.specLength; ++j)
    {
        localSolarSpectrumData.SetAt(j, m_window.fraunhoferRef.m_data->GetAt(j));
    }

    //----------------------------------------------------------------

    // Copy the measured spectrum to vMeas
    CVector vMeas;
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
        m_log.Error(context, "Failed to evaluate shift: could not initialize spline object for the solar spectrum.");
        return(1);
    }

    // Chech the options for the column value
    const double concentrationMultiplier = (m_window.fitType == FIT_TYPE::FIT_POLY) ? -1.0 : 1.0;
    solarSpec->FixParameter(CReferenceSpectrumFunction::CONCENTRATION, concentrationMultiplier * solarSpec->GetAmplitudeScale());
    solarSpec->ReleaseParameter(CReferenceSpectrumFunction::SHIFT);

    // Fix the squeeze
    solarSpec->FixParameter(CReferenceSpectrumFunction::SQUEEZE, (TFitData)1.0);

    // Add the solar-reference to the fit
    cRefSum.AddReference(*solarSpec); // <-- at last add the reference to the summation object

    // Link the shifts of the 'normal' cross sections to the shift of the solar spectrum
    for (size_t ii = 0; ii < m_window.reference.size(); ++ii)
    {
        // Link the shift and squeeze to the solar-reference
        m_ref[ii]->ReleaseParameter(CReferenceSpectrumFunction::SHIFT);
        m_ref[ii]->ReleaseParameter(CReferenceSpectrumFunction::SQUEEZE);
        solarSpec->LinkParameter(CReferenceSpectrumFunction::SHIFT, *m_ref[ii], CReferenceSpectrumFunction::SHIFT);
        solarSpec->LinkParameter(CReferenceSpectrumFunction::SQUEEZE, *m_ref[ii], CReferenceSpectrumFunction::SQUEEZE);

        cRefSum.AddReference(*m_ref[ii]); // <-- at last add the reference to the summation object
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
        if (!cFirstFit.Minimize())
        {
            m_log.Error(context, "Failed to evaluate shift: fit failed.");
            return 1;
        }

        // finalize the fitting process. This will calculate the error measurements and other statistical stuff
        cFirstFit.FinishMinimize();

        // get the basic fit results
        // long stepNum = (long)cFirstFit.GetFitSteps();
        shiftResult.chi2 = (double)cFirstFit.GetChiSquare();
        // unsigned long speciesNum = (unsigned long)m_window.nRef;

        SaveResidual(cFirstFit);

        // finally get the fit-result
        assert(std::abs((double)solarSpec->GetModelParameter(CReferenceSpectrumFunction::CONCENTRATION) - concentrationMultiplier) < 1e-6);
        assert(std::abs((double)solarSpec->GetModelParameterError(CReferenceSpectrumFunction::CONCENTRATION)) < 1e-6);

        shiftResult.shift = (double)solarSpec->GetModelParameter(CReferenceSpectrumFunction::SHIFT);
        shiftResult.shiftError = (double)solarSpec->GetModelParameterError(CReferenceSpectrumFunction::SHIFT);
        shiftResult.squeeze = (double)solarSpec->GetModelParameter(CReferenceSpectrumFunction::SQUEEZE);
        shiftResult.squeezeError = (double)solarSpec->GetModelParameterError(CReferenceSpectrumFunction::SQUEEZE);

        // also verify that the setup did what we expected out of it...
        for (size_t ii = 0; ii < m_window.reference.size(); ++ii)
        {
            assert(std::abs(shiftResult.shift - m_ref[ii]->GetModelParameter(CReferenceSpectrumFunction::SHIFT)) < 1e-3);
            assert(std::abs(shiftResult.squeeze - m_ref[ii]->GetModelParameter(CReferenceSpectrumFunction::SQUEEZE)) < 1e-3);
        }

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

        m_log.Error(context, "Failed to evaluate shift: a fit exception occurred.");

        return (1);
    }
}

std::string CEvaluationBase::GetReferenceName(size_t referenceIndex) const
{
    if (referenceIndex < m_window.reference.size())
    {
        return m_window.reference[referenceIndex].m_specieName; // user supplied reference
    }
    else if (referenceIndex >= this->m_ref.size())
    {
        return "N/A";
    }
    else if (m_ref[referenceIndex] == this->m_skyReference)
    {
        return "Sky";
    }
    else if (m_ref[referenceIndex] == this->m_intensitySpacePolynomial)
    {
        return "IntensitySpacePolynomial";
    }
    else if (m_ref[referenceIndex] == this->m_ringSpectrum)
    {
        return "Ring";
    }
    else if (m_ref[referenceIndex] == this->m_ringSpectrumLambda4)
    {
        return "Ring * lambda^4";
    }

    return "N/A"; // unknown.
}
}