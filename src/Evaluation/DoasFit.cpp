#include <SpectralEvaluation/Evaluation/DoasFit.h>
#include <SpectralEvaluation/Evaluation/FitWindow.h>

#include <SpectralEvaluation/Fit/Vector.h>
#include <SpectralEvaluation/Fit/DiscreteFunction.h>
#include <SpectralEvaluation/Fit/PolynomialFunction.h>
#include <SpectralEvaluation/Fit/ReferenceSpectrumFunction.h>
#include <SpectralEvaluation/Fit/SimpleDOASFunction.h>
#include <SpectralEvaluation/Fit/StandardFit.h>
#include <SpectralEvaluation/Fit/StandardMetricFunction.h>

#include <sstream>

namespace novac
{

// ----------------------------- Free functions, to help with the fit below -----------------------------
MathFit::CVector Generate(int first, int last)
{
    assert(last > first);
    MathFit::CVector result(last - first);
    for (int i = 0; i < last - first; ++i)
    {
        result.SetAt(i, static_cast<MathFit::TFitData>(first + i));
    }

    return result;
}

void SaveResidual(MathFit::CStandardFit& cFirstFit, DoasResult& result)
{
    const auto& res = cFirstFit.GetResiduum();
    result.residual.resize(res.GetSize());
    for (int ii = 0; ii < res.GetSize(); ++ii)
    {
        result.residual[ii] = res.GetAt(ii);
    }
}

void SavePolynomial(MathFit::CPolynomialFunction& fittedPolynomial, int fitLow, int fitHigh, DoasResult& result)
{
    result.polynomialValues.resize(fitHigh - fitLow);
    for (int ii = 0; ii < fitHigh - fitLow; ++ii)
    {
        result.polynomialValues[ii] = fittedPolynomial.GetValue(static_cast<MathFit::TFitData>(ii + fitLow));
    }
}

MathFit::CReferenceSpectrumFunction* DefaultReferenceSpectrumFunction()
{
    auto newRef = new MathFit::CReferenceSpectrumFunction();

    // reset all reference's parameters
    newRef->ResetLinearParameter();
    newRef->ResetNonlinearParameter();

    // enable amplitude normalization. This should normally be done in order to avoid numerical
    // problems during fitting.
    newRef->SetNormalize(true);

    return newRef;
}

// Helper class, for storing the references.
class DoasReferenceSetup
{
public:
    // The CReferenceSpectrumFunctions are used in the evaluation process to model
    // the reference spectra for the different species that are being fitted.
    //  The vector must hold pointers to the references, as these cannot be copied...
    // TODO: Review the structure here.
    std::vector<MathFit::CReferenceSpectrumFunction*> m_ref;

    /// The name of each reference.
    std::vector<std::string> name;

    // The scale factor for the retrieved values, can be either +1.0 or -1.0 depending on the fit-type.
    double columnScaleFactor = -1.0;
};

DoasFit::DoasFit()
{
}

DoasFit::~DoasFit()
{
    DeallocateReferenceSetup();
}

void DoasFit::DeallocateReferenceSetup()
{
    DoasReferenceSetup* setup = static_cast<DoasReferenceSetup*>(m_referenceSetup);

    if (setup != nullptr)
    {
        for (MathFit::CReferenceSpectrumFunction* reference : setup->m_ref)
        {
            delete reference;
            reference = nullptr;
        }
        delete setup;
    }

    m_referenceSetup = nullptr;
}

void DoasFit::Setup(const CFitWindow& setup)
{
    m_fitLow = setup.fitLow;
    m_fitHigh = setup.fitHigh;
    m_polynomialOrder = setup.polyOrder;
    m_name = setup.name;

    DeallocateReferenceSetup();
    DoasReferenceSetup* newReferenceSetup = new DoasReferenceSetup();
    newReferenceSetup->columnScaleFactor = (setup.fitType == FIT_TYPE::FIT_POLY) ? -1.0 : +1.0;

    // 1) Create the references
    for (int refIdx = 0; refIdx < setup.nRef; ++refIdx)
    {
        if (setup.ref[refIdx].m_data == nullptr)
        {
            throw std::invalid_argument("Error in setting up DOAS fit, reference is null.");
        }
        if (setup.ref[refIdx].m_data->GetSize() == 0)
        {
            throw std::invalid_argument("Error in setting up DOAS fit, reference does not contain any data.");
        }

        auto newRef = DefaultReferenceSpectrumFunction();

        // set the spectral data of the reference spectrum to the object. This also causes an internal
        // transformation of the spectral data into a B-Spline that will be used to interpolate the 
        // reference spectrum during shift and squeeze operations
        MathFit::CVector yValues;
        yValues.Copy(setup.ref[refIdx].m_data->m_crossSection.data(), setup.ref[refIdx].m_data->GetSize());

        auto tempXVec = Generate(0, setup.ref[refIdx].m_data->GetSize()); // the x-axis vector here is pixels.
        if (!newRef->SetData(tempXVec, yValues))
        {
            throw std::invalid_argument("Error in DOAS reference, failed to initialize spline object. Make sure that the reference is ok and try again.");
        }

        // Finally add this reference to the vector
        newReferenceSetup->m_ref.push_back(newRef);
        newReferenceSetup->name.push_back(setup.ref[refIdx].m_specieName);
    }

    // 2) Couple the references
    for (int refIdx = 0; refIdx < setup.nRef; ++refIdx)
    {
        // Check the options for the column value.
        //  Notice the multiplication with minus one here, this is done to keep the signs of everything compatible with DOASIS.
        switch (setup.ref[refIdx].m_columnOption)
        {
        case novac::SHIFT_TYPE::SHIFT_FIX:
            newReferenceSetup->m_ref[refIdx]->FixParameter(MathFit::CReferenceSpectrumFunction::CONCENTRATION, -1.0 * setup.ref[refIdx].m_columnValue * newReferenceSetup->m_ref[refIdx]->GetAmplitudeScale());
            break;
        case novac::SHIFT_TYPE::SHIFT_LINK:
            newReferenceSetup->m_ref[(int)setup.ref[refIdx].m_columnValue]->LinkParameter(MathFit::CReferenceSpectrumFunction::CONCENTRATION, *newReferenceSetup->m_ref[refIdx], MathFit::CReferenceSpectrumFunction::CONCENTRATION);
            break;
        case novac::SHIFT_TYPE::SHIFT_FREE:
            newReferenceSetup->m_ref[refIdx]->ReleaseParameter(MathFit::CReferenceSpectrumFunction::CONCENTRATION);
            break;
        default:
            throw std::invalid_argument("Invalid type of shift set for the column option for a reference in DoasFit");
        }

        // Check the options for the shift
        switch (setup.ref[refIdx].m_shiftOption)
        {
        case novac::SHIFT_TYPE::SHIFT_FIX:
            newReferenceSetup->m_ref[refIdx]->FixParameter(MathFit::CReferenceSpectrumFunction::SHIFT, setup.ref[refIdx].m_shiftValue);
            break;
        case novac::SHIFT_TYPE::SHIFT_LINK:
            newReferenceSetup->m_ref[(int)setup.ref[refIdx].m_shiftValue]->LinkParameter(MathFit::CReferenceSpectrumFunction::SHIFT, *newReferenceSetup->m_ref[refIdx], MathFit::CReferenceSpectrumFunction::SHIFT);
            break;
        case novac::SHIFT_TYPE::SHIFT_LIMIT:
            newReferenceSetup->m_ref[refIdx]->SetParameterLimits(MathFit::CReferenceSpectrumFunction::SHIFT, (MathFit::TFitData)setup.ref[refIdx].m_shiftValue, (MathFit::TFitData)setup.ref[refIdx].m_shiftMaxValue, 1);
            break;
        default:
            newReferenceSetup->m_ref[refIdx]->SetDefaultParameter(MathFit::CReferenceSpectrumFunction::SHIFT, (MathFit::TFitData)0.0);
            newReferenceSetup->m_ref[refIdx]->SetParameterLimits(MathFit::CReferenceSpectrumFunction::SHIFT, (MathFit::TFitData)-10.0, (MathFit::TFitData)10.0, (MathFit::TFitData)1e0);
            break; // TODO: Get these limits as parameters!
        }

        // Check the options for the squeeze
        switch (setup.ref[refIdx].m_squeezeOption)
        {
        case novac::SHIFT_TYPE::SHIFT_FIX:
            newReferenceSetup->m_ref[refIdx]->FixParameter(MathFit::CReferenceSpectrumFunction::SQUEEZE, setup.ref[refIdx].m_squeezeValue);
            break;
        case novac::SHIFT_TYPE::SHIFT_LINK:
            newReferenceSetup->m_ref[(int)setup.ref[refIdx].m_squeezeValue]->LinkParameter(MathFit::CReferenceSpectrumFunction::SQUEEZE, *newReferenceSetup->m_ref[refIdx], MathFit::CReferenceSpectrumFunction::SQUEEZE);
            break;
        case novac::SHIFT_TYPE::SHIFT_LIMIT:
            newReferenceSetup->m_ref[refIdx]->SetParameterLimits(MathFit::CReferenceSpectrumFunction::SQUEEZE, (MathFit::TFitData)setup.ref[refIdx].m_squeezeValue, (MathFit::TFitData)setup.ref[refIdx].m_squeezeMaxValue, 1e7);
            break;
        default:
            newReferenceSetup->m_ref[refIdx]->SetDefaultParameter(MathFit::CReferenceSpectrumFunction::SQUEEZE, (MathFit::TFitData)1.0);
            newReferenceSetup->m_ref[refIdx]->SetParameterLimits(MathFit::CReferenceSpectrumFunction::SQUEEZE, (MathFit::TFitData)0.98, (MathFit::TFitData)1.02, (MathFit::TFitData)1e0);
            break; // TODO: Get these limits as parameters!
        }
    }

    assert(newReferenceSetup->name.size() == newReferenceSetup->m_ref.size());

    // Set the member
    if (m_referenceSetup != nullptr)
    {
        auto currentSetup = static_cast<DoasReferenceSetup*>(m_referenceSetup);
        delete currentSetup;
    }
    m_referenceSetup = newReferenceSetup;
}

void ValidateDoasInputData(const double* measuredData, size_t measuredLength, const DoasReferenceSetup* referenceSetup)
{
    if (measuredData == nullptr)
    {
        throw std::invalid_argument("Cannot perform a DOAS fit on a null measured spectrum.");
    }
    if (measuredLength == 0)
    {
        throw std::invalid_argument("Cannot perform a DOAS fit if the measured and reference spectra have zero length.");
    }
    if (static_cast<int>(measuredLength) < 0)
    {
        throw std::invalid_argument("Invalid input to DOAS fit, overflow error in measured spectrum length.");
    }
    if (referenceSetup == nullptr)
    {
        throw std::invalid_argument("Invalid setup of DOAS fit, the references must be setup first.");
    }
    if (referenceSetup->m_ref.size() == 0)
    {
        throw std::invalid_argument("Invalid setup of DOAS fit, at least one reference must be setup before doing the DOAS fit.");
    }
}

void DoasFit::Run(const double* measuredData, size_t measuredLength, DoasResult& result)
{
    DoasReferenceSetup* referenceSetup = static_cast<DoasReferenceSetup*>(m_referenceSetup);

    ValidateDoasInputData(measuredData, measuredLength, referenceSetup);

    // Make a local copy of the data. TODO: Check if this actually is necessary anymore?!?
    std::vector<double> measArray(measuredData, measuredData + measuredLength);

    //----------------------------------------------------------------

    // Copy the measured spectrum to vMeas
    MathFit::CVector vMeas;
    vMeas.Copy(measArray.data(), static_cast<int>(measuredLength), 1);

    // To perform the fit we need to extract the wavelength (or pixel)
    //  information from the vXData-vector
    MathFit::CVector vXSec = Generate(m_fitLow, m_fitHigh); // the x-axis data of the fit, here in pixels

    ////////////////////////////////////////////////////////////////////////////
    // now we start building the model function needed for fitting.
    //
    // First we create a function object that represents our measured spectrum. Since we do not
    // need any interpolation on the measured data its enough to use a CDiscreteFunction object.
    MathFit::CDiscreteFunction dataTarget;

    // now set the data of the measured spectrum in regard to the wavelength information
    {
        // use channel base fitting.
        auto temp = Generate(0, static_cast<int>(measuredLength));
        dataTarget.SetData(temp, vMeas);
    }

    // since the DOAS model function consists of the sum of all reference spectra and a polynomial,
    // we first create a summation object
    MathFit::CSimpleDOASFunction cRefSum;

    // now we add the required CReferenceSpectrumFunction objects that actually represent the 
    // reference spectra used in the DOAS model function
    for (MathFit::CReferenceSpectrumFunction* reference : referenceSetup->m_ref)
    {
        cRefSum.AddReference(*reference); // <-- at last add the reference to the summation object
    }

    // create the additional polynomial with the correct order
    //	and add it to the summation object, too
    MathFit::CPolynomialFunction cPoly(m_polynomialOrder);
    cRefSum.AddReference(cPoly);

    // the last step in the model function will be to define how the difference between the measured data and the modeled
    // data will be determined. In this case we will use the CStandardMetricFunction which actually just calculate the difference
    // between the measured data and the modeled data channel by channel. The fit will try to minimize these differences.
    // So we create the metric object and set the measured spectrum function object and the DOAS model function object as parameters
    MathFit::CStandardMetricFunction cDiff(dataTarget, cRefSum);

    /////////////////////////////////////////////////////////////////
    // Now its time to create the fit object. The CStandardFit object will 
    // provide a combination of a linear Least Square Fit and a nonlinear Levenberg-Marquardt Fit, which
    // should be sufficient for most needs.
    MathFit::CStandardFit cFirstFit(cDiff);

    // don't forget to the the already extracted fit range to the fit object!
    // without a valid fit range you'll get an exception.
    cFirstFit.SetFitRange(vXSec);

    // limit the number of fit iteration to 5000. This can still take a long time! More convinient values are
    // between 100 and 1000
    cFirstFit.GetNonlinearMinimizer().SetMaxFitSteps(m_maximumNumberOfSteps);
    cFirstFit.GetNonlinearMinimizer().SetMinChiSquare(0.0001);

    try
    {
        // prepare everything for fitting
        cFirstFit.PrepareMinimize();

        // actually do the fitting
        if (!cFirstFit.Minimize())
        {
            throw DoasFitException("Failed to evaluate: fit " + m_name + " failed.");
        }

        // finalize the fitting process. This will calculate the error measurements and other statistical stuff
        cFirstFit.FinishMinimize();

        // Save the results of the fit.
        result.fitLow = m_fitLow;
        result.fitHigh = m_fitHigh;
        result.iterations = (long)cFirstFit.GetFitSteps();
        result.chiSquare = (double)cFirstFit.GetChiSquare();
        result.delta = cFirstFit.GetResiduum().Max() - cFirstFit.GetResiduum().Min();

        result.polynomialCoefficients.resize(1 + m_polynomialOrder);
        for (int tmpInt = 0; tmpInt <= m_polynomialOrder; ++tmpInt)
        {
            result.polynomialCoefficients[tmpInt] = (double)cPoly.GetCoefficient(tmpInt);
        }

        SaveResidual(cFirstFit, result);

        SavePolynomial(cPoly, m_fitLow, m_fitHigh, result);

        // Save the filtered measured spectrum
        result.measuredSpectrum = std::vector<double>(begin(measArray) + m_fitLow, begin(measArray) + m_fitHigh);

        // finally display the fit results for each reference spectrum including their appropriate error
        result.referenceResult.resize(referenceSetup->m_ref.size());
        for (size_t ii = 0; ii < referenceSetup->m_ref.size(); ii++)
        {
            result.referenceResult[ii].name = referenceSetup->name[ii];
            result.referenceResult[ii].column = referenceSetup->columnScaleFactor * (double)referenceSetup->m_ref[ii]->GetModelParameter(MathFit::CReferenceSpectrumFunction::CONCENTRATION);
            result.referenceResult[ii].columnError = (double)referenceSetup->m_ref[ii]->GetModelParameterError(MathFit::CReferenceSpectrumFunction::CONCENTRATION);
            result.referenceResult[ii].shift = (double)referenceSetup->m_ref[ii]->GetModelParameter(MathFit::CReferenceSpectrumFunction::SHIFT);
            result.referenceResult[ii].shiftError = (double)referenceSetup->m_ref[ii]->GetModelParameterError(MathFit::CReferenceSpectrumFunction::SHIFT);
            result.referenceResult[ii].squeeze = (double)referenceSetup->m_ref[ii]->GetModelParameter(MathFit::CReferenceSpectrumFunction::SQUEEZE);
            result.referenceResult[ii].squeezeError = (double)referenceSetup->m_ref[ii]->GetModelParameterError(MathFit::CReferenceSpectrumFunction::SQUEEZE);

            // Get the scaled reference as well
            result.referenceResult[ii].scaledValues.reserve(m_fitHigh - m_fitLow);
            for (int pixelIdx = m_fitLow; pixelIdx < m_fitHigh; ++pixelIdx)
            {
                result.referenceResult[ii].scaledValues.push_back(referenceSetup->m_ref[ii]->GetValue(static_cast<MathFit::TFitData>(pixelIdx)));
            }
        }

        return;
    }
    catch (MathFit::CFitException& e)
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
        std::stringstream message;
        message << "Failed to evaluate" + m_name + ": a fit exception occurred.";
        if (nullptr != e.mMessage && strlen(e.mMessage) > 0)
        {
            message << " Message: '" + std::string{ e.mMessage } + "'";
        }

        throw DoasFitException(message.str());
    }
}

}