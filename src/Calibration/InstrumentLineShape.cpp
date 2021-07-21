#include <SpectralEvaluation/Calibration/InstrumentLineShape.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Fit/StandardFit.h>
#include <SpectralEvaluation/Fit/StandardMetricFunction.h>
#include <SpectralEvaluation/Fit/GaussFunction.h>
#include <SpectralEvaluation/FitExtensions/AsymmetricGaussFunction.h>
#include <SpectralEvaluation/FitExtensions/SuperGaussFunction.h>
#include <SpectralEvaluation/Fit/CubicSplineFunction.h>
#include <numeric>

namespace novac
{

template<class T>
FUNCTION_FIT_RETURN_CODE FitFunction(MathFit::CVector& xData, MathFit::CVector& yData, T& functionToFit)
{
    try
    {
        MathFit::CCubicSplineFunction cubicSplineRepresentation{ xData, yData };

        MathFit::CStandardMetricFunction diff(cubicSplineRepresentation, functionToFit);

        MathFit::CStandardFit fit(diff);
        fit.SetFitRange(xData);
        fit.SetMaxFitSteps(500);
        fit.PrepareMinimize();

        if (!fit.Minimize())
        {
            return FUNCTION_FIT_RETURN_CODE::FIT_FAILURE;
        }
        fit.FinishMinimize();

        return FUNCTION_FIT_RETURN_CODE::SUCCESS;
    }
    catch (MathFit::CFitException& e)
    {
        std::cout << "Fit failed: " << e.mMessage << std::endl;

        return FUNCTION_FIT_RETURN_CODE::SUCCESS;
    }
}

FUNCTION_FIT_RETURN_CODE FitInstrumentLineShape(const CSpectrum& mercuryLine, GaussianLineShape& result)
{
    if (mercuryLine.m_length == 0)
    {
        return FUNCTION_FIT_RETURN_CODE::EMPTY_INPUT;
    }
    else if (mercuryLine.m_wavelength.size() != mercuryLine.m_length)
    {
        return FUNCTION_FIT_RETURN_CODE::MISSING_WAVELENGTH_CALIBRATION;
    }

    std::vector<double> localX{ mercuryLine.m_wavelength };
    std::vector<double> localY{ mercuryLine.m_data, mercuryLine.m_data + mercuryLine.m_length };

    MathFit::CGaussFunction gaussianToFit;
    FUNCTION_FIT_RETURN_CODE ret = FitGaussian(localX, localY, gaussianToFit);

    if (ret != FUNCTION_FIT_RETURN_CODE::SUCCESS)
    {
        return ret;
    }
    else
    {
        result.sigma = gaussianToFit.GetSigma();
        result.center = gaussianToFit.GetCenter();

        return FUNCTION_FIT_RETURN_CODE::SUCCESS;
    }
}

FUNCTION_FIT_RETURN_CODE FitInstrumentLineShape(const CSpectrum& mercuryLine, AsymmetricGaussianLineShape& result)
{
    if (mercuryLine.m_length == 0)
    {
        return FUNCTION_FIT_RETURN_CODE::EMPTY_INPUT;
    }
    else if (mercuryLine.m_wavelength.size() != mercuryLine.m_length)
    {
        return FUNCTION_FIT_RETURN_CODE::MISSING_WAVELENGTH_CALIBRATION;
    }

    std::vector<double> localY{ mercuryLine.m_data, mercuryLine.m_data + mercuryLine.m_length };

    // First create an initial estimation of the location, width and amplitude of the Gaussian
    MathFit::CGaussFunction simpleGaussian;
    CreateInitialEstimate(mercuryLine.m_wavelength, localY, simpleGaussian);


    const bool autoReleaseData = false;
    std::vector<double> localX{ mercuryLine.m_wavelength };
    MathFit::CVector xData{ &localX[0], mercuryLine.m_length, 1, autoReleaseData };
    MathFit::CVector yData{ &localY[0], mercuryLine.m_length, 1, autoReleaseData };


    MathFit::CAsymmetricGaussFunction gaussianToFit;
    gaussianToFit.SetCenter(simpleGaussian.GetCenter());
    gaussianToFit.SetScale(simpleGaussian.GetScale());
    gaussianToFit.SetSigmaLeft(simpleGaussian.GetSigma());
    gaussianToFit.SetSigmaRight(simpleGaussian.GetSigma());

    FUNCTION_FIT_RETURN_CODE ret = FitFunction(xData, yData, gaussianToFit);
    if (ret != FUNCTION_FIT_RETURN_CODE::SUCCESS)
    {
        return ret;
    }
    else
    {
        result.sigmaLeft = gaussianToFit.GetSigmaLeft();
        result.sigmaRight = gaussianToFit.GetSigmaRight();

        return FUNCTION_FIT_RETURN_CODE::SUCCESS;
    }
}

FUNCTION_FIT_RETURN_CODE FitInstrumentLineShape(const CSpectrum& mercuryLine, SuperGaussianLineShape& result)
{
    if (mercuryLine.m_length == 0)
    {
        return FUNCTION_FIT_RETURN_CODE::EMPTY_INPUT;
    }
    else if (mercuryLine.m_wavelength.size() != mercuryLine.m_length)
    {
        return FUNCTION_FIT_RETURN_CODE::MISSING_WAVELENGTH_CALIBRATION;
    }

    std::vector<double> localX{ mercuryLine.m_wavelength };
    std::vector<double> localY{ mercuryLine.m_data, mercuryLine.m_data + mercuryLine.m_length };

    const bool autoReleaseData = false;
    MathFit::CVector xData{ &localX[0], mercuryLine.m_length, 1, autoReleaseData };
    MathFit::CVector yData{ &localY[0], mercuryLine.m_length, 1, autoReleaseData };

    // First create an initial estimation of the location, width and amplitude of the Gaussian
    MathFit::CGaussFunction regularGaussian;
    CreateInitialEstimate(localX, localY, regularGaussian);

    // Copy the relevant parameters from the regular Gaussian and proceed with fitting the super-gaussian
    MathFit::CSuperGaussFunction superGaussian;
    superGaussian.SetCenter(regularGaussian.GetCenter());
    superGaussian.SetSigma(regularGaussian.GetSigma());

    FUNCTION_FIT_RETURN_CODE ret = FitFunction(xData, yData, superGaussian);
    if (ret != FUNCTION_FIT_RETURN_CODE::SUCCESS)
    {
        return ret;
    }
    else
    {
        result.sigma = superGaussian.GetSigma();
        result.P = superGaussian.GetPower();
        result.center = superGaussian.GetCenter();

        return FUNCTION_FIT_RETURN_CODE::SUCCESS;
    }
}

template<class T>
std::vector<double> GetFunctionValues(T& function, const std::vector<double>& x, double baseline = 0.0)
{
    std::vector<double> y(x.size());

    for (size_t ii = 0; ii < x.size(); ++ii)
    {
        y[ii] = baseline + function.GetValue(x[ii]);
    }

    return y;
}

std::vector<double> SampleInstrumentLineShape(const GaussianLineShape& lineShape, const std::vector<double>& x, double center, double amplitude, double baseline)
{
    MathFit::CGaussFunction gaussian;
    gaussian.SetCenter(center);
    gaussian.SetSigma(lineShape.sigma);
    gaussian.SetScale(amplitude);

    return GetFunctionValues(gaussian, x, baseline);
}

std::vector<double> SampleInstrumentLineShape(const AsymmetricGaussianLineShape& lineShape, const std::vector<double>& x, double center, double amplitude, double baseline)
{
    MathFit::CAsymmetricGaussFunction gaussian;
    gaussian.SetCenter(center);
    gaussian.SetSigmaLeft(lineShape.sigmaLeft);
    gaussian.SetSigmaRight(lineShape.sigmaRight);
    gaussian.SetScale(amplitude);

    return GetFunctionValues(gaussian, x, baseline);
}

std::vector<double> SampleInstrumentLineShape(const SuperGaussianLineShape& lineShape, const std::vector<double>& x, double center, double amplitude, double baseline)
{
    MathFit::CSuperGaussFunction superGauss;
    superGauss.SetCenter(center);
    superGauss.SetSigma(lineShape.sigma);
    superGauss.SetPower(lineShape.P);
    superGauss.SetScale(amplitude);

    return GetFunctionValues(superGauss, x, baseline);
}


}










