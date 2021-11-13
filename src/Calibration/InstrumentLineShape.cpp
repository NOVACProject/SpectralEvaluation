#include <SpectralEvaluation/Calibration/InstrumentLineShape.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Fit/StandardFit.h>
#include <SpectralEvaluation/Fit/StandardMetricFunction.h>
#include <SpectralEvaluation/Fit/GaussFunction.h>
#include <SpectralEvaluation/FitExtensions/AsymmetricGaussFunction.h>
#include <SpectralEvaluation/FitExtensions/SuperGaussFunction.h>
#include <SpectralEvaluation/Fit/CubicSplineFunction.h>
#include <SpectralEvaluation/VectorUtils.h>
#include <numeric>

namespace novac
{

double SuperGaussianLineShape::Fwhm() const { return 2.0 * w * std::pow(0.69314718056, 1.0 / k); }

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
    else if (mercuryLine.m_wavelength.size() != static_cast<size_t>(mercuryLine.m_length))
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
    else if (mercuryLine.m_wavelength.size() != static_cast<size_t>(mercuryLine.m_length))
    {
        return FUNCTION_FIT_RETURN_CODE::MISSING_WAVELENGTH_CALIBRATION;
    }

    std::vector<double> localY{ mercuryLine.m_data, mercuryLine.m_data + mercuryLine.m_length };

    // First create an initial estimation of the location, width and amplitude of the Gaussian
    MathFit::CGaussFunction simpleGaussian;
    CreateInitialEstimate(mercuryLine.m_wavelength, localY, simpleGaussian);


    const bool autoReleaseData = false;
    std::vector<double> localX{ mercuryLine.m_wavelength };
    MathFit::CVector xData{ &localX[0], static_cast<int>(mercuryLine.m_length), 1, autoReleaseData };
    MathFit::CVector yData{ &localY[0], static_cast<int>(mercuryLine.m_length), 1, autoReleaseData };


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
    else if (mercuryLine.m_wavelength.size() != static_cast<size_t>(mercuryLine.m_length))
    {
        return FUNCTION_FIT_RETURN_CODE::MISSING_WAVELENGTH_CALIBRATION;
    }

    std::vector<double> localX{ mercuryLine.m_wavelength };
    std::vector<double> localY{ mercuryLine.m_data, mercuryLine.m_data + mercuryLine.m_length };

    const bool autoReleaseData = false;
    MathFit::CVector xData{ &localX[0], static_cast<int>(mercuryLine.m_length), 1, autoReleaseData };
    MathFit::CVector yData{ &localY[0], static_cast<int>(mercuryLine.m_length), 1, autoReleaseData };

    // First create an initial estimation of the location, width and amplitude of the Gaussian
    MathFit::CGaussFunction regularGaussian;
    CreateInitialEstimate(localX, localY, regularGaussian);

    // Copy the relevant parameters from the regular Gaussian and proceed with fitting the super-gaussian
    MathFit::CSuperGaussFunction superGaussian;
    superGaussian.SetCenter(regularGaussian.GetCenter());
    superGaussian.SetW(regularGaussian.GetSigma());

    FUNCTION_FIT_RETURN_CODE ret = FitFunction(xData, yData, superGaussian);
    if (ret != FUNCTION_FIT_RETURN_CODE::SUCCESS)
    {
        return ret;
    }
    else
    {
        result.w = superGaussian.GetW();
        result.k = superGaussian.GetK();
        result.center = superGaussian.GetCenter();

        return FUNCTION_FIT_RETURN_CODE::SUCCESS;
    }
}

// TODO: Merge with overload above!
FUNCTION_FIT_RETURN_CODE FitInstrumentLineShape(const CCrossSectionData& mercuryLine, SuperGaussianLineShape& result)
{
    if (mercuryLine.GetSize() == 0)
    {
        return FUNCTION_FIT_RETURN_CODE::EMPTY_INPUT;
    }
    else if (mercuryLine.m_waveLength.size() != mercuryLine.GetSize())
    {
        return FUNCTION_FIT_RETURN_CODE::MISSING_WAVELENGTH_CALIBRATION;
    }

    std::vector<double> localX{ mercuryLine.m_waveLength };
    std::vector<double> localY{ mercuryLine.m_crossSection.data(), mercuryLine.m_crossSection.data() + mercuryLine.GetSize() };

    const bool autoReleaseData = false;
    MathFit::CVector xData{ &localX[0], static_cast<int>(mercuryLine.GetSize()), 1, autoReleaseData };
    MathFit::CVector yData{ &localY[0], static_cast<int>(mercuryLine.GetSize()), 1, autoReleaseData };

    // First create an initial estimation of the location, width and amplitude of the Gaussian
    MathFit::CGaussFunction regularGaussian;
    CreateInitialEstimate(localX, localY, regularGaussian);

    // Copy the relevant parameters from the regular Gaussian and proceed with fitting the super-gaussian
    MathFit::CSuperGaussFunction superGaussian;
    superGaussian.SetCenter(regularGaussian.GetCenter());
    superGaussian.SetW(regularGaussian.GetSigma());

    FUNCTION_FIT_RETURN_CODE ret = FitFunction(xData, yData, superGaussian);
    if (ret != FUNCTION_FIT_RETURN_CODE::SUCCESS)
    {
        return ret;
    }
    else
    {
        result.w = superGaussian.GetW();
        result.k = superGaussian.GetK();
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
    superGauss.SetW(lineShape.w);
    superGauss.SetK(lineShape.k);
    superGauss.SetScale(amplitude);

    return GetFunctionValues(superGauss, x, baseline);
}

CCrossSectionData SampleInstrumentLineShape(const GaussianLineShape& lineShape)
{
    // Use the fact that a Gaussian line shape is just a special case of a super-gaussian.
    SuperGaussianLineShape superGaussianLineShape;
    superGaussianLineShape.k = 2.0;
    superGaussianLineShape.w = lineShape.sigma * std::sqrt(2.0);

    return SampleInstrumentLineShape(superGaussianLineShape);
}

CCrossSectionData SampleInstrumentLineShape(const SuperGaussianLineShape& lineShape)
{
    const double amplitude = lineShape.k / (2.0 * lineShape.w * std::lgamma(1. / lineShape.k));

    MathFit::CSuperGaussFunction superGauss;
    superGauss.SetCenter(0.0);
    superGauss.SetW(lineShape.w);
    superGauss.SetK(lineShape.k);
    superGauss.SetScale(amplitude);

    CCrossSectionData result;

    const size_t length = 101;
    const double range = 3 * lineShape.Fwhm();
    const double xMin = -range * 0.5;

    result.m_waveLength.resize(length);
    for (size_t ii = 0; ii < length; ++ii)
    {
        result.m_waveLength[ii] = xMin + range * ii / (double)(length - 1);
    }

    result.m_crossSection = GetFunctionValues(superGauss, result.m_waveLength, 0.0);

    // Make sure that the area below the curve equals one.
    NormalizeArea(result.m_crossSection, result.m_crossSection);

    return result;
}

std::vector<double> PartialDerivative(const GaussianLineShape& lineShape, const std::vector<double>& x)
{
    const double center = 0.0;
    const double amplitude = 1.0;
    const double baseline = 0.0;
    const double delta = 0.01;

    MathFit::CGaussFunction originalLineShape;
    originalLineShape.SetCenter(center);
    originalLineShape.SetSigma(lineShape.sigma + delta);
    originalLineShape.SetScale(amplitude);
    const auto x0 = GetFunctionValues(originalLineShape, x, baseline);

    MathFit::CGaussFunction modifiedLineShape;
    modifiedLineShape.SetCenter(center);
    modifiedLineShape.SetSigma(lineShape.sigma + delta);
    modifiedLineShape.SetScale(amplitude);
    const auto x1 = GetFunctionValues(modifiedLineShape, x, baseline);

    std::vector<double> result;
    result.resize(x0.size());

    for (size_t ii = 0; ii < result.size(); ++ii)
    {
        result[ii] = (x1[ii] - x0[ii]) / delta;
    }

    return result;
}

void Setup(const SuperGaussianLineShape& src, MathFit::CSuperGaussFunction& dst)
{
    const double center = 0.0;
    const double amplitude = 1.0;

    dst.SetCenter(center);
    dst.SetW(src.w);
    dst.SetK(src.k);
    dst.SetScale(amplitude);
}

double GetParameterValue(const SuperGaussianLineShape& lineShape, int parameterIdx)
{
    if (parameterIdx == 0)
    {
        return lineShape.w;
    }
    else
    {
        return lineShape.k;
    }
}

void SetParameterValue(SuperGaussianLineShape& lineShape, int parameterIdx, double value)
{
    if (parameterIdx == 0)
    {
        lineShape.w = value;
    }
    else
    {
        lineShape.k = value;
    }
}

/// <summary>
/// Normalizes the curve such that the area under abs(data) is one,
///     assuming that there is no baseline offset in the data.
/// </summary>
void NormalizeAbsArea(std::vector<double>& data)
{
    if (data.size() == 0)
    {
        return;
    }

    const double sumOfValues = SumAbs(data);

    for (size_t ii = 0; ii < data.size(); ++ii)
    {
        data[ii] = data[ii] / sumOfValues;
    }

    assert(std::abs(SumAbs(data) - 1.0) < 0.1);
}

std::vector<double> PartialDerivative(const SuperGaussianLineShape& lineShape, const std::vector<double>& x, int parameter)
{
    const double baseline = 0.0;
    const double delta = 0.01;

    MathFit::CSuperGaussFunction originalLineShape;
    Setup(lineShape, originalLineShape);
    auto x0 = GetFunctionValues(originalLineShape, x, baseline);
    NormalizeArea(x0, x0);

    MathFit::CSuperGaussFunction forwardLineShape;
    Setup(lineShape, forwardLineShape);
    if (parameter == 0)
    {
        forwardLineShape.SetW(lineShape.w + delta);
    }
    else
    {
        forwardLineShape.SetK(lineShape.k + delta);
    }
    auto xF = GetFunctionValues(forwardLineShape, x, baseline);
    NormalizeArea(xF, xF);

    std::vector<double> result;
    result.resize(x0.size());

    for (size_t ii = 0; ii < result.size(); ++ii)
    {
        result[ii] = (xF[ii] - x0[ii]) / (delta);
    }

    return result;
}

}










