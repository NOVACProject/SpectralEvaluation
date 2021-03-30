#include <SpectralEvaluation/Math/FunctionFit.h>
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
ILF_RETURN_CODE FitFunction(MathFit::CVector& xData, MathFit::CVector& yData, T& functionToFit)
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
            return ILF_RETURN_CODE::FIT_FAILURE;
        }
        fit.FinishMinimize();

        return ILF_RETURN_CODE::SUCCESS;
    }
    catch (MathFit::CFitException& e)
    {
        std::cout << "Fit failed: " << e.mMessage << std::endl;

        return ILF_RETURN_CODE::SUCCESS;
    }
}

void CreateInitialEstimate(const std::vector<double>& x, const std::vector<double>& y, MathFit::CGaussFunction& result)
{
    // Find the amplitude of the Gaussian
    size_t idxOfMaximumAmplitude = 0;
    double maximumAmplitude = y[0];
    for (size_t ii = 1; ii < y.size(); ++ii)
    {
        if (y[ii] > maximumAmplitude)
        {
            maximumAmplitude = y[ii];
            idxOfMaximumAmplitude = ii;
        }
    }
    result.SetScale(maximumAmplitude);

    if (idxOfMaximumAmplitude <= 2 || idxOfMaximumAmplitude >= y.size() - 2)
    {
        return; // the algorithm below will not work for this case, skip
    }

    // Estimate the Full Width at Half Maximum (FWHM) by finding the points where the amplitude has dropped by half
    size_t halfAmplitudeLeft = idxOfMaximumAmplitude;
    for (size_t ii = idxOfMaximumAmplitude - 1; ii > 1; --ii)
    {
        if (y[ii] < maximumAmplitude * 0.5)
        {
            halfAmplitudeLeft = ii;
            break;
        }
    }
    size_t halfAmplitudeRight = idxOfMaximumAmplitude;
    for (size_t ii = idxOfMaximumAmplitude + 1; ii < y.size(); ++ii)
    {
        if (y[ii] < maximumAmplitude * 0.5)
        {
            halfAmplitudeRight = ii;
            break;
        }
    }
    const double fwhm = x[halfAmplitudeRight] - x[halfAmplitudeLeft];
    const double estimatedSigma = fwhm / 2.35482;

    result.SetSigma(estimatedSigma);

    // Estimate the center as the center-of-mass of the entire dataset
    double sumOfWeights = 0.0;
    double weightedSum = 0.0;
    for (size_t ii = 0; ii < y.size(); ++ii)
    {
        weightedSum += x[ii] * y[ii];
        sumOfWeights += y[ii];
    }
    const double centerOfMass = weightedSum / sumOfWeights;
    result.SetCenter(centerOfMass);
}

ILF_RETURN_CODE FitGaussian(std::vector<double>& x, std::vector<double>& y, MathFit::CGaussFunction& gaussian)
{
    const bool autoReleaseData = false;
    MathFit::CVector xData{ &x[0], static_cast<int>(x.size()), 1, autoReleaseData };
    MathFit::CVector yData{ &y[0], static_cast<int>(x.size()), 1, autoReleaseData };

    // First create an initial estimation of the location, width and amplitude of the Gaussian
    CreateInitialEstimate(x, y, gaussian);

    return FitFunction(xData, yData, gaussian);
}

ILF_RETURN_CODE FitGaussian(std::vector<double>& y, MathFit::CGaussFunction& gaussian)
{
    std::vector<double> x;
    x.resize(y.size());
    std::iota(begin(x), end(x), 0);

    return FitGaussian(x, y, gaussian);
}

}
