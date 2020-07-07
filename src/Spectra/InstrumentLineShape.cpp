#include <SpectralEvaluation/Spectra/InstrumentLineShape.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Fit/StandardFit.h>
#include <SpectralEvaluation/Fit/StandardMetricFunction.h>
#include <SpectralEvaluation/Fit/GaussFunction.h>
#include <SpectralEvaluation/FitExtensions/AsymmetricGaussFunction.h>
#include <SpectralEvaluation/FitExtensions/SuperGaussFunction.h>
#include <SpectralEvaluation/Fit/CubicSplineFunction.h>

namespace Calibration
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


ILF_RETURN_CODE FitInstrumentLineShape(const CSpectrum& mercuryLine, GaussianLineShape& result)
{
    if (mercuryLine.m_length == 0)
    {
        return ILF_RETURN_CODE::EMPTY_INPUT;
    }
    else if (mercuryLine.m_wavelength.size() != mercuryLine.m_length)
    {
        return ILF_RETURN_CODE::MISSING_WAVELENGTH_CALIBRATION;
    }

    std::vector<double> localX{ mercuryLine.m_wavelength };
    std::vector<double> localY{ mercuryLine.m_data, mercuryLine.m_data + mercuryLine.m_length };

    const bool autoReleaseData = false;
    MathFit::CVector xData{ &localX[0], mercuryLine.m_length, 1, autoReleaseData };
    MathFit::CVector yData{ &localY[0], mercuryLine.m_length, 1, autoReleaseData };

    MathFit::CGaussFunction gaussianToFit;

    ILF_RETURN_CODE ret = FitFunction(xData, yData, gaussianToFit);
    if (ret != ILF_RETURN_CODE::SUCCESS)
    {
        return ret;
    }
    else
    {
        result.sigma = gaussianToFit.GetSigma();

        return ILF_RETURN_CODE::SUCCESS;
    }
}

ILF_RETURN_CODE FitInstrumentLineShape(const CSpectrum& mercuryLine, AsymmetricGaussianLineShape& result)
{
    if (mercuryLine.m_length == 0)
    {
        return ILF_RETURN_CODE::EMPTY_INPUT;
    }
    else if (mercuryLine.m_wavelength.size() != mercuryLine.m_length)
    {
        return ILF_RETURN_CODE::MISSING_WAVELENGTH_CALIBRATION;
    }

    std::vector<double> localX{ mercuryLine.m_wavelength };
    std::vector<double> localY{ mercuryLine.m_data, mercuryLine.m_data + mercuryLine.m_length };

    const bool autoReleaseData = false;
    MathFit::CVector xData{ &localX[0], mercuryLine.m_length, 1, autoReleaseData };
    MathFit::CVector yData{ &localY[0], mercuryLine.m_length, 1, autoReleaseData };

    MathFit::CAsymmetricGaussFunction gaussianToFit;

    ILF_RETURN_CODE ret = FitFunction(xData, yData, gaussianToFit);
    if (ret != ILF_RETURN_CODE::SUCCESS)
    {
        return ret;
    }
    else
    {
        result.sigmaLeft = gaussianToFit.GetSigmaLeft();
        result.sigmaRight = gaussianToFit.GetSigmaRight();

        return ILF_RETURN_CODE::SUCCESS;
    }
}


ILF_RETURN_CODE FitInstrumentLineShape(const CSpectrum& mercuryLine, SuperGaussianLineShape& result)
{
    if (mercuryLine.m_length == 0)
    {
        return ILF_RETURN_CODE::EMPTY_INPUT;
    }
    else if (mercuryLine.m_wavelength.size() != mercuryLine.m_length)
    {
        return ILF_RETURN_CODE::MISSING_WAVELENGTH_CALIBRATION;
    }

    std::vector<double> localX{ mercuryLine.m_wavelength };
    std::vector<double> localY{ mercuryLine.m_data, mercuryLine.m_data + mercuryLine.m_length };

    const bool autoReleaseData = false;
    MathFit::CVector xData{ &localX[0], mercuryLine.m_length, 1, autoReleaseData };
    MathFit::CVector yData{ &localY[0], mercuryLine.m_length, 1, autoReleaseData };

    // Start by fitting a "regular" Gaussian function
    MathFit::CGaussFunction regularGaussian;
    ILF_RETURN_CODE ret = FitFunction(xData, yData, regularGaussian);
    if (ret != ILF_RETURN_CODE::SUCCESS)
    {
        return ret;
    }

    // Copy the relevant parameters from the regular Gaussian and proceed with fitting the super-gaussian
    MathFit::CSuperGaussFunction superGaussian;
    superGaussian.SetCenter(regularGaussian.GetCenter());
    superGaussian.SetSigma(regularGaussian.GetSigma());

    ret = FitFunction(xData, yData, superGaussian);
    if (ret != ILF_RETURN_CODE::SUCCESS)
    {
        return ret;
    }
    else
    {
        result.sigma = superGaussian.GetSigma();
        result.P = superGaussian.GetPower();

        return ILF_RETURN_CODE::SUCCESS;
    }
}

}










