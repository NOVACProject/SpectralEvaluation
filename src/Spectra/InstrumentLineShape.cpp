#include <SpectralEvaluation/Spectra/InstrumentLineShape.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Fit/StandardFit.h>
#include <SpectralEvaluation/Fit/StandardMetricFunction.h>
#include <SpectralEvaluation/Fit/GaussFunction.h>
#include <SpectralEvaluation/Fit/CubicSplineFunction.h>

namespace Evaluation
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


}
