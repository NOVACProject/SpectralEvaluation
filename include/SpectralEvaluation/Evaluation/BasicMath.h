// BasicMath.h: interface for the CBasicMath class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BASICMATH_H__1DEB20E2_5D81_11D4_866C_00E098701FA6__INCLUDED_)
#define AFX_BASICMATH_H__1DEB20E2_5D81_11D4_866C_00E098701FA6__INCLUDED_

namespace MathFit
{
class CVector;
}

#include "../Fit/FitException.h"

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

using namespace MathFit;

class CBasicMath
{
public:
    void FFT(double* fData, double* fReal, double* fImaginary, int iLength);

    void InverseFFT(double* fReal, double* fImaginary, double* fData, int iLength);

    bool GaussFit(double* fXData, int iNXValues, double* fYData, double& fCenter, double& fSigma, double& fScale);

    bool PolynomialFit(double* fXData, int iNXValues, double* fYData, double* fCoeff, int iOrder);

    void CrossCorrelate(double* fFirst, int iLengthFirst, double* fSec, int iLengthSec);

    void FillGauss(double* fData, int iSize, double fA, double fSigma, double fScale = 0);

    void FillRandom(double* fData, int iLength, double fVariance);

    enum { NOWEIGHT = 0, SCANWEIGHT, TIMEWEIGHT };

    double FitRes2MicroGrammPerCubicMeter(double fFitResult, double fLightPathLength, double fMolecularWeight);

    double FitRes2ppb(double fFitResult, double fLightPathLength, double fTemperature, double fPreasure);

    bool ShiftAndSqueeze(CVector& vXData, CVector& vYData, double fOrigin, double fShift, double fSqueeze);

    void Div(double* fFirst, int iSize, double fConst);

    void Mul(double* fFirst, int iSize, double fConst);

    void Sub(double* fFirst, int iSize, double fConst);

    void Add(double* fFirst, int iSize, double fConst);

    void Div(double* fFirst, const double* fSec, int iSize, double fFactor = 0.0);

    void Mul(double* fFirst, const double* fSec, int iSize, double fFactor = 0.0);

    void Sub(double* fFirst, const double* fSec, int iSize, double fFactor = 0.0);

    void Add(double* fFirst, const double* fSec, int iSize, double fFactor = 0.0);

    void Reverse(double* fData, int iSize);

    void Convolute(double* fFirst, int iSize, double* fCore, int iCoreSize);

    void Invert(double* fData, int iSize);

    void Reciprocal(double* fData, int iSize);

    void NormalizeEnergy(double* fData, int iSize, double fEnergy);

    void NormalizeAmplitude(double* fData, int iSize, double fMaxAmplitude);

    void Zero(double* fData, int iSize, double fZeroLimit);

    void BiasAdjust(double* fData, int iSize);

    double CalcExposureTime(double fSaturation, double fEOffsetExpTime, double fCurrentExpTime, double fCurrentAvg, double fBackExpTime, double fBackAvg, double fIntTime);

    double* CalcMeasuredSpec(double* fRes, double* fMea, int iMeaScans, double fMeaExpTime, double* fLamp, int iLampScans, double fLampExpTime, double* fBack, int iBackScans, double fBackExpTime, double fOffset, double fOffsetExpTime, int iSize);

    double* Log(double* fData, int iSize);

    double* Delog(double* fData, int iSize);

    double* HighPassBinomial(double* fData, int iSize, int iNIterations);

    double* LowPassBinomial(double* fData, int iSize, int iNIterations);

    CBasicMath();

    virtual ~CBasicMath();

private:
    void DFourier(double data[], unsigned long nn, int isign);

    static bool mDoNotUseMathLimits;
};

#endif // !defined(AFX_BASICMATH_H__1DEB20E2_5D81_11D4_866C_00E098701FA6__INCLUDED_)

