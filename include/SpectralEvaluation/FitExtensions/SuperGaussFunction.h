/**
* Contains a defintion of an super-Gaussian function object.
*
* @author		Mattias Johansson
* @version		2021/09/04
*/

#pragma once

#include <SpectralEvaluation/Fit/ParamFunction.h>

namespace MathFit
{
/**
* This object represents a super-Gaussian function.
* The super-Gaussian function is a generalization of the Gaussian function but allows for a higher order power:
*   \begin{verbatim}f(x)=s*exp(-[(x-a)/w]^k)\end{verbatim}
* Where P=2 yields a regular Gaussian function.
*/
class CSuperGaussFunction : public IParamFunction
{
private:
    const int LinearParamIdx_Scale = 0;

    const int NonLinearParamIdx_Center = 0;
    const int NonLinearParamIdx_w = 1;
    const int NonLinearParamIdx_k = 2;

public:
    /**
    * Creates the object and sets the default parameter values.
    *
    * @param bNormAmp If TRUE the area under the function will always be one.
    */
    CSuperGaussFunction(bool bNormAmp = false)
        : fSqrPI((TFitData)sqrt(2 * MATHFIT_PI))
    {
        // set default parameters
        ResetNonlinearParameter();

        // we can only scale the amplitude, if we should not keep the area size fixed
        mNormAmp = bNormAmp;
        ResetLinearParameter();
    }

    /**
    * Returns the value of the reference spectrum at the given data point.
    *
    * @param fXValue    The X data point.
    *
    * @return   The value of the reference at the given data point.
    */
    virtual TFitData GetValue(TFitData fXValue)
    {
        const TFitData q = std::abs((fXValue - GetCenter()) / GetW());
        return GetScale() * (TFitData)exp(-std::pow(q, GetK()));
    }

    /**
    * Calculates the function values at a set of given data points.
    *
    * @param vXValues       A vector object containing the X values at which the function has to be evaluated.
    * @param vYTargetVector A vector object which receives the resulting function values.
    *
    * @return   A reference to the Y vector object.
    *
    * @see  GetValue
    */
    virtual CVector& GetValues(CVector& vXValues, CVector& vYTargetVector)
    {
        const int iXSize = vXValues.GetSize();
        const double center = GetCenter();
        const double scale = GetScale();
        const double w = GetW();
        const double k = GetK();

        for (int i = 0; i < iXSize; i++)
        {
            const TFitData q = std::abs((vXValues.GetAt(i) - center) / w);
            vYTargetVector.SetAt(i, scale * (TFitData)exp(-std::pow(q, k)));
        }

        return vYTargetVector;
    }

    /**
    * Returns the first derivative of the reference spectrum at the given data point.
    * The first derivative is given by
    *
    * \begin{verbatim}f'(x)=-s*(x-a)/o^2*exp(-1/2*(x-a)^2/o^2)\end{verbatim}
    *
    * where
    \begin{verbatim}
    c := concentration
    w := squeeze
    v := shift
    \end{verbatim}
    *
    * @param fXValue The X data point at which the first derivation is needed.
    *
    * @return The slope of the reference at the given data point.
    */
    virtual TFitData GetSlope(TFitData fXValue)
    {
        TFitData fVal = GetValue(fXValue);
        fVal *= -(fXValue - GetCenter());
        fVal *= GetK();
        fVal /= std::pow(GetW(), GetK());

        return fVal;
    }

    /**
    * Calculates the first derivative of the function at a set of given data points.
    *
    * @param vXValues       A vector object containing the X values at which the function has to be evaluated.
    * @param vSlopeVector   A vector object which receives the resulting function values.
    *
    * @return A reference to the slope vector object.
    *
    * @see GetSlope
    */
    virtual CVector& GetSlopes(CVector& vXValues, CVector& vSlopeVector)
    {
        const int iXSize = vXValues.GetSize();
        const double factor = GetK() / (std::pow(GetW(), GetK()));

        for (int i = 0; i < iXSize; i++)
        {
            TFitData fVal = GetValue(vXValues.GetAt(i));
            fVal *= -(vXValues.GetAt(i) - GetCenter());
            vSlopeVector.SetAt(i, fVal * factor);
        }

        return vSlopeVector;
    }

    /**
    * Returns the basis function of the specified linear parameter.
    * A basis function is defined as the term by which the linear parameter is multiplied.
    *
    * @param fXValue	The data point at which the basis function should be determined.
    * @param iParamID	The index within the linear parameter vector of the linear parameter.
    * @param bFixedID	If TRUE the given parameter ID is the parameter ID without all fixed parameter.
    *
    * @return	The basis function in regard to the given linear parameter.
    */
    virtual TFitData GetLinearBasisFunction(TFitData fXValue, int /*iParamID*/, bool /*bFixedID = true*/)
    {
        // get the function value without the scale factor
        const TFitData val = std::pow(std::abs((fXValue - GetCenter()) / GetW()), GetK());

        return (TFitData)exp(-val);
    }

    /**
    * Returns the center parameter of the Gauss Bell function.
    *
    * @return The center parameter of the Gauss Bell function.
    */
    TFitData GetCenter()
    {
        return mNonlinearParams.GetAllParameter().GetAt(NonLinearParamIdx_Center);
    }

    /**
    * Returns the width parameter (w).
    */
    TFitData GetW()
    {
        return mNonlinearParams.GetAllParameter().GetAt(NonLinearParamIdx_w);
    }

    /**
    * Returns the power parameter (k)
    */
    TFitData GetK()
    {
        return mNonlinearParams.GetAllParameter().GetAt(NonLinearParamIdx_k);
    }

    /**
    * Returns the amplitude scale factor.
    * If the \Ref{mNormAmp} flag is set, the amplitude is scaled so that
    * the area size of the function will always be one. Otherwise an
    * appropriate amplitude scale factor is used.
    *
    * @return	The amplitude scale factor.
    */
    TFitData GetScale()
    {
        if (!mNormAmp)
        {
            return mLinearParams.GetAllParameter().GetAt(LinearParamIdx_Scale);
        }
        else
        {
            return GetK() / (2.0 * GetW() * std::lgamma(1.0 / GetK()));
        }
    }

    /**
    * Sets the scale value of the function.
    *
    * @param fScale	The new scale value.
    */
    void SetScale(TFitData fScale)
    {
        mLinearParams.SetParameter(LinearParamIdx_Scale, fScale);
    }

    /**
    * Sets the center value of the function.
    *
    * @param fCenter The new center value.
    */
    void SetCenter(TFitData fCenter)
    {
        mNonlinearParams.SetParameter(NonLinearParamIdx_Center, fCenter);
    }

    /**
    * Sets the w value.
    */
    void SetW(TFitData newValue)
    {
        mNonlinearParams.SetParameter(NonLinearParamIdx_w, newValue);
    }

    /**
    * Sets the power value.
    */
    void SetK(TFitData newValue)
    {
        mNonlinearParams.SetParameter(NonLinearParamIdx_k, newValue);
    }

    /**
    * Resets the linear parameters to default values.
    * The default implementation sets every parameter to zero.
    */
    virtual void ResetLinearParameter()
    {
        mLinearParams.SetSize(1);

        if (mNormAmp)
        {
            mLinearParams.FixParameter(LinearParamIdx_Scale, 1);
        }
        else
        {
            mLinearParams.ReleaseParameter(LinearParamIdx_Scale);
        }
        mLinearParams.SetParameter(LinearParamIdx_Scale, 1);
    }

    /**
    * Resets the nonlinear parameters to default values.
    * The default implementation sets every parameter to zero.
    */
    virtual void ResetNonlinearParameter()
    {
        mNonlinearParams.SetSize(3);
        mNonlinearParams.SetParameter(NonLinearParamIdx_Center, 0.0); // default: center = 0.0
        mNonlinearParams.SetParameter(NonLinearParamIdx_w, 1.0);  // default: w = 1.0
        mNonlinearParams.SetParameter(NonLinearParamIdx_k, 2.0);  // default: k = 2.0 (standard Gaussian)
    }

private:
    /**
    * Indicates wheter the area size normalization is active or not.
    */
    bool mNormAmp;
    /**
    * Holds a needed scale factor.
    */
    const TFitData fSqrPI;
};
}

