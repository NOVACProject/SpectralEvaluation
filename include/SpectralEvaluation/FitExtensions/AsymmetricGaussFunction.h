/**
* Contains a defintion of an asymmetrical Gauss function object.
*
* @author		Mattias Johansson
* @version		2020/06/29
*/

#pragma once

#include "../Fit/ParamFunction.h"

namespace MathFit
{
/**
* This object represents an asymmetrical Gauss function.
  This function is composed of two Gaussian functions with one width (sigma) for x < center and one width for x >= center
*/
class CAsymmetricGaussFunction : public IParamFunction
{
private:
    const int LinearParamIdx_Scale = 0;

    const int NonLinearParamIdx_Center = 0;
    const int NonLinearParamIdx_SigmaNegative = 1;
    const int NonLinearParamIdx_SigmaPositive = 2;

public:
    /**
    * Creates the object and sets the default parameter values.
    *
    * @param bNormAmp	If TRUE the area under the function will always be one.
    */
    CAsymmetricGaussFunction(bool bNormAmp = false) : fSqrPI((TFitData)sqrt(2 * MATHFIT_PI))
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
        TFitData fVal = fXValue - GetCenter();
        if (fVal < 0.0)
        {
            fVal *= fVal;
            fVal /= (GetSigmaLeft() * GetSigmaLeft());
        }
        else
        {
            fVal *= fVal;
            fVal /= (GetSigmaRight() * GetSigmaRight());
        }
        fVal *= -0.5;

        return GetScale() * (TFitData)exp(fVal);
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
        const TFitData sigma2Left = (TFitData)-0.5 / (GetSigmaLeft() * GetSigmaLeft());
        const TFitData sigma2Right = (TFitData)-0.5 / (GetSigmaRight() * GetSigmaRight());
        const int iXSize = vXValues.GetSize();
        const double center = GetCenter();
        const double scale = GetScale();

        for (int i = 0; i < iXSize; i++)
        {
            TFitData delta = vXValues.GetAt(i) - center;
            if (delta < 0.0)
            {
                delta *= delta;
                vYTargetVector.SetAt(i, scale * (TFitData)exp(delta * sigma2Left));
            }
            else
            {
                delta *= delta;
                vYTargetVector.SetAt(i, scale * (TFitData)exp(delta * sigma2Right));
            }
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
    * @param fXValue	The X data point at which the first derivation is needed.
    *
    * @return	The slope of the reference at the given data point.
    */
    virtual TFitData GetSlope(TFitData fXValue)
    {
        TFitData fVal = GetValue(fXValue);

        fVal *= -(fXValue - GetCenter());
        fVal /= GetSigma(fXValue) * GetSigma(fXValue);

        return fVal;
    }

    /**
    * Calculates the first derivative of the function at a set of given data points.
    *
    * @param vXValues		A vector object containing the X values at which the function has to be evaluated.
    * @param vSlopeVector	A vector object which receives the resulting function values.
    *
    * @return	A reference to the slope vector object.
    *
    * @see GetSlope
    */
    virtual CVector& GetSlopes(CVector& vXValues, CVector& vSlopeVector)
    {
        GetValues(vXValues, vSlopeVector);

        const TFitData sigma2Left = GetSigmaLeft() * GetSigmaLeft();
        const TFitData sigma2Right = GetSigmaRight() * GetSigmaRight();

        const int iXSize = vXValues.GetSize();
        const double center = GetCenter();

        for (int i = 0; i < iXSize; i++)
        {
            const double diff = vXValues.GetAt(i) - center;
            if (diff < 0)
            {
                vSlopeVector.SetAt(i, -vSlopeVector.GetAt(i) * (vXValues.GetAt(i) - center) / sigma2Left);
            }
            else
            {
                vSlopeVector.SetAt(i, -vSlopeVector.GetAt(i) * (vXValues.GetAt(i) - center) / sigma2Right);
            }
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
        TFitData fVal = fXValue - GetCenter();
        if (fVal < 0.0)
        {
            fVal *= fVal;
            fVal /= (GetSigmaLeft() * GetSigmaLeft());
        }
        else
        {
            fVal *= fVal;
            fVal /= (GetSigmaRight() * GetSigmaRight());
        }
        fVal *= -0.5;

        return (TFitData)exp(fVal);
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
    * Returns the sigma parameter on the left side of the center
    */
    TFitData GetSigmaLeft()
    {
        return mNonlinearParams.GetAllParameter().GetAt(NonLinearParamIdx_SigmaNegative);
    }

    /**
    * Returns the sigma parameter on the right side of the center
    */
    TFitData GetSigmaRight()
    {
        return mNonlinearParams.GetAllParameter().GetAt(NonLinearParamIdx_SigmaPositive);
    }

    /**
    * Returns the sigma parameter at the requested x-axis point
    */
    TFitData GetSigma(double x)
    {
        if (x < GetCenter())
        {
            return mNonlinearParams.GetAllParameter().GetAt(NonLinearParamIdx_SigmaNegative);
        }
        else
        {
            return mNonlinearParams.GetAllParameter().GetAt(NonLinearParamIdx_SigmaPositive);
        }
    }

    /**
    * Returns the amplitude scal factor.
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
            TFitData fVal = fSqrPI * 0.5 * (GetSigmaLeft() + GetSigmaRight());
            return 1 / fVal;
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
    * Sets the sigma value for x < center.
    */
    void SetSigmaLeft(TFitData newValue)
    {
        mNonlinearParams.SetParameter(NonLinearParamIdx_SigmaNegative, newValue);
    }

    /**
    * Sets the sigma value for x >= center
    */
    void SetSigmaRight(TFitData newValue)
    {
        mNonlinearParams.SetParameter(NonLinearParamIdx_SigmaPositive, newValue);
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
        mNonlinearParams.SetParameter(NonLinearParamIdx_Center, 0);
        mNonlinearParams.SetParameter(NonLinearParamIdx_SigmaNegative, 1);
        mNonlinearParams.SetParameter(NonLinearParamIdx_SigmaPositive, 1);
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

