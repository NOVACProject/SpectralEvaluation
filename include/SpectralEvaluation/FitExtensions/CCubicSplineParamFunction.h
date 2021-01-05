/**
* Contains a defintion of an cubic spline which also implements IParamFunction.
* This is a later addition to DOASIS after the original import.
* In current version of DOASIS the CCubicSplineFunction also implements IParamFunction
*
* @author		Mattias Johansson
* @version		2020/09/06
*/

#pragma once

#include "../Fit/ParamFunction.h"

namespace MathFit
{
/**
  This class implements a function that uses cubic spline interpolation
  <seealso cref="!:Press et al., Numerical Recipes, p. 113 ff." />
*/
class CCubicSplineParamFunction : public IParamFunction
{
private:
    CVector mY2ndDerivates;
    CVector mH;
    CVector mSlopeInvariant;
    CVector mDeltaHSquareLow;
    CVector mDeltaHSquareHigh;
    CVector mSlopeDeltaHSquareLow;
    CVector mSlopeDeltaHSquareHigh;

public:

    CCubicSplineParamFunction()
    {
    }

    /// Create a cubic spline object using the given data.
    ///             
    ///              @param vXValues			The vector containing the X values.
    ///              @param vYValues			The vector containing the Y values in regard to the X values.
    CCubicSplineParamFunction(CVector& vXValues, CVector& vYValues)
    {
        this->SetData(vXValues, vYValues);
    }

    /// Create a cubic spline object using the given data.
    ///             
    ///              @param vXValues			The vector containing the X values.
    ///              @param vYValues			The vector containing the Y values in regard to the X values.
    ///              @param vError				The vector containing the errors of the Y values. This vector will not be interpolated!
    CCubicSplineParamFunction(CVector& vXValues, CVector& vYValues, CVector& vError)
    {
        this->SetData(vXValues, vYValues, vError);
    }

    /// @param vXValues			The vector containing the X values.
    ///              @param vYValues			The vector containing the Y values in regard to the X values.
    ///             
    ///              @return	TRUE if successful, FALSE otherwise
    bool SetData(CVector& vXValues, CVector& vYValues, CVector& vError) override
    {
        return IParamFunction::SetData(vXValues, vYValues, vError) && this->InitializeSpline();
    }

    /// @param vXValues			The vector containing the X values.
    ///              @param vYValues			The vector containing the Y values in regard to the X values.
    ///              @param vError			The vector containing the errors of the Y values. This vector will not be interpolated!
    ///             
    ///              @return	TRUE if successful, FALSE otherwise
    bool SetData(CVector& vXValues, CVector& vYValues) override
    {
        return IParamFunction::SetData(vXValues, vYValues) && this->InitializeSpline();
    }

    /// Returns the value of the cubic spline at the given X value.
    ///             
    ///              @param fXValue	The X value at which to evaluate the spline
    ///             
    ///              @return	The evaluated spline value.
    double GetValue(double fXValue) override
    {
        return this->EvaluateSpline(fXValue);
    }

    /// Calculates the function values at a set of given data points.
    ///             
    ///              @param vXValues			A vector object containing the X values at which the function has to be evaluated.
    ///              @param vYTargetVector	A vector object which receives the resulting function values.
    ///             
    ///              @return	A reference to the Y vector object
    CVector& GetValues(CVector& vXValues, CVector& vYTargetVector) override
    {
        if (vXValues.GetSize() != vYTargetVector.GetSize())
        {
            vYTargetVector.SetSize(vXValues.GetSize());
        }
        return this->EvaluateSplineVector(vXValues, vYTargetVector);
    }

    /// Returns the first derivative of the spline at the given data point.
    ///             
    ///              @param fXValue	The X value at which the slope is needed.
    ///             
    ///              @return	The slope of the B-Spline at the given data point.
    double GetSlope(double fXValue) override
    {
        return this->SlopeSpline(fXValue);
    }

    /// Calculates the first derivative of the function at a set of given data points.
    ///             
    ///              @param vXValues		A vector object containing the X values at which the function has to be evaluated.
    ///              @param vSlopeVector	A vector object which receives the resulting function values.
    ///             
    ///              @return	A reference to the slope vector object.
    CVector& GetSlopes(CVector& vXValues, CVector& vSlopeVector) override
    {
        if (vXValues.GetSize() != vSlopeVector.GetSize())
        {
            vSlopeVector.SetSize(vXValues.GetSize());
        }
        return this->SlopeSplineVector(vXValues, vSlopeVector);
    }

private:

    bool InitializeSpline()
    {
        int size = this->mXData.GetSize();
        if (size <= 3)
            return false;
        if (this->mXData[0] > this->mXData[this->mXData.GetSize() - 1])
        {
            this->mXData.Reverse();
            this->mYData.Reverse();
            this->mError.Reverse();
        }
        this->mY2ndDerivates.SetSize(size);
        CVector vector(size - 1);
        this->mY2ndDerivates[0] = 0.0;
        this->mY2ndDerivates[size - 1] = 0.0;
        vector[0] = 0.0;
        double num1 = vector[0];
        double num2 = this->mXData[0];
        double num3 = this->mXData[1];
        double num4 = num3 - num2;
        double d1 = this->mYData[0];
        if (std::isnan(d1) || std::isinf(d1))
            d1 = 0.0;
        double d2 = this->mYData[1];
        if (std::isnan(d2) || std::isinf(d2))
            d2 = 0.0;
        double num5 = d2 - d1;
        double num6 = this->mY2ndDerivates[0];
        for (int index = 1; index < size - 1; ++index)
        {
            double num7 = num3;
            num3 = this->mXData[index + 1];
            double num8 = num4;
            num4 = num3 - num7;
            double num9 = num8 + num4;
            double num10 = num8 / num9;
            double num11 = num10 * num6 + 2.0;
            num6 = (num10 - 1.0) / num11;
            this->mY2ndDerivates[index] = num6;
            double num12 = d2;
            d2 = this->mYData[index + 1];
            if (std::isnan(d2) || std::isinf(d2))
                d2 = num12;
            double num13 = num5;
            num5 = d2 - num12;
            num1 = (6.0 * (num5 / num4 - num13 / num8) / num9 - num10 * num1) / num11;
            vector[index] = num1;
        }
        double num14 = this->mY2ndDerivates[size - 1];
        for (int index = size - 2; index >= 0; --index)
        {
            num14 = this->mY2ndDerivates[index] * num14 + vector[index];
            this->mY2ndDerivates[index] = num14;
        }
        this->mH.SetSize(size);
        this->mSlopeInvariant.SetSize(size);
        this->mDeltaHSquareLow.SetSize(size);
        this->mDeltaHSquareHigh.SetSize(size);
        this->mSlopeDeltaHSquareLow.SetSize(size);
        this->mSlopeDeltaHSquareHigh.SetSize(size);
        this->mH[0] = 0.0;
        this->mSlopeInvariant[0] = 0.0;
        this->mDeltaHSquareLow[0] = 0.0;
        this->mDeltaHSquareHigh[0] = 0.0;
        this->mSlopeDeltaHSquareLow[0] = 0.0;
        this->mSlopeDeltaHSquareHigh[0] = 0.0;
        for (int index = 1; index < size; ++index)
        {
            double num7 = this->mXData[index] - this->mXData[index - 1];
            this->mH[index] = num7;
            double num8 = num7 / 6.0;
            this->mSlopeInvariant[index] = (this->mYData[index] - this->mYData[index - 1]) / num7;
            double num9 = num7 * num7 / 6.0;
            this->mDeltaHSquareLow[index] = this->mY2ndDerivates[index - 1] * num9;
            this->mDeltaHSquareHigh[index] = this->mY2ndDerivates[index] * num9;
            this->mSlopeDeltaHSquareLow[index] = this->mY2ndDerivates[index - 1] * num8;
            this->mSlopeDeltaHSquareHigh[index] = this->mY2ndDerivates[index] * num8;
        }
        return true;
    }

    double EvaluateSpline(double fXValue)
    {
        MATHFIT_ASSERT(mY2ndDerivates.GetSize() >= 3);

        int index1 = this->mXData.FindIndex(fXValue, CVector::EIndexConditions::LESSEQUAL);
        if (index1 < 0)
            index1 = 0;
        int index2 = index1 + 1;
        if (index2 >= this->mXData.GetSize())
        {
            index2 = this->mXData.GetSize() - 1;
            index1 = index2 - 1;
        }
        double num1 = this->mXData[index2] - this->mXData[index1];
        double num2 = std::max(std::min((this->mXData[index2] - fXValue) / num1, 1.0), 0.0);
        double num3 = 1.0 - num2;
        return num2 * this->mYData[index1] + num3 * this->mYData[index2] + ((num2 * num2 * num2 - num2) * this->mY2ndDerivates[index1] + (num3 * num3 * num3 - num3) * this->mY2ndDerivates[index2]) * num1 * num1 / 6.0;
    }

    CVector& EvaluateSplineVector(CVector& vXData, CVector& vYData)
    {
        MATHFIT_ASSERT(mY2ndDerivates.GetSize() >= 3);

        bool flag = false;
        if (vXData[0] > vXData[vXData.GetSize() - 1])
        {
            flag = true;
            vXData.Reverse();
        }
        int index1 = this->mXData.FindIndex(vXData[0], CVector::EIndexConditions::LESSEQUAL);
        int num1 = this->mXData.GetSize() - 1;
        if (index1 < 0)
            index1 = 0;
        int index2 = index1 + 1;
        if (index2 > num1)
        {
            index2 = num1;
            index1 = index2 - 1;
        }
        int index3 = 0;
        int size = vXData.GetSize();
        do
        {
            double num2 = this->mXData[index2];
            double num3 = this->mYData[index1];
            double num4 = this->mYData[index2];
            double num5 = this->mH[index2];
            double num6 = this->mDeltaHSquareLow[index2];
            double num7 = this->mDeltaHSquareHigh[index2];
            double num8;
            for (num8 = vXData[index3]; num8 < num2 || index2 >= num1; num8 = vXData[index3])
            {
                double num9 = std::max(std::min((num2 - num8) / num5, 1.0), 0.0);
                double num10 = 1.0 - num9;
                double num11 = num9 * num3 + num10 * num4 + ((num9 * num9 * num9 - num9) * num6 + (num10 * num10 * num10 - num10) * num7);
                vYData[index3++] = num11;
                if (index3 >= size)
                {
                    if (flag)
                    {
                        vXData.Reverse();
                        vYData.Reverse();
                    }
                    return vYData;
                }
            }
            do
            {
                ++index2;
                if (index2 > num1)
                {
                    index2 = num1;
                    break;
                }
            } while (this->mXData[index2] < num8);
            index1 = index2 - 1;
        } while (index3 < size);
        if (flag)
        {
            vXData.Reverse();
            vYData.Reverse();
        }
        return vYData;
    }

    double SlopeSpline(double fXValue)
    {
        // check wheter we have a valid spline
        MATHFIT_ASSERT(mY2ndDerivates.GetSize() >= 3);

        int index1 = this->mXData.FindIndex(fXValue, CVector::EIndexConditions::LESSEQUAL);
        if (index1 < 0)
            index1 = 0;
        int index2 = index1 + 1;
        if (index2 >= this->mXData.GetSize())
        {
            index2 = this->mXData.GetSize() - 1;
            index1 = index2 - 1;
        }
        double num1 = this->mXData[index2] - this->mXData[index1];
        double num2 = std::max(std::min((this->mXData[index2] - fXValue) / num1, 1.0), 0.0);
        double num3 = 1.0 - num2;
        return (this->mYData[index2] - this->mYData[index1]) / num1 + ((3.0 * num3 * num3 - 1.0) * this->mY2ndDerivates[index2] - (3.0 * num2 * num2 - 1.0) * this->mY2ndDerivates[index1]) * num1 / 6.0;
    }

    CVector& SlopeSplineVector(CVector& vXData, CVector& vYData)
    {
        MATHFIT_ASSERT(mY2ndDerivates.GetSize() >= 3);

        bool flag = false;
        if (vXData[0] > vXData[vXData.GetSize() - 1])
        {
            flag = true;
            vXData.Reverse();
        }
        int num1 = this->mXData.FindIndex(vXData[0], CVector::EIndexConditions::LESSEQUAL);
        int num2 = this->mXData.GetSize() - 1;
        if (num1 < 0)
            num1 = 0;
        int index1 = num1 + 1;
        int num3;
        if (index1 > num2)
        {
            index1 = num2;
            num3 = index1 - 1;
        }
        int index2 = 0;
        int size = vXData.GetSize();
        do
        {
            double num4 = this->mXData[index1];
            double num5 = this->mH[index1];
            double num6 = this->mSlopeInvariant[index1];
            double num7 = this->mSlopeDeltaHSquareLow[index1];
            double num8 = this->mSlopeDeltaHSquareHigh[index1];
            double num9;
            for (num9 = vXData[index2]; num9 < num4 || index1 >= num2; num9 = vXData[index2])
            {
                double num10 = std::max(std::min((num4 - num9) / num5, 1.0), 0.0);
                double num11 = 1.0 - num10;
                double num12 = num6 + ((3.0 * num11 * num11 - 1.0) * num8 - (3.0 * num10 * num10 - 1.0) * num7);
                vYData[index2++] = num12;
                if (index2 >= size)
                {
                    if (flag)
                    {
                        vXData.Reverse();
                        vYData.Reverse();
                    }
                    return vYData;
                }
            }
            do
            {
                ++index1;
                if (index1 > num2)
                {
                    index1 = num2;
                    break;
                }
            } while (this->mXData[index1] < num9);
            num3 = index1 - 1;
        } while (index2 < size);
        if (flag)
        {
            vXData.Reverse();
            vYData.Reverse();
        }
        return vYData;
    }

    /// <summary>
    /// Required for the ParamFunction interface. Throws NotImplemented exception always
    /// </summary>
    double GetLinearBasisFunction(double /*fXValue*/, int /*iParamID*/, bool /*bFixedID*/) override
    {
        throw std::exception("The method GetLinearBasisFunction in CCubicSplineParamFunction is not implemented");
    }
};

}

