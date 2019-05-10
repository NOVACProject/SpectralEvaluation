/**
 * Contains a defintion of a function that makes all parameters nonlinear.
 *
 * @author		\URL[Stefan Kraus]{http://stefan@00kraus.de} @ \URL[IWR, Image Processing Group]{http://klimt.iwr.uni-heidelberg.de}
 * @version		1.0 @ 2002/05/02
 */
#if !defined(NONLINEARPARAMETERFUNCTION_H_020502)
#define NONLINEARPARAMETERFUNCTION_H_020502

#include "ParamFunction.h"

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#ifdef _MSC_VER
#pragma warning (push, 3)
#endif

namespace MathFit
{
	/**
	* This object represents a function that makes all parameters nonlinear.
	*
	* @author		\URL[Stefan Kraus]{http://stefan@00kraus.de} @ \URL[IWR, Image Processing Group]{http://klimt.iwr.uni-heidelberg.de}
	* @version		1.0 @ 2002/05/02
	*/
	class CNonlinearParameterFunction : public IParamFunction
	{
	public:
		/**
		* Creates the object and associates a function object.
		*
		* @param ipfOperand	The function object, which's parameters should be nonlinear.
		*/
		CNonlinearParameterFunction(IParamFunction& ipfOperand) : mOperand(ipfOperand)
		{
			BuildLinearParameter();
			BuildNonlinearParameter();
		}

		/**
		* Returns the value of the operand function object.
		*
		* @param fXValue	The X data point.
		*
		* @return	The value of the function at the given data point.
		*/
		virtual TFitData GetValue(TFitData fXValue)
		{
			return mOperand.GetValue(fXValue);
		}

		/**
		* Returns the first derivative of the function at the given data point.
		*
		* @param fXValue	The X data point at which the first derivation is needed.
		*
		* @return	The slope of the function at the given data point.
		*/
		virtual TFitData GetSlope(TFitData fXValue)
		{
			return mOperand.GetSlope(fXValue);
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
		virtual TFitData GetLinearBasisFunction(TFitData /*fXValue*/, int /*iParamID*/, bool /*bFixedID = true*/)
		{
			throw(EXCEPTION(CNotImplementedException));
		}

		/**
		* Sets the new set of nonlinear parameters.
		* The given parameter vector is split up into the pieces needed by the references.
		*
		* @param vParam	The vector containing the new unfixed parameters.
		*
		* @return	TRUE if successful, FALSE otherwise.
		*/
		virtual bool SetNonlinearParameter(CVector& vParam)
		{
			// keep internal copy
			mNonlinearParams.SetParameters(vParam);

			// split up vev�ctor into linear and nonlinear parts.
			const int iLinSize = mOperand.GetLinearParameter().GetSize();
			if(iLinSize > 0)
			{
				auto temp = vParam.SubVector(0, iLinSize);
				mOperand.SetLinearParameter(temp);
			}

			const int iNonlinSize = mOperand.GetNonlinearParameter().GetSize();
			if(iNonlinSize > 0)
			{
				auto temp = vParam.SubVector(iLinSize, iNonlinSize);
				mOperand.SetNonlinearParameter(temp);
			}

			return true;
		}

		/**
		* Sets the covariance matrix of the nonlinear parameters.
		*
		* @param mCovar	The covariance matrix.
		*/
		virtual void SetNonlinearCovarMatrix(CMatrix& mCovar)
		{
			mNonlinearParams.SetCovarMatrix(mCovar);

			// split up vev�ctor into linear and nonlinear parts.
			const int iLinSize = mOperand.GetLinearParameter().GetSize();
			if(iLinSize > 0)
			{
				auto temp = mCovar.SubMatrix(0, 0, iLinSize, iLinSize);
				mOperand.SetLinearCovarMatrix(temp);
			}

			const int iNonlinSize = mOperand.GetNonlinearParameter().GetSize();
			if(iNonlinSize > 0)
			{
				auto temp = mCovar.SubMatrix(iLinSize, iLinSize, iNonlinSize, iNonlinSize);
				mOperand.SetNonlinearCovarMatrix(temp);
			}
		}

		/**
		* Sets the correlation matrix of the nonlinear parameters.
		*
		* @param mCorrel	The correlation matrix.
		*/
		virtual void SetNonlinearCorrelMatrix(CMatrix& mCorrel)
		{
			mNonlinearParams.SetCorrelMatrix(mCorrel);

			// split up vev�ctor into linear and nonlinear parts.
			const int iLinSize = mOperand.GetLinearParameter().GetSize();
			if(iLinSize > 0)
			{
				auto temp = mCorrel.SubMatrix(0, 0, iLinSize, iLinSize);
				mOperand.SetLinearCorrelMatrix(temp);
			}

			const int iNonlinSize = mOperand.GetNonlinearParameter().GetSize();
			if(iNonlinSize > 0)
			{
				auto temp = mCorrel.SubMatrix(iLinSize, iLinSize, iNonlinSize, iNonlinSize);
				mOperand.SetNonlinearCorrelMatrix(temp);
			}
		}

		/**
		* Sets the error of the nonlinear parameters.
		*
		* @param vError	The error vector.
		*/
		virtual void SetNonlinearError(CVector& vError)
		{
			mNonlinearParams.SetError(vError);

			// split up vev�ctor into linear and nonlinear parts.
			const int iLinSize = mOperand.GetLinearParameter().GetSize();
			if(iLinSize > 0)
			{
				auto temp = vError.SubVector(0, iLinSize);
				mOperand.SetLinearError(temp);
			}

			const int iNonlinSize = mOperand.GetNonlinearParameter().GetSize();
			if(iNonlinSize > 0)
			{
				auto temp = vError.SubVector(iLinSize, iNonlinSize);
				mOperand.SetNonlinearError(temp);
			}
		}

		/**
		* Resets the linear parameters to default values.
		* The default implementation sets every parameter to zero.
		*/
		virtual void ResetLinearParameter()
		{
			mOperand.ResetLinearParameter();
			BuildLinearParameter();
		}

		/**
		* Resets the nonlinear parameters to default values.
		* The default implementation sets every parameter to zero.
		*/
		virtual void ResetNonlinearParameter()
		{
			mOperand.ResetNonlinearParameter();
			BuildNonlinearParameter();
		}

		virtual TFitData GetNonlinearPenalty(TFitData fChiSquare)
		{
			return mOperand.GetNonlinearPenalty(fChiSquare) + mOperand.GetLinearPenalty(fChiSquare);
		}

	private:
		/**
		* Builds the linear parameter vector from the operands parameters.
		*/
		void BuildLinearParameter()
		{
			// no linear parameters, since everything gets nonlinear
			mLinearParams.SetSize(0);
		}

		/**
		* Builds the nonlinear parameter vector from the operands parameters.
		*/
		void BuildNonlinearParameter()
		{
			// append the two parameter vectors
			CVector vBuf;

			vBuf.Append(mOperand.GetLinearParameter());
			vBuf.Append(mOperand.GetNonlinearParameter());
			mNonlinearParams.SetSize(vBuf.GetSize());
			mNonlinearParams.SetParameters(vBuf);
		}

		/**
		* Hold the logarithm's parameter function.
		*/
		IParamFunction& mOperand;
	};
}

#pragma warning (pop)
#endif