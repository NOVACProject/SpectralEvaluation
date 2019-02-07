# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)

set(SPECTRUM_FIT_HEADERS
    ${CMAKE_CURRENT_LIST_DIR}/ApertureFunction.h
    ${CMAKE_CURRENT_LIST_DIR}/BandedMatrix.h
    ${CMAKE_CURRENT_LIST_DIR}/BSpline.h
    ${CMAKE_CURRENT_LIST_DIR}/BSplineD.h
    ${CMAKE_CURRENT_LIST_DIR}/BSplineF.h
    ${CMAKE_CURRENT_LIST_DIR}/BSplineImpl.h
    ${CMAKE_CURRENT_LIST_DIR}/ConstFunction.h
    ${CMAKE_CURRENT_LIST_DIR}/ConvoluteFunction.h
    ${CMAKE_CURRENT_LIST_DIR}/ConvolutionCoreFunction.h
    ${CMAKE_CURRENT_LIST_DIR}/CubicBSplineFunction.h
    ${CMAKE_CURRENT_LIST_DIR}/CubicSplineFunction.h
    ${CMAKE_CURRENT_LIST_DIR}/DataSet.h
    ${CMAKE_CURRENT_LIST_DIR}/DiscreteFunction.h
    ${CMAKE_CURRENT_LIST_DIR}/DivFunction.h
    ${CMAKE_CURRENT_LIST_DIR}/DOASVector.h
    ${CMAKE_CURRENT_LIST_DIR}/ExpFunction.h
    ${CMAKE_CURRENT_LIST_DIR}/Fit.h
    ${CMAKE_CURRENT_LIST_DIR}/FitBasic.h
    ${CMAKE_CURRENT_LIST_DIR}/FitException.h
    ${CMAKE_CURRENT_LIST_DIR}/Function.h
    ${CMAKE_CURRENT_LIST_DIR}/GaussFunction.h
    ${CMAKE_CURRENT_LIST_DIR}/LeastSquareFit.h
    ${CMAKE_CURRENT_LIST_DIR}/LevenbergMarquardtFit.h
    ${CMAKE_CURRENT_LIST_DIR}/LinearInterpolateFunction.h
    ${CMAKE_CURRENT_LIST_DIR}/LnFunction.h
    ${CMAKE_CURRENT_LIST_DIR}/Matrix.h
    ${CMAKE_CURRENT_LIST_DIR}/Minimizer.h
    ${CMAKE_CURRENT_LIST_DIR}/MessageLog.h
    ${CMAKE_CURRENT_LIST_DIR}/MulFunction.h
    ${CMAKE_CURRENT_LIST_DIR}/NegateFunction.h
    ${CMAKE_CURRENT_LIST_DIR}/NonlinearParameterFunction.h
    ${CMAKE_CURRENT_LIST_DIR}/ParameterLinkItem.h
    ${CMAKE_CURRENT_LIST_DIR}/ParameterVector.h
    ${CMAKE_CURRENT_LIST_DIR}/ParamFunction.h
    ${CMAKE_CURRENT_LIST_DIR}/PolynomialFunction.h
    ${CMAKE_CURRENT_LIST_DIR}/ReferenceSpectrumFunction.h
    ${CMAKE_CURRENT_LIST_DIR}/SimpleDOASFunction.h
    ${CMAKE_CURRENT_LIST_DIR}/StandardFit.h
    ${CMAKE_CURRENT_LIST_DIR}/StandardMetricFunction.h
    ${CMAKE_CURRENT_LIST_DIR}/StatisticVector.h
    ${CMAKE_CURRENT_LIST_DIR}/SumFunction.h
    ${CMAKE_CURRENT_LIST_DIR}/Vector.h
    PARENT_SCOPE)
    
    
set(SPECTRUM_FIT_SOURCES
    PARENT_SCOPE)

    