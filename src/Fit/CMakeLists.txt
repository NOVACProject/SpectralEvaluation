# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)

set(SPECTRUM_FIT_HEADERS
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/ApertureFunction.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/BandedMatrix.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/BSpline.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/BSplineD.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/BSplineF.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/BSplineImpl.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/ConstFunction.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/ConvoluteFunction.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/ConvolutionCoreFunction.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/CubicBSplineFunction.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/CubicSplineFunction.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/DataSet.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/DiscreteFunction.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/DivFunction.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/DOASVector.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/ExpFunction.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/Fit.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/FitBasic.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/FitException.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/Function.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/GaussFunction.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/LeastSquareFit.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/LevenbergMarquardtFit.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/LinearInterpolateFunction.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/LnFunction.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/Matrix.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/Minimizer.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/MessageLog.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/MulFunction.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/NegateFunction.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/NonlinearParameterFunction.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/ParameterLinkItem.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/ParameterVector.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/ParamFunction.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/PolynomialFunction.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/ReferenceSpectrumFunction.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/SimpleDOASFunction.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/StandardFit.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/StandardMetricFunction.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/StatisticVector.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/SumFunction.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Fit/Vector.h
    PARENT_SCOPE)
    
    
set(SPECTRUM_FIT_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/MessageLog.cpp
    PARENT_SCOPE)

    