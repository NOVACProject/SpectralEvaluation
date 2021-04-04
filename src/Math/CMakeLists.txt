# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)

set(SPECTRUM_MATH_HEADERS
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Math/FFT.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Math/FunctionFit.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Math/PolynomialFit.h
    PARENT_SCOPE)


set(SPECTRUM_MATH_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/FFT.cpp
    ${CMAKE_CURRENT_LIST_DIR}/FunctionFit.cpp
    ${CMAKE_CURRENT_LIST_DIR}/PolynomialFit.cpp
    PARENT_SCOPE)
