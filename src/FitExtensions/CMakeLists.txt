# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)

set(SPECTRUM_FIT_EXTENSION_HEADERS
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/FitExtensions/AsymmetricGaussFunction.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/FitExtensions/SuperGaussFunction.h
    PARENT_SCOPE)
    

    