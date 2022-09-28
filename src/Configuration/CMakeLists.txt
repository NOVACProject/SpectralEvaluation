# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)

set(SPECTRUM_CONFIGURATION_HEADERS
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Configuration/DarkSettings.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Configuration/RatioEvaluationSettings.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Configuration/SkySettings.h
    PARENT_SCOPE)


set(SPECTRUM_CONFIGURATION_SOURCES
    PARENT_SCOPE)
