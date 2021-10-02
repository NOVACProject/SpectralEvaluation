# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)

set(SPECTRUM_FLUX_HEADERS
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Flux/Flux.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Flux/PlumeInScanProperty.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Flux/ScanFluxResult.h
    PARENT_SCOPE)


set(SPECTRUM_FLUX_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/Flux.cpp
    ${CMAKE_CURRENT_LIST_DIR}/PlumeInScanProperty.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ScanFluxResult.cpp
    PARENT_SCOPE)
