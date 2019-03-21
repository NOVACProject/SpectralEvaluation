# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)

set(SPECTRUM_CLASS_HEADERS
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Air.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/DateTime.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/GpsData.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Spectra/Grid.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Spectra/ReferenceSpectrumConvolution.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Spectra/Spectrum.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Spectra/SpectrumInfo.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Spectra/SpectrometerModel.h
    PARENT_SCOPE)


set(SPECTRUM_CLASS_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/../DateTime.cpp
    ${CMAKE_CURRENT_LIST_DIR}/../GpsData.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ReferenceSpectrumConvolution.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Spectrum.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SpectrumInfo.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SpectrometerModel.cpp
    PARENT_SCOPE)
