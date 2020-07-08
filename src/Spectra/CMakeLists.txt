# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)

set(SPECTRUM_CLASS_HEADERS
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Air.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/DateTime.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/GPSData.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/EstimatedValue.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Spectra/Grid.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Spectra/Scattering.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Spectra/Spectrum.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Spectra/SpectrumInfo.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Spectra/SpectrometerModel.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Spectra/SpectrumUtils.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Spectra/WavelengthRange.h
    PARENT_SCOPE)


set(SPECTRUM_CLASS_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/../DateTime.cpp
    ${CMAKE_CURRENT_LIST_DIR}/../GPSData.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Scattering.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Spectrum.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SpectrumInfo.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SpectrometerModel.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SpectrumUtils.cpp
    PARENT_SCOPE)
