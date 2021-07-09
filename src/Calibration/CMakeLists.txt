# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)

set(SPECTRUM_CALIBRATION_HEADERS
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Calibration/InstrumentLineShape.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Calibration/ReferenceSpectrumConvolution.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Calibration/WavelengthCalibration.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Calibration/WavelengthCalibrationByRansac.h
    PARENT_SCOPE)


set(SPECTRUM_CALIBRATION_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/InstrumentLineShape.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ReferenceSpectrumConvolution.cpp
    ${CMAKE_CURRENT_LIST_DIR}/WavelengthCalibration.cpp
    ${CMAKE_CURRENT_LIST_DIR}/WavelengthCalibrationByRansac.cpp
    PARENT_SCOPE)
