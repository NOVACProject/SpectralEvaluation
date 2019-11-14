# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)

set(SPECTRUM_EVALUATION_HEADERS
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Evaluation/BasicMath.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Evaluation/DarkSpectrum.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Evaluation/EvaluationBase.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Evaluation/EvaluationResult.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Evaluation/FitParameter.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Evaluation/FitWindow.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Evaluation/CrossSectionData.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Evaluation/PlumeSpectrumSelector.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Evaluation/Ratio.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Evaluation/RatioEvaluation.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Evaluation/ReferenceFile.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Evaluation/ReferenceFitResult.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Evaluation/ScanEvaluationBase.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Evaluation/WavelengthFit.h
    PARENT_SCOPE)


set(SPECTRUM_EVALUATION_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/BasicMath.cpp
    ${CMAKE_CURRENT_LIST_DIR}/BasicScanEvaluationResult.cpp
    ${CMAKE_CURRENT_LIST_DIR}/DarkSpectrum.cpp
    ${CMAKE_CURRENT_LIST_DIR}/EvaluationBase.cpp
    ${CMAKE_CURRENT_LIST_DIR}/EvaluationResult.cpp
    ${CMAKE_CURRENT_LIST_DIR}/FitWindow.cpp
    ${CMAKE_CURRENT_LIST_DIR}/CrossSectionData.cpp
    ${CMAKE_CURRENT_LIST_DIR}/PlumeSpectrumSelector.cpp
    ${CMAKE_CURRENT_LIST_DIR}/RatioEvaluation.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ReferenceFile.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ReferenceFitResult.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ScanEvaluationBase.cpp
    ${CMAKE_CURRENT_LIST_DIR}/WavelengthFit.cpp
    PARENT_SCOPE)
