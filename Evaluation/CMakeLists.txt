# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)

set(SPECTRUM_EVALUATION_HEADERS
    ${CMAKE_CURRENT_LIST_DIR}/BasicMath.h
    ${CMAKE_CURRENT_LIST_DIR}/EvaluationBase.h
    ${CMAKE_CURRENT_LIST_DIR}/EvaluationResult.h
    ${CMAKE_CURRENT_LIST_DIR}/FitParameter.h
    ${CMAKE_CURRENT_LIST_DIR}/FitWindow.h
    ${CMAKE_CURRENT_LIST_DIR}/CrossSectionData.h
    ${CMAKE_CURRENT_LIST_DIR}/RatioEvaluation.h
    ${CMAKE_CURRENT_LIST_DIR}/ReferenceFile.h
    ${CMAKE_CURRENT_LIST_DIR}/ReferenceFitResult.h
    ${CMAKE_CURRENT_LIST_DIR}/ScanEvaluationBase.h
    PARENT_SCOPE)


set(SPECTRUM_EVALUATION_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/BasicMath.cpp
    ${CMAKE_CURRENT_LIST_DIR}/EvaluationBase.cpp
    ${CMAKE_CURRENT_LIST_DIR}/EvaluationResult.cpp
    ${CMAKE_CURRENT_LIST_DIR}/FitWindow.cpp
    ${CMAKE_CURRENT_LIST_DIR}/CrossSectionData.cpp
    ${CMAKE_CURRENT_LIST_DIR}/RatioEvaluation.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ReferenceFile.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ReferenceFitResult.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ScanEvaluationBase.cpp
    PARENT_SCOPE)
