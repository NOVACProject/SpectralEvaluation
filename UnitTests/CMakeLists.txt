# This is the main project of Novac Post Processing Program (NovacPPP)

cmake_minimum_required (VERSION 3.6)

## -------------------- SpectralEvaluationTests -------------------------

add_executable(SpectralEvaluationTests
    ${CMAKE_CURRENT_LIST_DIR}/main.cpp
    ${CMAKE_CURRENT_LIST_DIR}/catch.hpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTest_ScanFileHandler.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTest_SpectrumIO.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_Air.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_Convolution.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_CrossSectionData.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_DateTime.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_Evaluation.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_FindPlume.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_FitWindow.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_RatioEvaluation.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_Utils.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_VectorUtils.cpp
)

target_include_directories(SpectralEvaluationTests PRIVATE ${SPECTRALEVAUATION_INCLUDE_DIRS})
target_link_libraries(SpectralEvaluationTests PRIVATE NovacSpectralEvaluation)

IF(MSVC)
    target_compile_definitions(SpectralEvaluationTests PRIVATE -D_CRT_SECURE_NO_WARNINGS)
    target_compile_options(SpectralEvaluationTests PRIVATE /W4 /WX /sdl)
ELSE()
    target_compile_options(SpectralEvaluationTests PRIVATE -Wall -std=c++11)
ENDIF()


