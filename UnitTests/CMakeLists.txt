# This is the main project of Novac Post Processing Program (NovacPPP)

cmake_minimum_required (VERSION 3.6)

## -------------------- SpectralEvaluationTests -------------------------

add_executable(SpectralEvaluationTests
    ${CMAKE_CURRENT_LIST_DIR}/main.cpp
    ${CMAKE_CURRENT_LIST_DIR}/catch.hpp
    ${CMAKE_CURRENT_LIST_DIR}/TestData.h
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_Correspondence.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_File.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_InstrumentCalibrationInStdFile.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_InstrumentLineshapeCalibrationController.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_InstrumentLineShapeEstimationFromDoas.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_InstrumentLineShapeEstimationFromKeypointDistance.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_PlumeSpectrumSelector.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_ScanEvaluationLogFileHandler.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_ScanFileHandler.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_SpectrumIO.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_StdFile.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_WavelengthCalibration.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_WavelengthCalibrationController.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_Air.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_Convolution.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_Correspondence.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_CrossSectionData.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_CVector.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_DateTime.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_Evaluation.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_Fft.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_FindPlume.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_FitExtensions.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_FitWindow.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_FitInstrumentLineShapeFunction.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_InstrumentLineShape.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_InstrumentLineShapeEstimation.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_InstrumentLineshapeEstimationFromDoas.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_Interpolation.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_PolynomialFit.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_RatioEvaluation.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_SpectrometerModel.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_SpectrumUtils.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_Utils.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_VectorUtils.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_WavelengthCalibration.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_XmlUtils.cpp
)

target_include_directories(SpectralEvaluationTests PRIVATE ${SPECTRALEVAUATION_INCLUDE_DIRS})
target_link_libraries(SpectralEvaluationTests PRIVATE NovacSpectralEvaluation)

IF(MSVC)
    target_compile_definitions(SpectralEvaluationTests PRIVATE -D_CRT_SECURE_NO_WARNINGS)
    target_compile_options(SpectralEvaluationTests PRIVATE /W4 /WX /sdl)
ELSE()
    target_compile_options(SpectralEvaluationTests PRIVATE -Wall -std=c++14)
ENDIF()


