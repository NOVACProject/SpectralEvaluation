# This is the main project of Novac Post Processing Program (NovacPPP)

cmake_minimum_required (VERSION 3.6)

## -------------------- SpectralEvaluationTests -------------------------

add_executable(SpectralEvaluationTests
    ${CMAKE_CURRENT_LIST_DIR}/main.cpp
    ${CMAKE_CURRENT_LIST_DIR}/catch.hpp
    ${CMAKE_CURRENT_LIST_DIR}/TestData.h
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_Correspondence.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_DarkSpectrum.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_DoasFit.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_EvaluationBase.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_File.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_FitWindowFileHandler.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_InstrumentCalibrationInStdFile.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_InstrumentLineshapeCalibrationController.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_InstrumentLineShapeEstimationFromDoas.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_InstrumentLineShapeEstimationFromKeypointDistance.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_PlumeSpectrumSelector.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_RatioEvaluation.cpp
    ${CMAKE_CURRENT_LIST_DIR}/IntegrationTests_RatioCalculationController.cpp
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
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_EvaluationResult.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_Fft.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_File.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_FitExtensions.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_FitWindow.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_FitInstrumentLineShapeFunction.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_GpsData.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_GuardedList.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_InstrumentLineShape.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_InstrumentLineShapeEstimation.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_InstrumentLineshapeEstimationFromDoas.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_Interpolation.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_PolynomialFit.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_PlumeInScanProperty.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_SpectrometerModel.cpp
    ${CMAKE_CURRENT_LIST_DIR}/UnitTest_SpectrumMath.cpp
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
    target_compile_options(SpectralEvaluationTests PRIVATE /W4 /WX /sdl /MP)
ELSE()
    # Make sure both the Windows and the Linux versions have binaries which end up in the same directory
    set_target_properties(SpectralEvaluationTests
        PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release"
    )

    target_compile_options(SpectralEvaluationTests PRIVATE -Wall -std=c++14 -fopenmp)
ENDIF()


