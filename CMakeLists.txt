# This is the main project of Novac Post Processing Program (NovacPPP)

cmake_minimum_required (VERSION 3.6)

# Output the executables to the '/bin' directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/bin)

project(NovacSpectralEvaluation)


## -------------------- NovacSpectralEvaluation -------------------------

set(SPECTRALEVAUATION_INCLUDE_DIRS ${CMAKE_CURRENT_LIST_DIR} CACHE STRING "")

# Create a library from the SpectralEvaluation classes
add_subdirectory(Evaluation)
add_subdirectory(Fit)
add_subdirectory(Flux)
add_subdirectory(Spectra)

add_library(NovacSpectralEvaluation
    ${SPECTRUM_EVALUATION_HEADERS}
    ${SPECTRUM_EVALUATION_SOURCES}
    ${SPECTRUM_FIT_HEADERS}
    ${SPECTRUM_FIT_SOURCES}
    ${SPECTRUM_FLUX_HEADERS}
    ${SPECTRUM_FLUX_SOURCES}
    ${SPECTRUM_CLASS_HEADERS}
    ${SPECTRUM_CLASS_SOURCES}
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Utils.h
)

#Organize the files into folders, for ease of use
source_group("Header Files\\Evaluation"             FILES ${SPECTRUM_EVALUATION_HEADERS})
source_group("Header Files\\Fit"                    FILES ${SPECTRUM_FIT_HEADERS})
source_group("Header Files\\Flux"                   FILES ${SPECTRUM_FLUX_HEADERS})
source_group("Header Files\\Spectra"                FILES ${SPECTRUM_CLASS_HEADERS})

source_group("Source Files\\Evaluation"             FILES ${SPECTRUM_EVALUATION_SOURCES})
source_group("Source Files\\Fit"                    FILES ${SPECTRUM_FIT_SOURCES})
source_group("Source Files\\Flux"                   FILES ${SPECTRUM_FLUX_SOURCES})
source_group("Source Files\\Spectra"                FILES ${SPECTRUM_CLASS_SOURCES})


## -------------------- SpectralEvaluationTests -------------------------

add_subdirectory(UnitTests)

IF(MSVC)
    target_compile_definitions(NovacSpectralEvaluation PRIVATE -D_CRT_SECURE_NO_WARNINGS)
    target_compile_options(NovacSpectralEvaluation PRIVATE /W4 /WX /sdl)
ELSE()
    target_compile_options(NovacSpectralEvaluation PRIVATE -Wall -std=c++11)
ENDIF()


