# This is the main project of Novac SpectralEvaluation library which is an extract
#  of the shared classes and function used by the Novac software.

cmake_minimum_required (VERSION 3.6)

get_directory_property(hasParent PARENT_DIRECTORY)

if(NOT hasParent)
    # Output the executables to the '/bin' directory
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/bin)
    project(NovacSpectralEvaluation)
endif()

FIND_PACKAGE(OpenMP)
IF(OPENMP_FOUND)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

## -------------------- NovacSpectralEvaluation -------------------------

set(SPECTRALEVAUATION_INCLUDE_DIRS ${CMAKE_CURRENT_LIST_DIR}/include CACHE STRING "")

# Create a library from the SpectralEvaluation classes
add_subdirectory(src/Calibration)
add_subdirectory(src/Configuration)
add_subdirectory(src/DialogControllers)
add_subdirectory(src/Evaluation)
add_subdirectory(src/Fit)
add_subdirectory(src/FitExtensions)
add_subdirectory(src/File)
add_subdirectory(src/Flux)
add_subdirectory(src/Math)
add_subdirectory(src/Spectra)

add_library(NovacSpectralEvaluation
    ${SPECTRUM_CALIBRATION_HEADERS}
    ${SPECTRUM_CALIBRATION_SOURCES}
    ${SPECTRUM_EVALUATION_HEADERS}
    ${SPECTRUM_EVALUATION_SOURCES}
    ${SPECTRUM_FIT_HEADERS}
    ${SPECTRUM_FIT_SOURCES}
    ${SPECTRUM_FIT_EXTENSION_HEADERS}
    ${SPECTRUM_FLUX_HEADERS}
    ${SPECTRUM_FLUX_SOURCES}
    ${SPECTRUM_CLASS_HEADERS}
    ${SPECTRUM_CLASS_SOURCES}
    ${SPECTRUM_FILE_HEADERS}
    ${SPECTRUM_FILE_SOURCES}
    ${SPECTRUM_MATH_HEADERS}
    ${SPECTRUM_MATH_SOURCES}
    ${SPECTRUM_CONFIGURATION_HEADERS}
    ${SPECTRUM_CONFIGURATION_SOURCES}
    ${DIALOGCONTROLLERS_HEADERS}
    ${DIALOGCONTROLLERS_SOURCES}
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Air.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/DateTime.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Definitions.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/EstimatedValue.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Geometry.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/GPSData.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Interpolation.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Log.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Statistics.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/StringUtils.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/Units.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/VectorUtils.h
    ${CMAKE_CURRENT_SOURCE_DIR}/src/DateTime.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/EstimatedValue.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Geometry.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/GPSData.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Interpolation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Log.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/Statistics.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/StringUtils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/src/VectorUtils.cpp
)

target_include_directories(NovacSpectralEvaluation PUBLIC ${SPECTRALEVAUATION_INCLUDE_DIRS})

#Organize the files into folders, for ease of use
source_group("Header Files\\Calibration"            FILES ${SPECTRUM_CALIBRATION_HEADERS})
source_group("Header Files\\DialogControls"         FILES ${DIALOGCONTROLLERS_HEADERS})
source_group("Header Files\\Configuration"          FILES ${SPECTRUM_CONFIGURATION_HEADERS})
source_group("Header Files\\Evaluation"             FILES ${SPECTRUM_EVALUATION_HEADERS})
source_group("Header Files\\File"                   FILES ${SPECTRUM_FILE_HEADERS})
source_group("Header Files\\Fit"                    FILES ${SPECTRUM_FIT_HEADERS})
source_group("Header Files\\Fit"                    FILES ${SPECTRUM_FIT_EXTENSION_HEADERS})
source_group("Header Files\\Flux"                   FILES ${SPECTRUM_FLUX_HEADERS})
source_group("Header Files\\Math"                   FILES ${SPECTRUM_MATH_HEADERS})
source_group("Header Files\\Spectra"                FILES ${SPECTRUM_CLASS_HEADERS})


source_group("Source Files\\Calibration"            FILES ${SPECTRUM_CALIBRATION_SOURCES})
source_group("Source Files\\DialogControls"         FILES ${DIALOGCONTROLLERS_SOURCES})
source_group("Source Files\\Configuration"          FILES ${SPECTRUM_CONFIGURATION_SOURCES})
source_group("Source Files\\Evaluation"             FILES ${SPECTRUM_EVALUATION_SOURCES})
source_group("Source Files\\File"                   FILES ${SPECTRUM_FILE_SOURCES})
source_group("Source Files\\Fit"                    FILES ${SPECTRUM_FIT_SOURCES})
source_group("Source Files\\Flux"                   FILES ${SPECTRUM_FLUX_SOURCES})
source_group("Source Files\\Math"                   FILES ${SPECTRUM_MATH_SOURCES})
source_group("Source Files\\Spectra"                FILES ${SPECTRUM_CLASS_SOURCES})


## -------------------- SpectralEvaluationTests -------------------------

add_subdirectory(UnitTests)

IF(MSVC)
    target_compile_definitions(NovacSpectralEvaluation PRIVATE -D_CRT_SECURE_NO_WARNINGS)
    target_compile_options(NovacSpectralEvaluation PRIVATE /W4 /WX /sdl)
ELSE()
    target_compile_options(NovacSpectralEvaluation PRIVATE -Wall -std=c++14 -fopenmp)
ENDIF()


