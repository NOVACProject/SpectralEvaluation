# This is the main project of Novac Post Processing Program (NovacPPP)

cmake_minimum_required (VERSION 3.6)

set(PppLib_INCLUDE_DIRS ${CMAKE_CURRENT_LIST_DIR}/include CACHE STRING "PPPLib headers")

# Add the different components
add_library(PPPLib
    ${PppLib_INCLUDE_DIRS}/PPPLib/CArray.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/CCriticalSection.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/CFileUtils.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/CFtpUtils.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/CList.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/CSingleLock.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/CStdioFile.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/CString.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/CStringTokenizer.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/Measurement.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/PPPLib.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/ThreadUtils.h
    ${PppLib_INCLUDE_DIRS}/PPPLib/VolcanoInfo.h
    ${PPPLIB_SPECTRA_HEADERS}
    ${SPECTRUM_FIT_HEADERS}
    ${SPECTRUM_EVALUATION_HEADERS}
    ${SPECTRUM_CLASS_HEADERS}
    ${SPECTRUM_FLUX_HEADERS}
    ${SPECTRUM_FILE_HEADERS}
    ${SPECTRUM_FILE_SOURCES}
    ${SPECTRUM_CONFIGURATION_HEADERS}
    ${SPECTRUM_CONFIGURATION_SOURCES}

    ${CMAKE_CURRENT_LIST_DIR}/src/CFileUtils.cpp 
    ${CMAKE_CURRENT_LIST_DIR}/src/CFtpUtils.cpp 
    ${CMAKE_CURRENT_LIST_DIR}/src/CStdioFile.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/CString.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/CStringTokenizer.cpp
    ${CMAKE_CURRENT_LIST_DIR}/src/VolcanoInfo.cpp
    ${PPPLIB_SPECTRA_SOURCES}
    ${SPECTRUM_EVALUATION_SOURCES}
    ${SPECTRUM_CLASS_SOURCES}
    ${SPECTRUM_FLUX_SOURCES}
)

target_include_directories (PPPLib PUBLIC ${PppLib_INCLUDE_DIRS} ${SPECTRALEVAUATION_INCLUDE_DIRS})
target_link_libraries(PPPLib PUBLIC NovacSpectralEvaluation)

IF(WIN32)
    target_compile_options(PPPLib PRIVATE /W4 /WX /sdl)
    target_compile_definitions(PPPLib PRIVATE _CRT_SECURE_NO_WARNINGS)
ELSE()
    target_compile_options(PPPLib PRIVATE -Wall -std=c++14)
ENDIF()

source_group("Header Files\\Configuration"          FILES ${SPECTRUM_CONFIGURATION_HEADERS})
source_group("Header Files\\Spectra"                FILES ${PPPLIB_SPECTRA_HEADERS} ${SPECTRUM_CLASS_HEADERS})
source_group("Header Files\\SpectrumFormat"         FILES ${})
source_group("Header Files\\SpectrumEvaluation"     FILES ${SPECTRUM_FIT_HEADERS})
source_group("Header Files\\Evaluation"             FILES ${SPECTRUM_EVALUATION_HEADERS})
source_group("Header Files\\Flux"                   FILES ${SPECTRUM_FLUX_HEADERS})
source_group("Header Files\\File"                   FILES ${SPECTRUM_FILE_HEADERS})


source_group("Source Files\\Configuration"          FILES ${SPECTRUM_CONFIGURATION_SOURCES})
source_group("Source Files\\Spectra"                FILES ${PPPLIB_SPECTRA_SOURCES} ${SPECTRUM_CLASS_SOURCES})
source_group("Source Files\\SpectrumFormat"         FILES ${})
source_group("Source Files\\SpectrumEvaluation"     FILES ${SPECTRUM_FIT_SOURCES})
source_group("Source Files\\Evaluation"             FILES ${SPECTRUM_EVALUATION_SOURCES})
source_group("Source Files\\Flux"                   FILES ${SPECTRUM_FLUX_SOURCES})
source_group("Source Files\\File"                   FILES ${SPECTRUM_FILE_SOURCES})
