# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)

set(SPECTRUM_CLASS_HEADERS
    ${CMAKE_CURRENT_LIST_DIR}/DateTime.h
    ${CMAKE_CURRENT_LIST_DIR}/GpsData.h
    ${CMAKE_CURRENT_LIST_DIR}/Spectrum.h
    ${CMAKE_CURRENT_LIST_DIR}/SpectrumInfo.h
    ${CMAKE_CURRENT_LIST_DIR}/SpectrometerModel.h
    PARENT_SCOPE)


set(SPECTRUM_CLASS_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/DateTime.cpp
    ${CMAKE_CURRENT_LIST_DIR}/GpsData.cpp
    ${CMAKE_CURRENT_LIST_DIR}/Spectrum.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SpectrumInfo.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SpectrometerModel.cpp
    PARENT_SCOPE)
