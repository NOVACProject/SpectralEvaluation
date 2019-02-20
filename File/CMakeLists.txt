# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)

set(SPECTRUM_FILE_HEADERS
    ${CMAKE_CURRENT_LIST_DIR}/File.h
    ${CMAKE_CURRENT_LIST_DIR}/MKPack.h
    ${CMAKE_CURRENT_LIST_DIR}/ScanFileHandler.h
    ${CMAKE_CURRENT_LIST_DIR}/SpectrumIO.h
    PARENT_SCOPE)


set(SPECTRUM_FILE_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/File.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MKPack.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ScanFileHandler.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SpectrumIO.cpp
    PARENT_SCOPE)
