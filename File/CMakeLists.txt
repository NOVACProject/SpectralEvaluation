# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)

set(SPECTRUM_FILE_HEADERS
    ${CMAKE_CURRENT_LIST_DIR}/File.h
    PARENT_SCOPE)


set(SPECTRUM_FILE_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/File.cpp
    PARENT_SCOPE)
