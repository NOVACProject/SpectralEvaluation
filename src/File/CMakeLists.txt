# Adding the Common helper files

cmake_minimum_required (VERSION 3.6)

set(SPECTRUM_FILE_HEADERS
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/File/File.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/File/MKPack.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/File/ScanFileHandler.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/File/ScanEvaluationLogFileHandler.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/File/SpectrumIO.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/File/STDFile.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/File/TXTFile.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/File/XmlUtil.h
    PARENT_SCOPE)


set(SPECTRUM_FILE_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/File.cpp
    ${CMAKE_CURRENT_LIST_DIR}/MKPack.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ScanFileHandler.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ScanEvaluationLogFileHandler.cpp
    ${CMAKE_CURRENT_LIST_DIR}/SpectrumIO.cpp
    ${CMAKE_CURRENT_LIST_DIR}/STDFile.cpp
    ${CMAKE_CURRENT_LIST_DIR}/TXTFile.cpp
    ${CMAKE_CURRENT_LIST_DIR}/XmlUtil.cpp
    PARENT_SCOPE)
