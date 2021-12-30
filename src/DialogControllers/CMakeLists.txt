# Adding the DialogControllers

cmake_minimum_required (VERSION 3.6)

set(DIALOGCONTROLLERS_HEADERS
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/DialogControllers/InstrumentCalibrationInput.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/DialogControllers/InstrumentLineshapeCalibrationController.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/DialogControllers/NovacProgramWavelengthCalibrationController.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/DialogControllers/ReferenceCreationController.h
    ${SPECTRALEVAUATION_INCLUDE_DIRS}/SpectralEvaluation/DialogControllers/WavelengthCalibrationController.h
    PARENT_SCOPE)

set(DIALOGCONTROLLERS_SOURCES
    ${CMAKE_CURRENT_LIST_DIR}/InstrumentLineshapeCalibrationController.cpp
    ${CMAKE_CURRENT_LIST_DIR}/NovacProgramWavelengthCalibrationController.cpp
    ${CMAKE_CURRENT_LIST_DIR}/ReferenceCreationController.cpp
    ${CMAKE_CURRENT_LIST_DIR}/WavelengthCalibrationController.cpp
    PARENT_SCOPE)
