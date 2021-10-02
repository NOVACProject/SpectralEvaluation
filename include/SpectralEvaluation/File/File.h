#pragma once

#include <string>
#include <vector>

// This file contains basic file input and output operations used throughout the novac software

namespace novac
{

class CSpectrum;
class CCrossSectionData;
class InstrumentCalibration;

/** Checks if the provided file name is an existing and readable file
    by attempting to open the file for reading. */
bool IsExistingFile(const std::string& fullFileName);

/** Reads cross section data from the provded file.
        This is able to read ASCII data files with data in one or two columns.
        If the file contains two columns then the first column will be saved in result.m_waveLength
            and the second column will be saved in result.m_crossSection.
        If the file only contains one column then if:
            1) saveAsWavelength == true: the column data will be saved in result.m_waveLength.
            2) saveAsWavelength == false: the column data will be saved in result.m_crossSection.
        If the file contains more than two columns, then the file read will fail.
        @param fullFilePath The full path to the file to read.  */
bool ReadCrossSectionFile(const std::string& fullFilePath, CCrossSectionData& result, bool saveAsWavelength = false);

/** Saves the provied cross section data to a space-separated file */
bool SaveCrossSectionFile(const std::string& fullFilePath, const CCrossSectionData& data);

/** Saves the provied spectrum to a space-separated file in cross-section file format */
bool SaveCrossSectionFile(const std::string& fullFilePath, const CSpectrum& data);

/** Saves the provided vector with values to a file, with one data point per line. */
bool SaveDataToFile(const std::string& fullFilePath, const std::vector<double>& data);

/** Attempts to read a CSpectrum from file. This can read .std or .txt data files */
bool ReadSpectrum(const std::string& fullFilePath, CSpectrum& result);

/** This is a helper method which makes sure that the provided filename contains a file suffix.
    If no suffix exists, then the default suffix is appended.
    NOTICE: the provided suffix should not start with a period. */
std::string EnsureFilenameHasSuffix(const std::string& fullFilePath, const std::string& defaultSuffix);

/** Returns the file extension (suffix) from the provided full file path */
std::string GetFileExtension(const std::string& fullFilePath);

/** Saves the full instrument calibration to a single file using the extended STD-format.
    This requires that at least that 'pixelToWavelengthMapping', 'instrumentLineShape', and 'instrumentLineShapeGrid' are defined. */
bool SaveInstrumentCalibration(const std::string& fullFilePath, const InstrumentCalibration& calibration);

/** Reads a saved instrument line shape and pixel-to-wavelength mapping from an instrument calibration file in extended STD-format.
*   TODO: Try to replace with the overload below taking an InstrumentCalibration instead.
    @return true on success. */
bool ReadInstrumentCalibration(const std::string& fullFilePath, CSpectrum& instrumentLineShape, std::vector<double>& pixelToWavelengthMapping);

/** Reads an instrument calibration from a file in the extended STD-format.
    The read in data will be filled into 'result'.
    @return true on success. */
bool ReadInstrumentCalibration(const std::string& fullFilePath, InstrumentCalibration& result);

}
