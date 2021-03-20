#pragma once

#include <string>
#include <vector>

namespace novac
{

class CSpectrum;
class CCrossSectionData;

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

}
