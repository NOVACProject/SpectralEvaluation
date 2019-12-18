#pragma once

#include <string>

namespace Evaluation
{
    class CCrossSectionData;
}

class CSpectrum;

namespace FileIo
{
    /** Reads cross section data from the provded file.
        This is able to read ASCII data files with data in one or two columns.
        If the file contains two columns then the first column will be saved in result.m_waveLength
            and the second column will be saved in result.m_crossSection.
        If the file only contains one column then if:
            1) saveAsWavelength == true: the column data will be saved in result.m_waveLength.
            2) saveAsWavelength == false: the column data will be saved in result.m_crossSection.
        If the file contains more than two columns, then the file read will fail.
        @param fullFilePath The full path to the file to read.  */
    bool ReadCrossSectionFile(const std::string& fullFilePath, Evaluation::CCrossSectionData& result, bool saveAsWavelength = false);

    /** Saves the provied cross section data to a space-separated file */
    bool SaveCrossSectionFile(const std::string& fullFilePath, const Evaluation::CCrossSectionData& data);

    /** Attempts to read a CSpectrum from file. This can read .std or .txt data files */
    bool ReadSpectrum(const std::string& fullFilePath, CSpectrum& result);
}
