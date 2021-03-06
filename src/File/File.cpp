#include <SpectralEvaluation/StringUtils.h>
#include <SpectralEvaluation/File/File.h>
#include <SpectralEvaluation/File/STDFile.h>
#include <SpectralEvaluation/File/TXTFile.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <cstring>

namespace FileIo
{
    // @return the number of columns read.
    int ReadFromFile(FILE* f, const char* formatStr, double& col1, double& col2)
    {
        int nCols = fscanf(f, formatStr, &col1, &col2);
        if (nCols != 0)
        {
            return nCols;
        }
        return 0;
    }

    int GetFileFormat(const char* string, char* format)
    {
        if (nullptr == string || nullptr == format)
        {
            return 0;
        }
        if (0 == strlen(string))
        {
            return 0;
        }

        double col1 = 0, col2 = 0;

        if (2 == sscanf(string, "%lf %lf", &col1, &col2))
        {
            sprintf(format, "%s", "%lf %lf");
            return 2;
        }
        if (2 == sscanf(string, "%lf,%lf", &col1, &col2))
        {
            sprintf(format, "%s", "%lf,%lf");
            return 2;
        }
        if (2 == sscanf(string, "%lf;%lf", &col1, &col2))
        {
            sprintf(format, "%s", "%lf;%lf");
            return 2;
        }
        if (1 == sscanf(string, "%lf", &col1))
        {
            sprintf(format, "%s", "%lf");
            return 1;
        }

        return 0;
    }

    bool ReadCrossSectionFromStdFile(const std::string& fullFilePath, Evaluation::CCrossSectionData& result, bool saveAsWavelength)
    {
        CSpectrum spectrum;
        if (!SpectrumIO::CSTDFile::ReadSpectrum(spectrum, fullFilePath))
        {
            return false;
        }

        if (saveAsWavelength)
        {
            result.m_waveLength = std::vector<double>(spectrum.m_data, spectrum.m_data + spectrum.m_length);
        }
        else
        {
            result.m_crossSection = std::vector<double>(spectrum.m_data, spectrum.m_data + spectrum.m_length);
        }

        return true;
    }

    bool ReadCrossSectionFile(const std::string& fullFilePath, Evaluation::CCrossSectionData& result, bool saveAsWavelength)
    {
        const std::string fileEnding = Right(fullFilePath, 4);
        if (EqualsIgnoringCase(fileEnding, ".std"))
        {
            return ReadCrossSectionFromStdFile(fullFilePath, result, saveAsWavelength);
        }

        FILE* f = fopen(fullFilePath.c_str(), "r");
        if (nullptr == f)
        {
            return false;
        }

        // Get the file-format and number of columns from the file.
        char tempBuffer[8192];
        char format[256];
        if (nullptr == fgets(tempBuffer, 8191, f))
        {
            fclose(f);
            return false;
        }

        const int numColumns = GetFileFormat(tempBuffer, format);

        fseek(f, 0, SEEK_SET);

        // Make some space
        result.m_crossSection.reserve(2050);
        result.m_waveLength.reserve(2050);

        while (!feof(f))
        {
            double col1 = 0.0;
            double col2 = 0.0;
            int nCols = ReadFromFile(f, format, col1, col2);

            if (nCols != numColumns)
            {
                break;
            }

            if (nCols == 1)
            {
                if (saveAsWavelength)
                    result.m_waveLength.push_back(col1);
                else
                    result.m_crossSection.push_back(col1);
            }
            else if (nCols == 2)
            {
                result.m_waveLength.push_back(col1);
                result.m_crossSection.push_back(col2);
            }
        }

        fclose(f);

        return true;
    }

    bool SaveCrossSectionFile(const std::string& fullFilePath, const Evaluation::CCrossSectionData& data)
    {
        if (data.m_crossSection.size() == 0 && data.m_waveLength.size() == 0)
        {
            return false; // cannot save an empty reference.
        }

        FILE* f = fopen(fullFilePath.c_str(), "w");
        if (nullptr == f)
        {
            return false;
        }

        const bool twoColumns = data.m_waveLength.size() == data.m_crossSection.size();
        const size_t length = data.m_crossSection.size();

        for (size_t ii = 0; ii < length; ++ii)
        {
            if (twoColumns)
            {
                fprintf(f, "%.9lf %.9le\n", data.m_waveLength[ii], data.m_crossSection[ii]);
            }
            else
            {
                fprintf(f, "%.9le\n", data.m_crossSection[ii]);
            }
        }

        fclose(f);

        return true;
    }

    bool ReadSpectrum(const std::string& fullFilePath, CSpectrum& result)
    {
        if (fullFilePath.size() < 5)
        {
            return false; // invalid filename
        }
        const std::string fileEnding = Right(fullFilePath, 4);

        if (EqualsIgnoringCase(fileEnding, ".std"))
        {
            return SpectrumIO::CSTDFile::ReadSpectrum(result, fullFilePath);
        }
        else if (EqualsIgnoringCase(fileEnding, ".txt"))
        {
            return SpectrumIO::CTXTFile::ReadSpectrum(result, fullFilePath);
        }
        return false; // unknown file format.
    }
}