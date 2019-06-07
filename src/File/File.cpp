#include <SpectralEvaluation/File/File.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <cstring>

namespace FileIo
{
    // @return the number of columsn read.
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
        if(nullptr == string || nullptr == format) 
            return 0;
        if(0 == strlen(string))
            return 0;

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

    bool ReadCrossSectionFile(const std::string& fullFilePath, Evaluation::CCrossSectionData& result, bool saveAsWavelength)
    {
        FILE* f = fopen(fullFilePath.c_str(), "r");
        if (nullptr == f)
        {
            return false;
        }

        // Get the file-format and number of columns from the file.
        char tempBuffer[8192];
        char format[256];
        if(nullptr == fgets(tempBuffer, 8191, f))
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

        // close the file before we return
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
        const size_t length   = data.m_crossSection.size();

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

}