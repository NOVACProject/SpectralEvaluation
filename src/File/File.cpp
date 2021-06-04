#include <SpectralEvaluation/StringUtils.h>
#include <SpectralEvaluation/File/File.h>
#include <SpectralEvaluation/File/STDFile.h>
#include <SpectralEvaluation/File/TXTFile.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <cstring>
#include <sstream>

namespace novac
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

bool ReadCrossSectionFromStdFile(const std::string& fullFilePath, CCrossSectionData& result, bool saveAsWavelength)
{
    CSpectrum spectrum;
    if (!CSTDFile::ReadSpectrum(spectrum, fullFilePath))
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

bool ReadCrossSectionFile(const std::string& fullFilePath, CCrossSectionData& result, bool saveAsWavelength)
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

bool SaveCrossSectionFile(const std::string& fullFilePath, const CCrossSectionData& data)
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

bool SaveCrossSectionFile(const std::string& fullFilePath, const CSpectrum& spectrum)
{
    if (spectrum.m_length == 0)
    {
        return false; // cannot save an empty reference.
    }

    FILE* f = fopen(fullFilePath.c_str(), "w");
    if (nullptr == f)
    {
        return false;
    }

    const bool twoColumns = spectrum.m_wavelength.size() == spectrum.m_length;
    const size_t length = spectrum.m_length;

    for (size_t ii = 0; ii < length; ++ii)
    {
        if (twoColumns)
        {
            fprintf(f, "%.9lf %.9le\n", spectrum.m_wavelength[ii], spectrum.m_data[ii]);
        }
        else
        {
            fprintf(f, "%.9le\n", spectrum.m_data[ii]);
        }
    }

    fclose(f);

    return true;
}

bool SaveDataToFile(const std::string& fullFilePath, const std::vector<double>& data)
{
    if (data.size() == 0)
    {
        return false; // cannot save an empty reference.
    }

    FILE* f = fopen(fullFilePath.c_str(), "w");
    if (nullptr == f)
    {
        return false;
    }

    for (size_t ii = 0; ii < data.size(); ++ii)
    {

        fprintf(f, "%.9le\n", data[ii]);
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
        return CSTDFile::ReadSpectrum(result, fullFilePath);
    }
    else if (EqualsIgnoringCase(fileEnding, ".txt"))
    {
        return CTXTFile::ReadSpectrum(result, fullFilePath);
    }
    return false; // unknown file format.
}

std::string EnsureFilenameHasSuffix(const std::string& fullFilePath, const std::string& defaultSuffix)
{
    if (fullFilePath.size() == 0)
    {
        return fullFilePath;
    }

    const size_t lastPeriod = fullFilePath.rfind('.');
    if (lastPeriod == fullFilePath.npos)
    {
        // no period found
        std::stringstream correctedFileName;
        correctedFileName << fullFilePath << "." << defaultSuffix;
        return correctedFileName.str();
    }

    const size_t lastForwardSlash = fullFilePath.rfind('\\');
    const size_t lastBackwardSlash = fullFilePath.rfind('/');

    if ((lastForwardSlash != fullFilePath.npos && lastPeriod > lastForwardSlash) ||
        (lastBackwardSlash != fullFilePath.npos && lastPeriod > lastBackwardSlash))
    {
        return fullFilePath;
    }
    else
    {
        // no period found _after_ the last path-separator character. Add the suffix.
        std::stringstream correctedFileName;
        correctedFileName << fullFilePath << "." << defaultSuffix;
        return correctedFileName.str();
    }
}

bool SaveInstrumentCalibration(const std::string& fullFilePath, const CSpectrum& instrumentLineShape, const std::vector<double>& pixelToWavelengthMapping)
{
    if (instrumentLineShape.m_length == 0 || instrumentLineShape.m_wavelength.size() == 0)
    {
        return false; // cannot save an empty reference.
    }

    FILE* f = fopen(fullFilePath.c_str(), "w");
    if (nullptr == f)
    {
        return false;
    }
    fprintf(f, "{\n");

    {
        fprintf(f, "\t\"InstrumentLineShape\": {\n");

        fprintf(f, "\t\t\"Wavelength\" : [");
        for (size_t ii = 0; ii < instrumentLineShape.m_length - 1; ++ii)
        {
            fprintf(f, "%.9lf, ", instrumentLineShape.m_wavelength[ii]);
        }
        fprintf(f, "%.9lf", instrumentLineShape.m_wavelength[instrumentLineShape.m_length - 1]);
        fprintf(f, "],\n");

        fprintf(f, "\t\t\"Intensity\" : [");
        for (size_t ii = 0; ii < instrumentLineShape.m_length - 1; ++ii)
        {
            fprintf(f, "%.9lf, ", instrumentLineShape.m_data[ii]);
        }
        fprintf(f, "%.9lf]},\n", instrumentLineShape.m_data[instrumentLineShape.m_length - 1]);
    }

    {
        fprintf(f, "\t\"PixelToWavelengthMapping\": [");
        for (size_t ii = 0; ii < pixelToWavelengthMapping.size() - 1; ++ii)
        {
            fprintf(f, "%.9lf, ", pixelToWavelengthMapping[ii]);
        }
        fprintf(f, "%.9lf]\n", pixelToWavelengthMapping.back());
    }

    fprintf(f, "}\n");
    fclose(f);
    return true;
}

bool ReadInstrumentCalibration(const std::string& fullFilePath, CSpectrum& instrumentLineShape, std::vector<double>& pixelToWavelengthMapping)
{

    return false;
}


}