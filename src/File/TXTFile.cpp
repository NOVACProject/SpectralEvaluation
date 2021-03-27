#include <SpectralEvaluation/File/TXTFile.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>

namespace novac
{
// TODO: move these functions to a common helper
int ReadFromFile(FILE* f, const char* formatStr, double& col1, double& col2);
int GetFileFormat(const char* string, char* format);

bool CTXTFile::ReadSpectrum(CSpectrum& spec, const std::string& fileName)
{
    // Open the file
    FILE* f = fopen(fileName.c_str(), "r");
    if (f == nullptr)
    {
        return false;
    }

    // Clear all the information in the spectrum
    spec.Clear();

    // allocate some space for the wavlength data (if there is any)
    std::vector<double> wavelengths(MAX_SPECTRUM_LENGTH);
    bool containsWavelengthData = false;

    // Get the file-format and number of columns from the file.
    char tempBuffer[8192];
    char format[256];
    if (nullptr == fgets(tempBuffer, 8191, f))
    {
        fclose(f);
        return false;
    }
    const int numColumns = novac::GetFileFormat(tempBuffer, format);
    fseek(f, 0, SEEK_SET);

    // Simply read the spectrum, one pixel at a time
    int length = 0;
    while (length < MAX_SPECTRUM_LENGTH)
    {
        double col1 = 0.0;
        double col2 = 0.0;
        int nCols = novac::ReadFromFile(f, format, col1, col2);

        if (nCols != numColumns)
        {
            break;
        }

        if (nCols == 1)
        {
            spec.m_data[length] = col1;
        }
        else if (nCols == 2)
        {
            containsWavelengthData = true;
            wavelengths[length] = col1;
            spec.m_data[length] = col2;
        }
        else
        {
            break;
        }

        ++length;
    }
    spec.m_length = length;

    if (containsWavelengthData)
    {
        spec.m_wavelength = std::vector<double>(begin(wavelengths), begin(wavelengths) + length);
    }
    else
    {
        spec.m_wavelength.clear();
    }

    // close the file before we return
    fclose(f);
    return true;
}

bool CTXTFile::WriteSpectrum(const CSpectrum* spec, const std::string& fileName)
{
    if (spec == nullptr)
    {
        return false;
    }
    return WriteSpectrum(*spec, fileName);
}

bool CTXTFile::WriteSpectrum(const CSpectrum& spec, const std::string& fileName)
{
    // Open the file
    FILE* f = fopen(fileName.c_str(), "w");
    if (nullptr == f)
    {
        return false;
    }

    // Write the spectrum, one data-point at a time
    if (spec.m_wavelength.size() == (size_t)spec.m_length)
    {
        // Include the wavelength data as a separate column
        for (int k = 0; k < spec.m_length; ++k)
        {
            if (fprintf(f, "%lf\t%lf\n", spec.m_wavelength[k], spec.m_data[k]) < 0)
            {
                // an error has occured, try to close the file
                fclose(f);
                return false;
            }
        }
    }
    else
    {
        for (int k = 0; k < spec.m_length; ++k)
        {
            if (fprintf(f, "%lf\n", spec.m_data[k]) < 0)
            {
                // an error has occured, try to close the file
                fclose(f);
                return false;
            }
        }
    }

    fclose(f); // <-- close the file before we return
    return true;
}

}
