#include <SpectralEvaluation/File/STDFile.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/StringUtils.h>
#include <algorithm>
#include <cstring>
#include <cmath>

#include <iostream>

namespace novac
{
constexpr const char* elevationAngleStr = "ElevationAngle = ";
constexpr const char* azimuthAngleStr = "AzimuthAngle = ";
constexpr const char* temperatureStr = "Temperature = ";
constexpr const char* wavelengthStr = "Wavelength = ";

constexpr const char* MarkerStr = "Marker = ";
constexpr const char* MathLowStr = "MathLow = ";
constexpr const char* MathHighStr = "MathHigh = ";
constexpr const char* MinChannelStr = "MinChannel = ";
constexpr const char* MaxChannelStr = "MaxChannel = ";
constexpr const char* CalibPolynomialOrderStr = "CalibPolynomialOrder = ";
constexpr const char* CalibPolynomialStr = "CalibrationPolynomial = ";

static bool AttemptParseDouble(const std::vector<char>& stringToParse, const char* propertyName, double& value)
{
    const char* pt = strstr(stringToParse.data(), propertyName);
    if (nullptr != pt)
    {
        if (sscanf(pt + strlen(propertyName), "%lf", &value) == 1)
        {
            return true;
        }
    }

    return false;
}

static bool AttemptParseInt(const std::vector<char>& stringToParse, const char* propertyName, int& value)
{
    const char* pt = strstr(stringToParse.data(), propertyName);
    if (nullptr != pt)
    {
        if (sscanf(pt + strlen(propertyName), "%d", &value) == 1)
        {
            return true;
        }
    }

    return false;
}

static bool AttemptParseVector(const std::vector<char>& stringToParse, const char* propertyName, std::vector<double>& values)
{
    if (nullptr != strstr(stringToParse.data(), propertyName))
    {
        const char* pt = strstr(stringToParse.data(), propertyName) + strlen(propertyName);
        if (pt == nullptr)
        {
            return false;
        }

        if (nullptr != strstr(pt, "])"))
        {
            pt = strstr(pt, "])") + strlen("])");
        }
        // read array of space separated values
        const char* nextSeparator = pt;
        double temporaryDouble = 0.0;
        while (nextSeparator != nullptr)
        {
            if (sscanf(nextSeparator, "%lf", &temporaryDouble) == 1)
            {
                values.push_back(temporaryDouble);
            }
            nextSeparator = strstr(nextSeparator + 1, " ");
        }

        return values.size() > 0;
    }

    return false;
}

// Split the vector-of-characters into two substrings, at the position of the 'separator' character.
static std::pair<std::string, std::string> Split(std::vector<char>& str, const char* separator)
{
    char* pt = strstr(str.data(), separator);
    if (pt == nullptr)
    {
        return std::make_pair("", "");
    }
    *pt = '\0';
    std::string first = str.data();
    std::string second = (pt + strlen(separator));
    return std::make_pair(first, second);
}

bool CSTDFile::ReadSpectrum(CSpectrum& spec, const std::string& fileName)
{
    CSTDFile::ExtendedFormatInformation extendedInformation;
    return ReadSpectrum(spec, fileName, extendedInformation);
}

bool CSTDFile::ReadSpectrum(CSpectrum& spec, const std::string& fileName, CSTDFile::ExtendedFormatInformation& extendedInformation)
{
    try
    {
        extendedInformation.Clear();
        spec.m_wavelength.clear();

        FILE* f = fopen(fileName.c_str(), "r");
        if (f == nullptr)
        {
            throw std::invalid_argument("failed to open file.");
        }
        const int bufSize = 1024;
        char buffer[bufSize];
        int tmpInt, tmpInt2, tmpInt3;
        double tmpDbl;

        // 0. Clear the spectrum's information
        spec.Clear();

        // 1. the "GDBGMNUP" string that identifies a STD-file
        if (nullptr == fgets(buffer, bufSize, f))
        {
            fclose(f);
            throw std::invalid_argument("Failed to read mandatory header.");
        }
        if (0 != strncmp("GDBGMNUP", buffer, strlen("GDBGMNUP")))
        {
            fclose(f);
            throw std::invalid_argument("File does not contain mandatory header.");
        }

        // 2. The version number (always 1)
        if (nullptr == fgets(buffer, bufSize, f))
        {
            fclose(f);
            throw std::invalid_argument("File does not contain mandatory version.");
        }

        // 3. The spectrum length
        if (nullptr == fgets(buffer, bufSize, f))
        {
            fclose(f);
            throw std::invalid_argument("File does not contain mandatory length.");
        }
        if (1 > sscanf(buffer, "%d", &tmpInt))
        {
            fclose(f);
            throw std::invalid_argument("Failed to read first value.");
        }
        spec.m_length = std::min((long)tmpInt, (long)MAX_SPECTRUM_LENGTH);

        // 4. The spectrum data
        for (int i = 0; i < spec.m_length; ++i)
        {
            if (nullptr == fgets(buffer, bufSize, f))
            {
                fclose(f);
                throw std::invalid_argument("File does not contain enough values.");
            }
            if (1 > sscanf(buffer, "%lf", &tmpDbl))
            {
                fclose(f);
                throw std::invalid_argument("Failed to parse data.");
            }
            spec.m_data[i] = tmpDbl;
        }

        // 5. The fileName (ignore)
        if (nullptr == fgets(buffer, bufSize, f))
        {
            fclose(f);
            throw std::invalid_argument("Failed to parse filename.");
        }

        // 6. The detector (ignore)
        if (nullptr == fgets(buffer, bufSize, f))
        {
            fclose(f);
            throw std::invalid_argument("Failed to parse spectrometer model.");
        }
        spec.m_info.m_specModelName = std::string(buffer);
        Trim(spec.m_info.m_specModelName, " \t\r\n"); // remove initial and trailing spaces and newline characters

        // 7. The spectrometer
        if (nullptr == fgets(buffer, bufSize, f))
        {
            fclose(f);
            throw std::invalid_argument("Failed to parse spectrometer name.");
        }
        spec.m_info.m_device = std::string(buffer);
        Trim(spec.m_info.m_device, " \t\r\n"); // remove initial and trailing spaces and newline characters

        // 8. The date
        if (nullptr == fgets(buffer, bufSize, f))
        {
            fclose(f);
            throw std::invalid_argument("Failed to read date.");
        }
        if (strstr(buffer, "/"))
        {
            // DOASIS sometimes writes the date as MM/DD/YYYY
            if (3 > sscanf(buffer, "%d/%d/%d", &tmpInt2, &tmpInt3, &tmpInt))
            {
                fclose(f);
                throw std::invalid_argument("Failed to parse date.");
            }
        }
        else
        {
            if (3 > sscanf(buffer, "%d.%d.%d", &tmpInt3, &tmpInt2, &tmpInt))
            {
                fclose(f);
                throw std::invalid_argument("Failed to parse date.");
            }
        }
        if (tmpInt > 31)
        {
            spec.m_info.m_startTime.year = (unsigned short)((tmpInt < 1900) ? tmpInt + 2000 : tmpInt);
            spec.m_info.m_startTime.month = (unsigned char)tmpInt2;
            spec.m_info.m_startTime.day = (unsigned char)tmpInt3;
        }
        else
        {
            spec.m_info.m_startTime.year = (unsigned short)((tmpInt3 < 1900) ? tmpInt3 + 2000 : tmpInt3);
            spec.m_info.m_startTime.month = (unsigned char)tmpInt2;
            spec.m_info.m_startTime.day = (unsigned char)tmpInt;
        }
        spec.m_info.m_stopTime.year = spec.m_info.m_startTime.year;
        spec.m_info.m_stopTime.month = spec.m_info.m_startTime.month;
        spec.m_info.m_stopTime.day = spec.m_info.m_startTime.day;

        // 9. The starttime
        if (nullptr == fgets(buffer, bufSize, f))
        {
            fclose(f);
            throw std::invalid_argument("Failed to read starttime.");
        }
        if (3 > sscanf(buffer, "%d:%d:%d", &tmpInt, &tmpInt2, &tmpInt3))
        {
        }
        else
        {
            spec.m_info.m_startTime.hour = (unsigned char)tmpInt;
            spec.m_info.m_startTime.minute = (unsigned char)tmpInt2;
            spec.m_info.m_startTime.second = (unsigned char)tmpInt3;
        }

        // 10. The stoptime
        if (nullptr == fgets(buffer, bufSize, f))
        {
            fclose(f);
            throw std::invalid_argument("Failed to read stop time.");
        }
        if (3 > sscanf(buffer, "%d:%d:%d", &tmpInt, &tmpInt2, &tmpInt3))
        {
            fclose(f);
            throw std::invalid_argument("Failed to parse stop time.");
        }
        spec.m_info.m_stopTime.hour = (unsigned char)tmpInt;
        spec.m_info.m_stopTime.minute = (unsigned char)tmpInt2;
        spec.m_info.m_stopTime.second = (unsigned char)tmpInt3;

        // 11. The start wavelength (ignore)
        if (nullptr == fgets(buffer, bufSize, f))
        {
            fclose(f);
            throw std::invalid_argument("Failed to read first wavelength.");
        }

        // 12. The stop wavelength (ignore)
        if (nullptr == fgets(buffer, bufSize, f))
        {
            fclose(f);
            throw std::invalid_argument("Failed to read last wavelength.");
        }

        // 13. The number of scans
        if (nullptr == fgets(buffer, bufSize, f))
        {
            fclose(f);
            throw std::invalid_argument("Failed to read number of scans.");
        }
        if (1 > sscanf(buffer, "SCANS %d", &tmpInt))
        {
            fclose(f);
            throw std::invalid_argument("Failed to parse number of scans.");
        }
        spec.m_info.m_numSpec = tmpInt;

        // 14. The integration time
        if (nullptr == fgets(buffer, bufSize, f))
        {
            fclose(f);
            throw std::invalid_argument("Failed to read integration time.");
        }
        if (1 > sscanf(buffer, "INT_TIME %lf", &tmpDbl))
        {
            fclose(f);
            throw std::invalid_argument("Failed to parse integration time.");
        }
        spec.m_info.m_exposureTime = (int)tmpDbl;

        // 15. The site
        if (nullptr == fgets(buffer, bufSize, f))
        {
            fclose(f);
            throw std::invalid_argument("Failed to read site.");
        }
        buffer[strlen(buffer) - 1] = 0; // <-- Remove the trailing newline character
        spec.m_info.m_name = std::string(buffer + 5);

        // 15. The longitude
        if (nullptr == fgets(buffer, bufSize, f))
        {
            fclose(f);
            throw std::invalid_argument("Failed to read longitude.");
        }
        if (1 > sscanf(buffer, "LONGITUDE %lf", &tmpDbl))
        {
            fclose(f);
            throw std::invalid_argument("Failed to parse longitude.");
        }
        spec.m_info.m_gps.m_longitude = tmpDbl;

        // 15. The latitude
        if (nullptr == fgets(buffer, bufSize, f))
        {
            fclose(f);
            throw std::invalid_argument("Failed to read latitude.");
        }
        if (1 > sscanf(buffer, "LATITUDE %lf", &tmpDbl))
        {
            fclose(f);
            throw std::invalid_argument("Failed to parse latitude.");
        }
        spec.m_info.m_gps.m_latitude = tmpDbl;

        // ----------- EXTENDED STD ------------------
        // - if the file is in the extended STD-format then we can continue here... -

        std::vector<char> szLine(65535);
        while (fgets(szLine.data(), (int)szLine.size(), f))
        {
            // Read in scanAngle
            if (AttemptParseDouble(szLine, elevationAngleStr, tmpDbl))
            {
                if (std::abs(tmpDbl) < 360.0)
                {
                    spec.m_info.m_scanAngle = (float)tmpDbl;
                }
                continue;
            }

            // Read in scanAngle2
            if (AttemptParseDouble(szLine, azimuthAngleStr, tmpDbl))
            {
                if (std::abs(tmpDbl) < 360.0)
                {
                    spec.m_info.m_scanAngle2 = (float)tmpDbl;
                }
                continue;
            }

            if (AttemptParseDouble(szLine, temperatureStr, tmpDbl))
            {
                if (std::abs(tmpDbl) < 100.0)
                {
                    spec.m_info.m_temperature = (float)tmpDbl;
                }
                continue;
            }

            if (AttemptParseDouble(szLine, MarkerStr, tmpDbl))
            {
                extendedInformation.Marker = tmpDbl;
                continue;
            }
            else if (AttemptParseInt(szLine, MathLowStr, tmpInt))
            {
                extendedInformation.MathLow = tmpInt;
                continue;
            }
            else if (AttemptParseInt(szLine, MathHighStr, tmpInt))
            {
                extendedInformation.MathHigh = tmpInt;
                continue;
            }
            else if (AttemptParseInt(szLine, MinChannelStr, tmpInt))
            {
                extendedInformation.MinChannel = tmpInt;
                continue;
            }
            else if (AttemptParseInt(szLine, MaxChannelStr, tmpInt))
            {
                extendedInformation.MaxChannel = tmpInt;
                continue;
            }

            // Calibration polynomial
            if (AttemptParseVector(szLine, CalibPolynomialStr, extendedInformation.calibrationPolynomial))
            {
                continue;
            }

            // Wavelength calibration
            if (AttemptParseVector(szLine, wavelengthStr, spec.m_wavelength))
            {
                continue;
            }

            // Unknown property
            if (nullptr != strstr(szLine.data(), " = "))
            {
                auto properties = Split(szLine, " = ");
                extendedInformation.additionalProperties.push_back(properties);
            }
        }

        // Get the intensity
        double numSpec_inv = 1 / (double)spec.NumSpectra();
        spec.m_info.m_peakIntensity = static_cast<float>(spec.MaxValue() * numSpec_inv);
        spec.m_info.m_offset = static_cast<float>(spec.GetOffset() * numSpec_inv);

        // If the intensity is saved in the ddmm.mmmm - format (which is what
        //	the gps-reciever sends out), convert it to the dd.dddddd format
        if (std::abs(spec.m_info.m_gps.m_latitude) > 180)
        {
            spec.m_info.m_gps.m_latitude = CGPSData::DoubleToAngle(spec.m_info.m_gps.m_latitude);
        }
        if (std::abs(spec.m_info.m_gps.m_longitude) > 180)
        {
            spec.m_info.m_gps.m_longitude = CGPSData::DoubleToAngle(spec.m_info.m_gps.m_longitude);
        }

        fclose(f);
        return true;
    }
    catch (std::exception& e)
    {
        std::cout << "Failed to read .std file '" << fileName << "'. Error: " << e.what() << std::endl;
        return false;
    }
}

/// <summary>
/// Writes the mandatory part of the STD file to the provided (opened) file handle.
/// </summary>
/// <param name="spec">The spectrum data to save</param>
/// <param name="fileName">The name of the file we're saving</param>
/// <param name="f">An already opened file handle to which the data will be saved.</param>
static void WriteBasicData(const CSpectrum& spec, const std::string& fileName, FILE* f)
{
    const CSpectrumInfo& info = spec.m_info;

    fprintf(f, "GDBGMNUP\n1\n");
    fprintf(f, "%ld\n", spec.m_length);
    for (long i = 0; i < spec.m_length; ++i)
    {
        if (std::abs(spec.m_data[i] - floor(spec.m_data[i])) > 1e-9)
        {
            fprintf(f, "%.9lf\n", spec.m_data[i]);
        }
        else if (std::abs(spec.m_data[i]) < 1e-5)
        {
            fprintf(f, "%.9e\n", spec.m_data[i]);
        }
        else
        {
            fprintf(f, "%.0lf\n", spec.m_data[i]);
        }
    }

    int index = ReverseFind(fileName, '\\');
    if (index != 0)
    {
        std::string leftSubString = Left(fileName, index);
        fprintf(f, "%s\n", leftSubString.c_str());
    }
    else
    {
        fprintf(f, "%s\n", fileName.c_str());
    }
    fprintf(f, "%s\n", info.m_specModelName.c_str()); // the spectrometer model
    fprintf(f, "%s\n", info.m_device.c_str()); // the spectrometer serial number

    if (info.m_startTime.year == info.m_startTime.month && info.m_startTime.month == info.m_startTime.day && info.m_startTime.day == 0)
    {
        fprintf(f, "01.01.01\n"); /* Default date */
    }
    else
    {
        fprintf(f, "%02d.%02d.%02d\n", info.m_startTime.year, info.m_startTime.month, info.m_startTime.day); // the date
    }
    fprintf(f, "%02d:%02d:%02d\n", info.m_startTime.hour, info.m_startTime.minute, info.m_startTime.second);
    fprintf(f, "%02d:%02d:%02d\n", info.m_stopTime.hour, info.m_stopTime.minute, info.m_stopTime.second);

    fprintf(f, "0.0\n0.0\n"); // the start and stop wavelengths

    fprintf(f, "SCANS %ld\n", info.m_numSpec);
    fprintf(f, "INT_TIME %ld\n", info.m_exposureTime);

    fprintf(f, "SITE %s\n", info.m_name.c_str());
    fprintf(f, "LONGITUDE %.6lf\n", info.m_gps.m_longitude);
    fprintf(f, "LATITUDE %.6lf\n", info.m_gps.m_latitude);
}

/// <summary>
/// Writes the extended (DOASIS specific) part of the STD file to the provided (opened) file handle
/// </summary>
/// <param name="spec">The spectrum data to save</param>
/// <param name="fileName">The name of the file we're saving</param>
/// <param name="f">An already opened file handle to which the data will be saved.</param>
void WriteExtendedData(const CSpectrum& spec, const std::string& fileName, const CSTDFile::ExtendedFormatInformation& extendedInformation, FILE* f)
{
    const CSpectrumInfo& info = spec.m_info;

    if (spec.m_wavelength.size() == static_cast<size_t>(spec.m_length))
    {
        // DOASIS Specific format of wavelength data
        fprintf(f, wavelengthStr);
        fprintf(f, "(System.Double[%ld])", spec.m_length);
        for (size_t ii = 0; ii < static_cast<size_t>(spec.m_length - 1); ++ii)
        {
            fprintf(f, "%lf ", spec.m_wavelength[ii]);
        }
        fprintf(f, "%lf\n", spec.m_wavelength.back());
    }

    if (extendedInformation.calibrationPolynomial.size() > 0)
    {
        fprintf(f, CalibPolynomialOrderStr);
        fprintf(f, "%lu\n", static_cast<unsigned long>(extendedInformation.calibrationPolynomial.size() - 1));

        fprintf(f, CalibPolynomialStr);
        fprintf(f, "(System.Double[%lu])", static_cast<unsigned long>(extendedInformation.calibrationPolynomial.size()));
        for (size_t ii = 0; ii < extendedInformation.calibrationPolynomial.size() - 1; ++ii)
        {
            fprintf(f, "%.9g ", extendedInformation.calibrationPolynomial[ii]);
        }
        fprintf(f, "%.9g\n", extendedInformation.calibrationPolynomial.back());
    }

    fprintf(f, "Author = \"\"\n");
    fprintf(f, "Average = %.2lf\n", spec.AverageValue());
    fprintf(f, "AzimuthAngle = 0\n");
    fprintf(f, "Delta = 0\n");
    fprintf(f, "DeltaRel = 0\n");
    fprintf(f, "Deviation = 0\n");
    fprintf(f, "Device = \"\"\n");
    fprintf(f, "ElevationAngle = %.2lf\n", info.m_scanAngle);
    fprintf(f, "ExposureTime = %ld\n", info.m_exposureTime);
    fprintf(f, "FileName = %s\n", fileName.c_str());
    fprintf(f, "FitHigh = 0\n");
    fprintf(f, "FitLow = 0\n");
    fprintf(f, "Gain = 0\n");
    fprintf(f, "IntegrationMethod = Sum\n");
    fprintf(f, "Latitude = %.6lf\n", info.m_gps.m_latitude);
    fprintf(f, "LightPath = 0\n");
    fprintf(f, "LightSource = \"\"\n");
    fprintf(f, "Longitude = %.6lf\n", info.m_gps.m_longitude);
    fprintf(f, "%s %lf\n", MarkerStr, extendedInformation.Marker);
    fprintf(f, "%s %ld\n", MathHighStr, extendedInformation.MathHigh >= 0 ? extendedInformation.MathHigh : spec.m_length);
    fprintf(f, "%s %d\n", MathLowStr, extendedInformation.MathLow >= 0 ? extendedInformation.MathLow : 0);
    fprintf(f, "Max = %.3lf\n", spec.MaxValue());
    fprintf(f, "%s %ld\n", MaxChannelStr, extendedInformation.MaxChannel >= 0 ? extendedInformation.MaxChannel : spec.m_length);
    fprintf(f, "Min = %.3lf\n", spec.MinValue());
    fprintf(f, "%s %d\n", MinChannelStr, extendedInformation.MinChannel >= 0 ? extendedInformation.MinChannel : 0);
    fprintf(f, "MultiChannelCounter = 0\n");
    fprintf(f, "Name = \"%s\"\n", info.m_name.c_str());
    fprintf(f, "NumScans = %ld\n", info.m_numSpec);
    fprintf(f, "OpticalDensity = 0\n");
    fprintf(f, "OpticalDensityCenter = %ld\n", spec.m_length / 2);
    fprintf(f, "OpticalDensityLeft = 0\n");
    fprintf(f, "OpticalDensityRight = %ld\n", spec.m_length - 1);
    fprintf(f, "Pressure = 0\n");
    fprintf(f, "Remark = \"\"\n");
    fprintf(f, "ScanGeometry = 0\n"); //(DoasCore.Math.ScanGeometry)SAZ: 137.41237083135 SZA: 31.5085943481828 LAZ: 298.523110145623 LAZ: 129.285101310559 Date: 1/5/2007 10:35:07 Lat.: 0 Lon.: 0\n");
    fprintf(f, "ScanMax = 0\n");
    fprintf(f, "Temperature = %.2f\n", info.m_temperature);
    fprintf(f, "Variance = 0\n");

    for (const auto& prop : extendedInformation.additionalProperties)
    {
        fprintf(f, "%s = %s\n", prop.first.c_str(), prop.second.c_str());
    }
}

bool CSTDFile::WriteSpectrum(const CSpectrum* spec, const std::string& fileName, int extendedFormat)
{
    return WriteSpectrum(*spec, fileName, extendedFormat);
}

bool CSTDFile::WriteSpectrum(const CSpectrum& spec, const std::string& fileName, int extendedFormat)
{
    FILE* f = fopen(fileName.c_str(), "w");
    if (f == nullptr)
    {
        return false;
    }

    WriteBasicData(spec, fileName, f);

    if (extendedFormat)
    {
        CSTDFile::ExtendedFormatInformation defaultExt;
        defaultExt.Marker = spec.m_length / 2.0;
        defaultExt.MinChannel = 0;
        defaultExt.MaxChannel = spec.m_length;
        defaultExt.MathLow = 0;
        defaultExt.MathHigh = spec.m_length;
        WriteExtendedData(spec, fileName, defaultExt, f);
    }

    fclose(f);

    return true;
}

bool CSTDFile::WriteSpectrum(const CSpectrum& spec, const std::string& fileName, const CSTDFile::ExtendedFormatInformation& extendedInformation)
{
    FILE* f = fopen(fileName.c_str(), "w");
    if (f == nullptr)
    {
        return false;
    }

    WriteBasicData(spec, fileName, f);

    WriteExtendedData(spec, fileName, extendedInformation, f);

    fclose(f);

    return true;
}

}

