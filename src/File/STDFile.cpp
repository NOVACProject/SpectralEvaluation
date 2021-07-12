#include <SpectralEvaluation/File/STDFile.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/StringUtils.h>
#include <algorithm>
#include <cstring>
#include <cmath>

#include <iostream>

namespace novac
{
constexpr char* elevationAngleStr = "ElevationAngle = ";
constexpr char* azimuthAngleStr = "AzimuthAngle = ";
constexpr char* temperatureStr = "Temperature = ";
constexpr char* wavelengthStr = "Wavelength = ";

bool CSTDFile::ReadSpectrum(CSpectrum& spec, const std::string& fileName)
{
    FILE* f = fopen(fileName.c_str(), "r");
    if (f == nullptr)
    {
        return false;
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
        fclose(f); return false;
    }
    if (0 != strncmp("GDBGMNUP\n", buffer, strlen(buffer)))
    {
        fclose(f); return false;
    }

    // 2. The version number (always 1)
    if (nullptr == fgets(buffer, bufSize, f))
    {
        fclose(f); return false;
    }

    // 3. The spectrum length
    if (nullptr == fgets(buffer, bufSize, f))
    {
        fclose(f); return false;
    }
    if (1 > sscanf(buffer, "%d", &tmpInt))
    {
        fclose(f); return false;
    }
    spec.m_length = std::min((long)tmpInt, (long)MAX_SPECTRUM_LENGTH);

    // 4. The spectrum data
    for (int i = 0; i < spec.m_length; ++i)
    {
        if (nullptr == fgets(buffer, bufSize, f))
        {
            fclose(f); return false;
        }
        if (1 > sscanf(buffer, "%lf", &tmpDbl))
        {
            fclose(f); return false;
        }
        spec.m_data[i] = tmpDbl;
    }

    // 5. The fileName (ignore)
    if (nullptr == fgets(buffer, bufSize, f))
    {
        fclose(f); return false;
    }

    // 6. The detector (ignore)
    if (nullptr == fgets(buffer, bufSize, f))
    {
        fclose(f); return false;
    }

    // 7. The spectrometer
    if (nullptr == fgets(buffer, bufSize, f))
    {
        fclose(f); return false;
    }
    buffer[strlen(buffer) - 1] = 0; // <-- Remove the trailing newline character
    spec.m_info.m_device = std::string(buffer);

    // 8. The date
    if (nullptr == fgets(buffer, bufSize, f))
    {
        fclose(f); return false;
    }
    if (strstr(buffer, "/"))
    {
        // DOASIS sometimes writes the date as MM/DD/YYYY
        if (3 > sscanf(buffer, "%d/%d/%d", &tmpInt2, &tmpInt3, &tmpInt))
        {
            fclose(f); return false;
        }
    }
    else {
        if (3 > sscanf(buffer, "%d.%d.%d", &tmpInt3, &tmpInt2, &tmpInt))
        {
            fclose(f); return false;
        }
    }
    if (tmpInt > 31)
    {
        spec.m_info.m_startTime.year = (unsigned short)((tmpInt < 1900) ? tmpInt + 2000 : tmpInt);
        spec.m_info.m_startTime.month = (unsigned char)tmpInt2;
        spec.m_info.m_startTime.day = (unsigned char)tmpInt3;
    }
    else {
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
        fclose(f); return false;
    }
    if (3 > sscanf(buffer, "%d:%d:%d", &tmpInt, &tmpInt2, &tmpInt3))
    {
    }
    else {
        spec.m_info.m_startTime.hour = (unsigned char)tmpInt;
        spec.m_info.m_startTime.minute = (unsigned char)tmpInt2;
        spec.m_info.m_startTime.second = (unsigned char)tmpInt3;
    }

    // 10. The stoptime
    if (nullptr == fgets(buffer, bufSize, f))
    {
        fclose(f); return false;
    }
    if (3 > sscanf(buffer, "%d:%d:%d", &tmpInt, &tmpInt2, &tmpInt3))
    {
        fclose(f); return false;
    }
    spec.m_info.m_stopTime.hour = (unsigned char)tmpInt;
    spec.m_info.m_stopTime.minute = (unsigned char)tmpInt2;
    spec.m_info.m_stopTime.second = (unsigned char)tmpInt3;

    // 11. The start wavelength (ignore)
    if (nullptr == fgets(buffer, bufSize, f))
    {
        fclose(f); return false;
    }

    // 12. The stop wavelength (ignore)
    if (nullptr == fgets(buffer, bufSize, f))
    {
        fclose(f); return false;
    }

    // 13. The number of scans
    if (nullptr == fgets(buffer, bufSize, f))
    {
        fclose(f); return false;
    }
    if (1 > sscanf(buffer, "SCANS %d", &tmpInt))
    {
        fclose(f); return false;
    }
    spec.m_info.m_numSpec = tmpInt;

    // 14. The integration time
    if (nullptr == fgets(buffer, bufSize, f))
    {
        fclose(f); return false;
    }
    if (1 > sscanf(buffer, "INT_TIME %lf", &tmpDbl))
    {
        fclose(f); return false;
    }
    spec.m_info.m_exposureTime = (int)tmpDbl;

    // 15. The site
    if (nullptr == fgets(buffer, bufSize, f))
    {
        fclose(f); return false;
    }
    buffer[strlen(buffer) - 1] = 0; // <-- Remove the trailing newline character
    spec.m_info.m_name = std::string(buffer + 5);

    // 15. The longitude
    if (nullptr == fgets(buffer, bufSize, f))
    {
        fclose(f); return false;
    }
    if (1 > sscanf(buffer, "LONGITUDE %lf", &tmpDbl))
    {
        fclose(f); return false;
    }
    spec.m_info.m_gps.m_longitude = tmpDbl;

    // 15. The latitude
    if (nullptr == fgets(buffer, bufSize, f))
    {
        fclose(f); return false;
    }
    if (1 > sscanf(buffer, "LATITUDE %lf", &tmpDbl))
    {
        fclose(f); return false;
    }
    spec.m_info.m_gps.m_latitude = tmpDbl;

    // ----------- EXTENDED STD ------------------
    // - if the file is in the extended STD-format then we can continue here... -

    std::vector<char> szLine(65535);
    while (fgets(szLine.data(), (int)szLine.size(), f))
    {
        // Read in scanAngle
        if (nullptr != strstr(szLine.data(), elevationAngleStr))
        {
            char* pt = strstr(szLine.data(), elevationAngleStr);
            if (sscanf(pt + strlen(elevationAngleStr), "%lf", &tmpDbl) == 1)
            {
                if (fabs(tmpDbl) < 360.0)
                {
                    spec.m_info.m_scanAngle = (float)tmpDbl;
                }
            }
        }

        // Read in scanAngle2
        if (nullptr != strstr(szLine.data(), azimuthAngleStr))
        {
            char* pt = strstr(szLine.data(), azimuthAngleStr);
            if (sscanf(pt + strlen(azimuthAngleStr), "%lf", &tmpDbl) == 1)
            {
                if (fabs(tmpDbl) < 360.0)
                {
                    spec.m_info.m_scanAngle2 = (float)tmpDbl;
                }
            }
        }

        // Read in the temperature
        if (nullptr != strstr(szLine.data(), temperatureStr))
        {
            char* pt = strstr(szLine.data(), temperatureStr);
            if (sscanf(pt + strlen(temperatureStr), "%lf", &tmpDbl) == 1)
            {
                if (fabs(tmpDbl) < 100.0)
                {
                    spec.m_info.m_temperature = (float)tmpDbl;
                }
            }
        }

        // Wavelength calibration
        if (nullptr != strstr(szLine.data(), wavelengthStr))
        {
            char* pt = strstr(szLine.data(), wavelengthStr) + strlen(wavelengthStr);
            if (nullptr != strstr(pt, "])"))
            {
                pt = strstr(pt, "])") + strlen("])");
            }
            // read array of space separated values
            char* nextSeparator = pt;
            while (nextSeparator != nullptr)
            {
                if (sscanf(nextSeparator, "%lf", &tmpDbl) == 1)
                {
                    spec.m_wavelength.push_back(tmpDbl);
                }
                nextSeparator = strstr(nextSeparator + 1, " ");
            }
        }
    }

    // Get the intensity
    float numSpec_inv = 1 / (float)spec.NumSpectra();
    spec.m_info.m_peakIntensity = (float)spec.MaxValue() * numSpec_inv;
    spec.m_info.m_offset = (float)spec.GetOffset() * numSpec_inv;

    // If the intensity is saved in the ddmm.mmmm - format (which is what
    //	the gps-reciever sends out), convert it to the dd.dddddd format
    if (fabs(spec.m_info.m_gps.m_latitude) > 180)
        spec.m_info.m_gps.m_latitude = CGPSData::DoubleToAngle(spec.m_info.m_gps.m_latitude);
    if (fabs(spec.m_info.m_gps.m_longitude) > 180)
        spec.m_info.m_gps.m_longitude = CGPSData::DoubleToAngle(spec.m_info.m_gps.m_longitude);

    fclose(f);
    return true;
}

/// <summary>
/// Writes the mandatory part of the STD file to the provided (opened) file handle.
/// </summary>
/// <param name="spec">The spectrum data to save</param>
/// <param name="fileName">The name of the file we're saving</param>
/// <param name="f">An already opened file handle to which the data will be saved.</param>
void WriteBasicData(const CSpectrum& spec, const std::string& fileName, FILE* f)
{
    const CSpectrumInfo& info = spec.m_info;

    fprintf(f, "GDBGMNUP\n1\n");
    fprintf(f, "%ld\n", spec.m_length);
    for (long i = 0; i < spec.m_length; ++i)
    {
        if (fabs(spec.m_data[i] - floor(spec.m_data[i])) > 1e-9)
        {
            fprintf(f, "%.9lf\n", spec.m_data[i]);
        }
        else if (std::abs(spec.m_data[i] < 1e-5))
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
    fprintf(f, ".......\n"); // the detector
    fprintf(f, "%s\n", info.m_device.c_str()); // the spectrometer

    if (info.m_startTime.year == info.m_startTime.month && info.m_startTime.month == info.m_startTime.day && info.m_startTime.day == 0)
    {
        fprintf(f, "01.01.01\n"); /* Default date */
    }
    else {
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
        fprintf(f, "CalibPolynomialOrder = %lld\n", extendedInformation.calibrationPolynomial.size() - 1);

        fprintf(f, "CalibPolynomial = ");
        fprintf(f, "(System.Double[%lld])", extendedInformation.calibrationPolynomial.size());
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
    fprintf(f, "Marker = %lf\n", extendedInformation.Marker);
    fprintf(f, "MathHigh = %ld\n", extendedInformation.MathHigh >= 0 ? extendedInformation.MathHigh : spec.m_length);
    fprintf(f, "MathLow = %ld\n", extendedInformation.MathLow >= 0 ? extendedInformation.MathLow : 0);
    fprintf(f, "Max = %.3lf\n", spec.MaxValue());
    fprintf(f, "MaxChannel = %ld\n", extendedInformation.MaxChannel >= 0 ? extendedInformation.MaxChannel : spec.m_length);
    fprintf(f, "Min = %.3lf\n", spec.MinValue());
    fprintf(f, "MinChannel = %ld\n", extendedInformation.MinChannel >= 0 ? extendedInformation.MinChannel : 0);
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

