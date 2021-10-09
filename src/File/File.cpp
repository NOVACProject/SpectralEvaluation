#include <SpectralEvaluation/StringUtils.h>
#include <SpectralEvaluation/Interpolation.h>
#include <SpectralEvaluation/File/File.h>
#include <SpectralEvaluation/File/STDFile.h>
#include <SpectralEvaluation/File/TXTFile.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Calibration/InstrumentCalibration.h>
#include <SpectralEvaluation/Calibration/InstrumentLineShape.h>
#include <SpectralEvaluation/Interpolation.h>
#include <SpectralEvaluation/VectorUtils.h>
#include <cstring>
#include <cmath>
#include <limits>
#include <memory>
#include <sstream>
#include <stdio.h>
#include <fstream>

#include <iostream>

namespace novac
{

std::pair<double, double> GetFwhm(const std::vector<double>& lineShape);

bool IsExistingFile(const std::string& fullFileName)
{
    FILE* f = fopen(fullFileName.c_str(), "r");
    if (f == nullptr)
    {
        return false;
    }
    fclose(f);

    return true;
}

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

    const bool twoColumns = spectrum.m_wavelength.size() == static_cast<size_t>(spectrum.m_length);
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

std::string GetFileExtension(const std::string& fullFilePath)
{
    if (fullFilePath.size() == 0)
    {
        return fullFilePath;
    }

    const size_t lastPeriod = fullFilePath.rfind('.');
    if (lastPeriod == fullFilePath.npos)
    {
        // no period found, hence no suffix
        return std::string();
    }

    const size_t lastForwardSlash = fullFilePath.rfind('\\');
    const size_t lastBackwardSlash = fullFilePath.rfind('/');

    if ((lastForwardSlash != fullFilePath.npos && lastPeriod > lastForwardSlash) ||
        (lastBackwardSlash != fullFilePath.npos && lastPeriod > lastBackwardSlash))
    {
        return fullFilePath.substr(lastPeriod, fullFilePath.size() - lastPeriod);
    }
    else
    {
        // no period found _after_ the last path-separator character, hence no suffix
        return std::string();
    }
}

size_t WavelengthToPixel(const std::vector<double>& pixelToWavelengthMapping, double wavelength)
{
    if (wavelength <= 0.0)
    {
        return 0;
    }
    else if (wavelength > pixelToWavelengthMapping.back())
    {
        return pixelToWavelengthMapping.size() - 1;
    }
    else
    {
        const double instrumentLineShapeFractionalCenterPixel = novac::GetFractionalIndex(pixelToWavelengthMapping, wavelength);
        return static_cast<size_t>(std::round(instrumentLineShapeFractionalCenterPixel));
    }
}

std::pair<std::string, std::string> FormatProperty(const char* name, double value)
{
    char formattedValue[128];
    sprintf(formattedValue, "%.9g", value);

    return std::make_pair(std::string(name), std::string(formattedValue));
}

bool TryGetParameter(const std::vector<std::pair<std::string, std::string>>& properties, const char* parameterName, double& parameterValue)
{
    for (const auto& p : properties)
    {
        if (p.first.compare(parameterName) == 0)
        {
            // found the correct parameter, try parse the value
            return 1 == sscanf(p.second.c_str(), "%lf", &parameterValue);
        }
    }

    return false;
}

bool SaveInstrumentCalibration(const std::string& fullFilePath, const InstrumentCalibration& calibration)
{
    if (calibration.pixelToWavelengthMapping.size() == 0)
    {
        return false;
    }
    if (calibration.instrumentLineShape.size() == 0 ||
        calibration.instrumentLineShapeGrid.size() == 0 ||
        calibration.instrumentLineShape.size() != calibration.instrumentLineShapeGrid.size())
    {
        return false;
    }

    novac::CSTDFile::ExtendedFormatInformation extendedFileInfo;

    size_t instrumentLineShapeSpectrumStartIdx = 0; // The pixel where the instrument line shape should start
    {
        size_t instrumentLineShapeCenterPixel = WavelengthToPixel(calibration.pixelToWavelengthMapping, calibration.instrumentLineShapeCenter);

        if (instrumentLineShapeCenterPixel <= calibration.instrumentLineShape.size() / 2)
        {
            instrumentLineShapeSpectrumStartIdx = 0;
        }
        else if (instrumentLineShapeCenterPixel >= calibration.pixelToWavelengthMapping.size() - calibration.instrumentLineShape.size())
        {
            instrumentLineShapeSpectrumStartIdx = calibration.pixelToWavelengthMapping.size() - calibration.instrumentLineShape.size();
        }
        else
        {
            instrumentLineShapeSpectrumStartIdx = instrumentLineShapeCenterPixel - calibration.instrumentLineShape.size() / 2;
        }
    }

    extendedFileInfo.MinChannel = static_cast<int>(instrumentLineShapeSpectrumStartIdx);
    extendedFileInfo.MaxChannel = static_cast<int>(instrumentLineShapeSpectrumStartIdx + calibration.instrumentLineShape.size());
    extendedFileInfo.MathLow = extendedFileInfo.MinChannel;
    extendedFileInfo.MathHigh = extendedFileInfo.MaxChannel;
    extendedFileInfo.Marker = calibration.instrumentLineShapeCenter;

    // If there is a parameterized instrument line shape function
    if (calibration.instrumentLineShapeParameter != nullptr)
    {
        const novac::SuperGaussianLineShape* superGaussian = dynamic_cast<novac::SuperGaussianLineShape*>(calibration.instrumentLineShapeParameter);
        if (superGaussian != nullptr)
        {
            extendedFileInfo.additionalProperties.push_back(FormatProperty("SuperGaussFitWidth", superGaussian->w));
            extendedFileInfo.additionalProperties.push_back(FormatProperty("SuperGaussFitPower", superGaussian->k));
        }
        else
        {
            const novac::GaussianLineShape* gaussian = dynamic_cast<novac::GaussianLineShape*>(calibration.instrumentLineShapeParameter);
            if (gaussian != nullptr)
            {
                extendedFileInfo.additionalProperties.push_back(FormatProperty("GaussFitSigma", gaussian->sigma));
            }
        }
    }

    // Additional information about the spectrum
    extendedFileInfo.calibrationPolynomial = calibration.pixelToWavelengthPolynomial;

    // Create the spectrum to save from the instrument line-shape + the 
    std::unique_ptr<novac::CSpectrum> spectrumToSave;
    {
        // Extend the measured spectrum line shape into a full spectrum
        std::vector<double> extendedPeak(calibration.pixelToWavelengthMapping.size(), 0.0);
        std::copy(calibration.instrumentLineShape.data(), calibration.instrumentLineShape.data() + calibration.instrumentLineShape.size(), extendedPeak.begin() + instrumentLineShapeSpectrumStartIdx);

        spectrumToSave = std::make_unique<novac::CSpectrum>(calibration.pixelToWavelengthMapping, extendedPeak);
    }

    // spectrumToSave->m_info = this->m_inputspectrumInformation;

    novac::CSTDFile::WriteSpectrum(*spectrumToSave, fullFilePath, extendedFileInfo);

    return true;
}

bool ReadInstrumentCalibration(const std::string& fullFilePath, CSpectrum& instrumentLineShape, std::vector<double>& pixelToWavelengthMapping)
{
    CSTDFile::ExtendedFormatInformation extendedFormatInformation;

    CSpectrum tmpSpectrum;
    if (!CSTDFile::ReadSpectrum(tmpSpectrum, fullFilePath, extendedFormatInformation))
    {
        return false;
    }

    if (extendedFormatInformation.MaxChannel <= extendedFormatInformation.MinChannel ||
        tmpSpectrum.m_wavelength.size() != static_cast<size_t>(tmpSpectrum.m_length))
    {
        // not a correct format of the data
        return false;
    }

    // assign the output data
    pixelToWavelengthMapping = tmpSpectrum.m_wavelength;
    instrumentLineShape.m_info = tmpSpectrum.m_info;
    instrumentLineShape.m_length = extendedFormatInformation.MaxChannel - extendedFormatInformation.MinChannel;
    memcpy(instrumentLineShape.m_data, tmpSpectrum.m_data + extendedFormatInformation.MinChannel, instrumentLineShape.m_length * sizeof(double));

    // Differentiate the wavelength wrt the peak
    instrumentLineShape.m_wavelength = std::vector<double>(begin(tmpSpectrum.m_wavelength) + extendedFormatInformation.MinChannel, begin(tmpSpectrum.m_wavelength) + extendedFormatInformation.MaxChannel);

    std::vector<double> specData{ instrumentLineShape.m_data, instrumentLineShape.m_data + instrumentLineShape.m_length };
    const double centerPixel = Centroid(specData);

    double centerWavelength = 0.0;
    if (!novac::LinearInterpolation(instrumentLineShape.m_wavelength, centerPixel, centerWavelength))
    {
        return false;
    }
    for (double& lambda : instrumentLineShape.m_wavelength)
    {
        lambda -= centerWavelength;
    }

    return true;
}

bool ReadInstrumentCalibration(const std::string& fullFilePath, InstrumentCalibration& result)
{
    result.Clear();

    CSpectrum tmpSpectrum;
    CSTDFile::ExtendedFormatInformation extendedFormatInformation;
    if (!CSTDFile::ReadSpectrum(tmpSpectrum, fullFilePath, extendedFormatInformation))
    {
        return false;
    }

    if (extendedFormatInformation.MaxChannel <= extendedFormatInformation.MinChannel ||
        tmpSpectrum.m_wavelength.size() != static_cast<size_t>(tmpSpectrum.m_length))
    {
        // not a correct format of the data
        return false;
    }

    // The pixel to wavelength mapping
    result.pixelToWavelengthMapping = tmpSpectrum.m_wavelength;
    result.pixelToWavelengthPolynomial = extendedFormatInformation.calibrationPolynomial;

    // Copy out the instrument line shape
    std::vector<double> instrumentLineShape;
    instrumentLineShape.resize(static_cast<size_t>(extendedFormatInformation.MaxChannel - extendedFormatInformation.MinChannel));
    memcpy(instrumentLineShape.data(), tmpSpectrum.m_data + extendedFormatInformation.MinChannel, instrumentLineShape.size() * sizeof(double));

    const double instrumentLineShapeCenterPixel = Centroid(instrumentLineShape);
    const std::pair<double, double> leftAndRightFwhmIdx = GetFwhm(instrumentLineShape);

    // Check if the instrument line shape is reasonable.
    if (leftAndRightFwhmIdx.first < (double)extendedFormatInformation.MinChannel ||
        leftAndRightFwhmIdx.second >(double)extendedFormatInformation.MaxChannel ||
        instrumentLineShapeCenterPixel < leftAndRightFwhmIdx.first ||
        instrumentLineShapeCenterPixel > leftAndRightFwhmIdx.second)
    {
        // We did at least manage to read the pixel-to-wavelength mapping.
        return true;
    }

    double wavelengthAtLeftFwhm = 0.0;
    double wavelengthAtRightFwhm = 0.0;
    if (!LinearInterpolation(result.pixelToWavelengthMapping, leftAndRightFwhmIdx.first, wavelengthAtLeftFwhm) ||
        !LinearInterpolation(result.pixelToWavelengthMapping, leftAndRightFwhmIdx.second, wavelengthAtRightFwhm))
    {
        // Interpolation failure, instrument line shape seems unreasonable.
        return true;
    }

    if (std::abs(wavelengthAtRightFwhm - wavelengthAtLeftFwhm) > 0.01 * std::abs(result.pixelToWavelengthMapping.back() - result.pixelToWavelengthMapping.front()))
    {
        // Instrument line shape has a too large full width at half maximum.
        return true;
    }

    result.instrumentLineShape = instrumentLineShape;

    // Get the instrument line shape grid and make sure that it is differentiated wrt the peak
    {
        result.instrumentLineShapeGrid = std::vector<double>(begin(tmpSpectrum.m_wavelength) + extendedFormatInformation.MinChannel, begin(tmpSpectrum.m_wavelength) + extendedFormatInformation.MaxChannel);

        double centerWavelength = 0.0;
        if (!novac::LinearInterpolation(result.instrumentLineShapeGrid, instrumentLineShapeCenterPixel, centerWavelength))
        {
            return false;
        }

        result.instrumentLineShapeCenter = centerWavelength;

        if (std::abs(centerWavelength) > std::numeric_limits<double>::epsilon())
        {
            for (double& lambda : result.instrumentLineShapeGrid)
            {
                lambda -= centerWavelength;
            }
        }
    }

    // Handle the reading of an parametrization of the instrument line shape.
    {
        double superGaussianWidth = 0.0;
        double superGaussianPower = 0.0;
        if (TryGetParameter(extendedFormatInformation.additionalProperties, "SuperGaussFitWidth", superGaussianWidth) &&
            TryGetParameter(extendedFormatInformation.additionalProperties, "SuperGaussFitPower", superGaussianPower))
        {
            result.instrumentLineShapeParameter = new novac::SuperGaussianLineShape(superGaussianWidth, superGaussianPower);
        }
        else
        {
            double gaussSigma = 0.0;
            if (TryGetParameter(extendedFormatInformation.additionalProperties, "GaussFitSigma", gaussSigma))
            {
                result.instrumentLineShapeParameter = new novac::GaussianLineShape(gaussSigma);
            }
        }
    }

    return true;
}

}