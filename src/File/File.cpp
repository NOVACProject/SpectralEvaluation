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
#include <assert.h>
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

    // General functions declared elsewhere.
    std::pair<double, double> GetFwhm(const std::vector<double>& lineShape);
    size_t WavelengthToPixel(const std::vector<double>& pixelToWavelengthMapping, double wavelength);;

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

        const size_t lastBackwardSlash = fullFilePath.rfind('\\');
        const size_t lastForwardSlash = fullFilePath.rfind('/');

        if ((lastBackwardSlash != fullFilePath.npos && lastPeriod > lastBackwardSlash) ||
            (lastForwardSlash != fullFilePath.npos && lastPeriod > lastForwardSlash))
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
            throw std::invalid_argument("Cannot save an instrument calibration with an empty pixel to wavelength mapping.");
        }
        if (calibration.instrumentLineShape.size() == 0 ||
            calibration.instrumentLineShapeGrid.size() == 0 ||
            calibration.instrumentLineShape.size() != calibration.instrumentLineShapeGrid.size())
        {
            throw std::invalid_argument("Cannot save an instrument calibration with an empty instrument line shape description.");
        }
        if (calibration.instrumentLineShapeGrid.back() < calibration.instrumentLineShapeGrid.front())
        {
            throw std::invalid_argument("Cannot save an instrument calibration where the instrument line shape is not sampled on an increasing grid.");
        }

        const double instrumentLineShapeGridRange = calibration.instrumentLineShapeGrid.back() - calibration.instrumentLineShapeGrid.front();
        const double instrumentWavelengthRange = calibration.pixelToWavelengthMapping.back() - calibration.pixelToWavelengthMapping.front();
        if (instrumentLineShapeGridRange > instrumentWavelengthRange)
        {
            throw std::invalid_argument("Cannot save an instrument calibration where the instrument line shape is wider than the wavelength range.");
        }
        if (calibration.instrumentLineShapeGrid.front() > 0.0)
        {
            // If the instrument line shape grid is absolute, then it must lie inside the wavelength range of the device
            if (calibration.instrumentLineShapeGrid.front() < calibration.pixelToWavelengthMapping.front() ||
                calibration.instrumentLineShapeGrid.back() > calibration.pixelToWavelengthMapping.back())
            {
                throw std::invalid_argument("Cannot save an instrument calibration where the instrument line shape lies outside of the device's wavelength range.");
            }
        }

        novac::CSTDFile::ExtendedFormatInformation extendedFileInfo;

        // Make sure that the instrument line shape is sampled on the same grid as the pixel-to-wavelength mapping.

        // First calculate the first and the last (absolute) wavelengths of the current grid.
        //  Notice that the incoming calibration.instrumentLineShapeGrid _may_ be relative or absolute, use an offset to differentiate the two cases.
        const double instrumentLineShapeGridToAbsoluteWavelength = (calibration.instrumentLineShapeGrid.front() > 0.0) ? 0.0 : calibration.instrumentLineShapeCenter;
        double instrumentLineShapeStartWavelength = instrumentLineShapeGridToAbsoluteWavelength + calibration.instrumentLineShapeGrid.front();
        double instrumentLineShapeEndWavelength = instrumentLineShapeGridToAbsoluteWavelength + calibration.instrumentLineShapeGrid.back();

        std::vector<double> resampledLineShape;
        if (instrumentLineShapeStartWavelength < calibration.pixelToWavelengthMapping.front())
        {
            instrumentLineShapeStartWavelength = calibration.pixelToWavelengthMapping.front();
            instrumentLineShapeEndWavelength = instrumentLineShapeStartWavelength + instrumentLineShapeGridRange;
        }
        else if (instrumentLineShapeEndWavelength > calibration.pixelToWavelengthMapping.back())
        {
            instrumentLineShapeEndWavelength = calibration.pixelToWavelengthMapping.back();
            instrumentLineShapeStartWavelength = instrumentLineShapeEndWavelength - instrumentLineShapeGridRange;
        }

        // The pixels where the instrument line shape should start and end
        const size_t instrumentLineShapeSpectrumStartPixel = WavelengthToPixel(calibration.pixelToWavelengthMapping, instrumentLineShapeStartWavelength);
        const size_t instrumentLineShapeSpectrumEndPixel = 1 + WavelengthToPixel(calibration.pixelToWavelengthMapping, instrumentLineShapeEndWavelength);
        const size_t instrumentLineShapeCenterPixel = (instrumentLineShapeSpectrumEndPixel + instrumentLineShapeSpectrumStartPixel) / 2; // the middle pixel

        // Perform the resampling
        {
            std::vector<double> newLineShapeGrid;
            if (calibration.instrumentLineShapeGrid.front() > 0.0)
            {
                // Absolute sampling grid
                for (size_t ii = instrumentLineShapeSpectrumStartPixel; ii < instrumentLineShapeSpectrumEndPixel; ++ii)
                {
                    newLineShapeGrid.push_back(calibration.pixelToWavelengthMapping[ii]);
                }
            }
            else
            {
                // Relative sampling grid
                const double wavelengthAtInstrumentLineShapeCenterPixel = calibration.pixelToWavelengthMapping[instrumentLineShapeCenterPixel]; // due to rounding, this is not exactly same as instrumentLineShapeCenter
                for (size_t ii = instrumentLineShapeSpectrumStartPixel; ii < instrumentLineShapeSpectrumEndPixel; ++ii)
                {
                    newLineShapeGrid.push_back(calibration.pixelToWavelengthMapping[ii] - wavelengthAtInstrumentLineShapeCenterPixel);
                }
            }

            novac::Resample(calibration.instrumentLineShapeGrid, calibration.instrumentLineShape, newLineShapeGrid, resampledLineShape);
        }

        extendedFileInfo.MinChannel = static_cast<int>(instrumentLineShapeSpectrumStartPixel);
        extendedFileInfo.MaxChannel = static_cast<int>(instrumentLineShapeSpectrumEndPixel);
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
            std::copy(resampledLineShape.data(), resampledLineShape.data() + resampledLineShape.size(), extendedPeak.begin() + instrumentLineShapeSpectrumStartPixel);

            spectrumToSave = std::make_unique<novac::CSpectrum>(calibration.pixelToWavelengthMapping, extendedPeak);
        }

        return novac::CSTDFile::WriteSpectrum(*spectrumToSave, fullFilePath, extendedFileInfo);
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

        double wavelengthAtLeftFwhm = 0.0;
        double wavelengthAtRightFwhm = 0.0;
        if (std::abs(leftAndRightFwhmIdx.second - leftAndRightFwhmIdx.first) < 2.0 ||
            !LinearInterpolation(result.pixelToWavelengthMapping, leftAndRightFwhmIdx.first, wavelengthAtLeftFwhm) ||
            !LinearInterpolation(result.pixelToWavelengthMapping, leftAndRightFwhmIdx.second, wavelengthAtRightFwhm))
        {
            // Interpolation failure, instrument line shape seems unreasonable.
            return true;
        }

        if (std::abs(wavelengthAtRightFwhm - wavelengthAtLeftFwhm) > 0.01 * std::abs(result.pixelToWavelengthMapping.back() - result.pixelToWavelengthMapping.front()))
        {
            // Instrument line shape has a too large full width at half maximum.
            std::cout << "The provided .std file does not contain an instrument line shape. Estimated fwhm is: " << std::abs(wavelengthAtRightFwhm - wavelengthAtLeftFwhm);
            std::cout << " and entire wavelength range of spectrum is " << std::abs(result.pixelToWavelengthMapping.back() - result.pixelToWavelengthMapping.front()) << "nm " << std::endl;
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