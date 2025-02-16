#undef max
#undef min

#include <sstream>
#include <chrono>

#include <SpectralEvaluation/DialogControllers/WavelengthCalibrationController.h>
#include <SpectralEvaluation/VectorUtils.h>
#include <SpectralEvaluation/StringUtils.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Spectra/SpectrometerModel.h>
#include <SpectralEvaluation/File/File.h>
#include <SpectralEvaluation/Calibration/Correspondence.h>
#include <SpectralEvaluation/Calibration/InstrumentCalibration.h>
#include <SpectralEvaluation/Calibration/WavelengthCalibration.h>
#include <SpectralEvaluation/Calibration/InstrumentLineShapeEstimation.h>
#include <SpectralEvaluation/Calibration/FraunhoferSpectrumGeneration.h>
#include <SpectralEvaluation/Calibration/WavelengthCalibrationByRansac.h>

std::string ToString(double value)
{
    char buffer[64];
    sprintf(buffer, "%0.3lf", value);
    return std::string(buffer);
}

std::vector<std::pair<std::string, std::string>> GetFunctionDescription(const novac::ParametricInstrumentLineShape* lineShapeFunction)
{
    std::vector<std::pair<std::string, std::string>> result;

    if (lineShapeFunction == nullptr)
    {
        return result;
    }
    else if (lineShapeFunction->Type() == novac::InstrumentLineShapeFunctionType::Gaussian)
    {
        const auto func = static_cast<const novac::GaussianLineShape*>(lineShapeFunction);
        result.push_back(std::make_pair("type", "Gaussian"));
        result.push_back(std::make_pair("sigma", ToString(func->sigma)));
        result.push_back(std::make_pair("fwhm", ToString(func->Fwhm())));
    }
    else if (lineShapeFunction->Type() == novac::InstrumentLineShapeFunctionType::SuperGaussian)
    {
        const auto func = static_cast<const novac::SuperGaussianLineShape*>(lineShapeFunction);
        result.push_back(std::make_pair("type", "Super Gaussian"));
        result.push_back(std::make_pair("w", ToString(func->w)));
        result.push_back(std::make_pair("k", ToString(func->k)));
        result.push_back(std::make_pair("fwhm", ToString(func->Fwhm())));
    }
    else
    {
        // unknown type...
    }

    return result;
}

WavelengthCalibrationController::WavelengthCalibrationController(novac::ILogger& log)
    : m_instrumentLineShapeFitOption(InstrumentLineShapeFitOption::None),
    m_instrumentLineShapeFitRegion(330.0, 350.0),
    m_calibrationDebug(0U),
    m_log(log)
{
}

WavelengthCalibrationController::~WavelengthCalibrationController()
{
    m_initialCalibration.release();
    m_resultingCalibration.release();
}

/// <summary>
/// Reads and parses m_initialLineShapeFile and saves the result to the provided settings
/// </summary>
void ReadInstrumentLineShape(const std::string& initialLineShapeFile, novac::InstrumentCalibration& calibration)
{
    novac::CCrossSectionData measuredInstrumentLineShape;
    if (!novac::ReadCrossSectionFile(initialLineShapeFile, measuredInstrumentLineShape))
    {
        throw std::invalid_argument("Cannot read the provided instrument lineshape file");
    }

    calibration.instrumentLineShape = measuredInstrumentLineShape.m_crossSection;
    calibration.instrumentLineShapeGrid = measuredInstrumentLineShape.m_waveLength;
}

void WavelengthCalibrationController::CreateGuessForInstrumentLineShape(const std::string& solarSpectrumFile, const novac::CSpectrum& measuredSpectrum, novac::InstrumentCalibration& calibration)
{
    // The user has not supplied an instrument-line-shape, create a guess for one.
    std::vector<std::pair<std::string, double>> noCrossSections;
    novac::FraunhoferSpectrumGeneration fraunhoferSpectrumGen{ solarSpectrumFile, noCrossSections };

    novac::InstrumentLineShapeEstimationFromKeypointDistance ilsEstimation{ calibration.pixelToWavelengthMapping };

    double resultInstrumentFwhm;
    novac::CCrossSectionData estimatedInstrumentLineShape;
    ilsEstimation.EstimateInstrumentLineShape(fraunhoferSpectrumGen, measuredSpectrum, estimatedInstrumentLineShape, resultInstrumentFwhm);

    Log("No initial instrument line shape could be found, estimated instrument line shape as Gaussian with fwhm of: ", resultInstrumentFwhm);

    calibration.instrumentLineShape = estimatedInstrumentLineShape.m_crossSection;
    calibration.instrumentLineShapeGrid = estimatedInstrumentLineShape.m_waveLength;
}


novac::SpectrometerModel WavelengthCalibrationController::GetModelForMeasurement(const std::string& deviceSerial) const
{
    if (m_spectrometerModel != nullptr)
    {
        return *(m_spectrometerModel);
    }

    auto model = novac::CSpectrometerDatabase::GetInstance().GuessModelFromSerial(deviceSerial);
    model.averagesSpectra = m_spectraAreAverages;

    return model;
}

void WavelengthCalibrationController::RunCalibration()
{
    m_errorMessage.clear();
    m_logMessages.clear();
    m_instrumentLineShapeParameterDescriptions.clear();

    auto startTime = std::chrono::steady_clock::now();

    novac::WavelengthCalibrationSettings settings;
    settings.highResSolarAtlas = m_solarSpectrumFile;

    novac::CSpectrum measuredSpectrum;
    ReadInput(measuredSpectrum);

    // Copy out the spectrum, such that the user can see it.
    m_calibrationDebug.measuredSpectrum = std::vector<double>(measuredSpectrum.m_data, measuredSpectrum.m_data + measuredSpectrum.m_length);

    auto spectrometerModel = GetModelForMeasurement(measuredSpectrum.m_info.m_device);

    // Check that this is a good measurement, according to our standards.
    CheckSpectrumQuality(measuredSpectrum, spectrometerModel);

    // Read the initial calibration
    m_initialCalibration = std::make_unique<novac::InstrumentCalibration>();
    if (EqualsIgnoringCase(novac::GetFileExtension(m_initialCalibrationFile), ".std"))
    {
        // This is a file which may contain either just the wavelength calibration _or_ both a wavelength calibration and an instrument line shape.
        if (!novac::ReadInstrumentCalibration(m_initialCalibrationFile, *m_initialCalibration))
        {
            throw std::invalid_argument("Failed to read the instrument calibration file");
        }

        Log("Read initial instrument calibration from: ", m_initialCalibrationFile);

        if (m_initialCalibration->instrumentLineShapeGrid.size() > 0)
        {
            const double fwhm = novac::GetFwhm(m_initialCalibration->instrumentLineShapeGrid, m_initialCalibration->instrumentLineShape);
            Log("Initial instrument line shape has a fwhm of: ", fwhm);
        }
        else
        {
            Log("Provided instrument calibration file does not contain an instrument line shape.");
        }
    }
    else
    {
        m_initialCalibration->pixelToWavelengthMapping = novac::GetPixelToWavelengthMappingFromFile(m_initialCalibrationFile);
        if (!novac::IsPossiblePixelToWavelengthCalibration(m_initialCalibration->pixelToWavelengthMapping))
        {
            throw std::invalid_argument("Error interpreting the provided initial pixel to wavelength calibration, the wavelengths are not monotonically increasing.");
        }

        Log("Read initial pixel to wavelength calibration from: ", m_initialCalibrationFile);
    }

    // Check the initial instrument line shape.
    if (m_initialLineShapeFile.size() > 0)
    {
        ReadInstrumentLineShape(m_initialLineShapeFile, *m_initialCalibration);
        Log("Read initial instrument line shape from: ", m_initialLineShapeFile);

        const double fwhm = novac::GetFwhm(m_initialCalibration->instrumentLineShapeGrid, m_initialCalibration->instrumentLineShape);
        Log("Initial instrument line shape has a fwhm of :", fwhm);
    }
    else if (m_initialCalibration->instrumentLineShape.size() == 0)
    {
        CreateGuessForInstrumentLineShape(m_solarSpectrumFile, measuredSpectrum, *m_initialCalibration);
    }

    settings.initialPixelToWavelengthMapping = m_initialCalibration->pixelToWavelengthMapping;
    settings.initialInstrumentLineShape.m_crossSection = m_initialCalibration->instrumentLineShape;
    settings.initialInstrumentLineShape.m_waveLength = m_initialCalibration->instrumentLineShapeGrid;

    // Normalize the initial instrument line shape
    Normalize(m_initialCalibration->instrumentLineShape);

    if (m_instrumentLineShapeFitOption == WavelengthCalibrationController::InstrumentLineShapeFitOption::SuperGaussian)
    {
        if (m_instrumentLineShapeFitRegion.low > m_instrumentLineShapeFitRegion.high)
        {
            std::stringstream msg;
            msg << "Invalid region for fitting the instrument line shape ";
            msg << "(" << m_instrumentLineShapeFitRegion.low << ", " << m_instrumentLineShapeFitRegion.high << ") [nm]. ";
            msg << "From must be smaller than To";
            throw std::invalid_argument(msg.str());
        }
        if (m_instrumentLineShapeFitRegion.low < settings.initialPixelToWavelengthMapping.front() ||
            m_instrumentLineShapeFitRegion.high > settings.initialPixelToWavelengthMapping.back())
        {
            std::stringstream msg;
            msg << "Invalid region for fitting the instrument line shape ";
            msg << "(" << m_instrumentLineShapeFitRegion.low << ", " << m_instrumentLineShapeFitRegion.high << ") [nm]. ";
            msg << "Region does not overlap initial pixel to wavelength calibration: ";
            msg << "(" << settings.initialPixelToWavelengthMapping.front() << ", " << settings.initialPixelToWavelengthMapping.back() << ") [nm]. ";
            throw std::invalid_argument(msg.str());
        }

        settings.estimateInstrumentLineShape = novac::InstrumentLineshapeEstimationOption::SuperGaussian;
        settings.estimateInstrumentLineShapeWavelengthRegion.first = m_instrumentLineShapeFitRegion.low;
        settings.estimateInstrumentLineShapeWavelengthRegion.second = m_instrumentLineShapeFitRegion.high;

        settings.crossSectionsForInstrumentLineShapeFitting.clear();
        if (m_crossSectionsForInstrumentLineShapeFitting.size() > 0)
        {
            settings.crossSectionsForInstrumentLineShapeFitting.push_back(m_crossSectionsForInstrumentLineShapeFitting);
        }
    }
    else
    {
        settings.estimateInstrumentLineShape = novac::InstrumentLineshapeEstimationOption::None;
    }

    // Clear the previous result
    m_resultingCalibration = std::make_unique<novac::InstrumentCalibration>();

    novac::WavelengthCalibrationSetup setup{ settings };
    auto result = setup.DoWavelengthCalibration(measuredSpectrum);

    // Copy out the result
    m_resultingCalibration->pixelToWavelengthMapping = result.pixelToWavelengthMapping;
    m_resultingCalibration->pixelToWavelengthPolynomial = result.pixelToWavelengthMappingCoefficients;

    if (result.estimatedInstrumentLineShape.GetSize() > 0)
    {
        m_resultingCalibration->instrumentLineShape = result.estimatedInstrumentLineShape.m_crossSection;
        m_resultingCalibration->instrumentLineShapeGrid = result.estimatedInstrumentLineShape.m_waveLength;
        m_resultingCalibration->instrumentLineShapeCenter = 0.5 * (m_instrumentLineShapeFitRegion.low + m_instrumentLineShapeFitRegion.high);
    }

    if (result.estimatedInstrumentLineShapeParameters != nullptr)
    {
        m_instrumentLineShapeParameterDescriptions = GetFunctionDescription(&(*result.estimatedInstrumentLineShapeParameters));
        m_resultingCalibration->instrumentLineShapeParameter = result.estimatedInstrumentLineShapeParameters->Clone();
    }

    // Also copy out some debug information, which makes it possible for the user to inspect the calibration
    {
        const auto& calibrationDebug = setup.GetLastCalibrationSetup();
        m_calibrationDebug = WavelengthCalibrationDebugState(calibrationDebug.allCorrespondences.size());
        m_calibrationDebug.initialPixelToWavelengthMapping = settings.initialPixelToWavelengthMapping;
        m_calibrationDebug.spectrumInfo = measuredSpectrum.m_info;

        for (size_t correspondenceIdx = 0; correspondenceIdx < calibrationDebug.allCorrespondences.size(); ++correspondenceIdx)
        {
            const auto& corr = calibrationDebug.allCorrespondences[correspondenceIdx]; // easier syntax.

            if (calibrationDebug.correspondenceIsInlier[correspondenceIdx])
            {
                m_calibrationDebug.inlierCorrespondencePixels.push_back(corr.measuredValue);
                m_calibrationDebug.inlierCorrespondenceWavelengths.push_back(corr.theoreticalValue);
                m_calibrationDebug.inlierCorrespondenceMeasuredIntensity.push_back(calibrationDebug.measuredKeypoints[corr.measuredIdx].intensity);
                m_calibrationDebug.inlierCorrespondenceFraunhoferIntensity.push_back(calibrationDebug.fraunhoferKeypoints[corr.theoreticalIdx].intensity);

                m_calibrationDebug.measuredSpectrumInlierKeypointPixels.push_back(calibrationDebug.measuredKeypoints[corr.measuredIdx].pixel);
                m_calibrationDebug.measuredSpectrumInlierKeypointIntensities.push_back(calibrationDebug.measuredKeypoints[corr.measuredIdx].intensity);

                m_calibrationDebug.fraunhoferSpectrumInlierKeypointWavelength.push_back(calibrationDebug.fraunhoferKeypoints[corr.theoreticalIdx].wavelength);
                m_calibrationDebug.fraunhoferSpectrumInlierKeypointIntensities.push_back(calibrationDebug.fraunhoferKeypoints[corr.theoreticalIdx].intensity);
            }
            else
            {
                m_calibrationDebug.outlierCorrespondencePixels.push_back(corr.measuredValue);
                m_calibrationDebug.outlierCorrespondenceWavelengths.push_back(corr.theoreticalValue);
            }
        }

        m_calibrationDebug.measuredSpectrum = std::vector<double>(calibrationDebug.measuredSpectrum->m_data, calibrationDebug.measuredSpectrum->m_data + calibrationDebug.measuredSpectrum->m_length);
        m_calibrationDebug.fraunhoferSpectrum = std::vector<double>(calibrationDebug.fraunhoferSpectrum->m_data, calibrationDebug.fraunhoferSpectrum->m_data + calibrationDebug.fraunhoferSpectrum->m_length);

        for (const auto& pt : calibrationDebug.measuredKeypoints)
        {
            m_calibrationDebug.measuredSpectrumKeypointPixels.push_back(pt.pixel);
            m_calibrationDebug.measuredSpectrumKeypointIntensities.push_back(pt.intensity);
        }

        for (const auto& pt : calibrationDebug.fraunhoferKeypoints)
        {
            m_calibrationDebug.fraunhoferSpectrumKeypointWavelength.push_back(pt.wavelength);
            m_calibrationDebug.fraunhoferSpectrumKeypointIntensities.push_back(pt.intensity);
        }
    }

    auto stopTime = std::chrono::steady_clock::now();
    {
        std::stringstream message;
        message << "Instrument calibration done in " << std::chrono::duration_cast<std::chrono::seconds>(stopTime - startTime).count() << " s";
        Log(message.str());
    }
}

std::unique_ptr<novac::InstrumentCalibration> WavelengthCalibrationController::GetFinalCalibration() const
{
    std::unique_ptr<novac::InstrumentCalibration> output;

    if (m_resultingCalibration->instrumentLineShape.size() > 0)
    {
        // The created calibration is complete, copy and return
        output = std::make_unique<novac::InstrumentCalibration>(*m_resultingCalibration);
    }
    else
    {
        // Create an instrument calibration, taking the instrument line shape from the initial setup 
        //  and the pixel-to-wavelength mapping from the calibration result
        output = std::make_unique<novac::InstrumentCalibration>();

        output->pixelToWavelengthMapping = m_resultingCalibration->pixelToWavelengthMapping;
        output->pixelToWavelengthPolynomial = m_resultingCalibration->pixelToWavelengthPolynomial;

        output->instrumentLineShape = m_initialCalibration->instrumentLineShape;
        output->instrumentLineShapeGrid = m_initialCalibration->instrumentLineShapeGrid;
        output->instrumentLineShapeCenter = m_initialCalibration->instrumentLineShapeCenter;

        if (m_initialCalibration->instrumentLineShapeParameter != nullptr)
        {
            output->instrumentLineShapeParameter = m_initialCalibration->instrumentLineShapeParameter->Clone();
        }
    }

    return output;
}

void WavelengthCalibrationController::SaveResultAsStd(const std::string& filename)
{
    if (m_resultingCalibration->instrumentLineShape.size() > 0)
    {
        // We have fitted both an instrument line shape and a pixel-to-wavelength mapping.
        //  I.e. m_resultingCalibration contains a full calibration to be saved.
        if (!novac::SaveInstrumentCalibration(filename, *m_resultingCalibration))
        {
            throw std::invalid_argument("Failed to save the resulting instrument calibration file");
        }
        return;
    }

    const auto mixedCalibration = GetFinalCalibration();
    if (!novac::SaveInstrumentCalibration(filename, *mixedCalibration))
    {
        throw std::invalid_argument("Failed to save the resulting instrument calibration file");
    }
}

void WavelengthCalibrationController::SaveResultAsClb(const std::string& filename)
{
    novac::SaveDataToFile(filename, m_resultingCalibration->pixelToWavelengthMapping);
}

void WavelengthCalibrationController::SaveResultAsSlf(const std::string& filename)
{
    novac::CCrossSectionData instrumentLineShape;
    instrumentLineShape.m_crossSection = m_resultingCalibration->instrumentLineShape;
    instrumentLineShape.m_waveLength = m_resultingCalibration->instrumentLineShapeGrid;

    novac::SaveCrossSectionFile(filename, instrumentLineShape);
}

double WavelengthCalibrationController::GetSpectrometerMaxIntensityForSingleReadout(const novac::CSpectrum& spectrum, const novac::SpectrometerModel& model) const
{
    if (m_spectrometerMaximumIntensityForSingleReadout > 0.0)
    {
        // The dynamic range of the spectrometer is explicitly set by the user
        return novac::GetMaximumSaturationRatioOfSpectrum(spectrum, m_spectrometerMaximumIntensityForSingleReadout);
    }

    // If the device is unknown, then the dynamic range is also unknown.
    if (model.modelName == "UNKNOWN")
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    // Derive the dynamic range from the spectrometer model.
    const double maximumIntensityForSingleReadout = model.FullDynamicRangeForSpectrum(spectrum.m_info);
    return novac::GetMaximumSaturationRatioOfSpectrum(spectrum, maximumIntensityForSingleReadout);
}

void WavelengthCalibrationController::CheckSpectrumQuality(const novac::CSpectrum& spectrum, const novac::SpectrometerModel& spectrometerModel) const
{
    const double maximumSaturationRatio = GetSpectrometerMaxIntensityForSingleReadout(spectrum, spectrometerModel);

    if (std::isnan(maximumSaturationRatio))
    {
        // unknown spectrometer model.
        return;
    }

    if (maximumSaturationRatio > 0.85)
    {
        std::stringstream message;
        message << "The provided sky spectrum seems to be saturated ";
        message << "(maximum intensity: " << spectrum.MaxValue(0, spectrum.m_length) << ", corresponding to : " << 100 * maximumSaturationRatio << "% of full range for a " << spectrometerModel.modelName << ").";
        message << "Calibration aborted";
        throw std::invalid_argument(message.str());
    }

    if (maximumSaturationRatio < 0.20)
    {
        std::stringstream message;
        message << "The provided sky spectrum seems to be too dark for the calibration to succeed ";
        message << "(maximum intensity: " << spectrum.MaxValue(0, spectrum.m_length) << ", corresponding to : " << 100 * maximumSaturationRatio << "% of full range for a " << spectrometerModel.modelName << ").";
        message << "Calibration aborted";
        throw std::invalid_argument(message.str());
    }
}

void WavelengthCalibrationController::ClearResult()
{
    m_calibrationDebug = WavelengthCalibrationController::WavelengthCalibrationDebugState(0U);
    m_errorMessage.clear();
    m_logMessages.clear();
    m_instrumentLineShapeParameterDescriptions.clear();
    m_resultingCalibration.reset();
    m_initialCalibration.reset();
}

void WavelengthCalibrationController::Log(const std::string& message)
{
    m_logMessages.push_back(message);
}

void WavelengthCalibrationController::Log(const std::string& message, double value)
{
    std::stringstream formattedMessage;
    formattedMessage << message << value;
    m_logMessages.push_back(formattedMessage.str());
}

void WavelengthCalibrationController::Log(const std::string& part1, const std::string& part2)
{
    std::string message = part1;
    message.append(part2);
    m_logMessages.push_back(message);
}

void MobileDoasWavelengthCalibrationController::ReadInput(novac::CSpectrum& measuredSpectrum)
{
    if (!novac::ReadSpectrum(m_inputSpectrumFile, measuredSpectrum))
    {
        throw std::invalid_argument("Cannot read the provided input spectrum file");
    }
    Log("Read measured spectrum: ", m_inputSpectrumFile);

    // subtract the dark-spectrum (if this has been provided)
    novac::CSpectrum darkSpectrum;
    if (novac::ReadSpectrum(m_darkSpectrumFile, darkSpectrum))
    {
        measuredSpectrum.Sub(darkSpectrum);
        Log("Read and subtracted dark spectrum: ", m_darkSpectrumFile);
    }
}

void InMemoryWavelengthCalibrationController::ReadInput(novac::CSpectrum& measuredSpectrum)
{
    measuredSpectrum = m_measuredSpectrum;
}

