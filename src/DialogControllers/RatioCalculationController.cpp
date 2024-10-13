#include <SpectralEvaluation/DialogControllers/RatioCalculationController.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/Evaluation/DoasFitPreparation.h>
#include <SpectralEvaluation/Evaluation/RatioEvaluation.h>
#include <SpectralEvaluation/Configuration/DarkSettings.h>
#include <SpectralEvaluation/File/XmlUtil.h>
#include <SpectralEvaluation/File/SpectrumIO.h>
#include <SpectralEvaluation/File/STDFile.h>
#include <SpectralEvaluation/File/File.h>
#include <SpectralEvaluation/Spectra/SpectrometerModel.h>

#include <fstream>
#include <sstream>
#include <cmath>

namespace novac
{
    std::string FormatDate(const CSpectrumInfo& spectrumInfo);
    std::string FormatTimestamp(const CDateTime& time);
}

RatioCalculationController::RatioCalculationController(novac::ILogger& log)
    : m_so2FitRange(314.8, 326.8), m_broFitRange(330.6, 352.8), m_log(log)
{
    InitializeToDefault();
}

void RatioCalculationController::InitializeToDefault()
{
    m_pakfiles.clear();

    // Default settings
    m_ratioEvaluationSettings.minimumPlumeCompleteness = 0.7;
    m_ratioEvaluationSettings.requireVisiblePlumeEdges = true;

    m_spectrometerModel.reset();

    // Insert the default species
    m_references.clear();
    m_references.push_back(ReferenceForRatioCalculation(StandardDoasSpecie::SO2, "SO2", "", true, true, false));
    m_references.push_back(ReferenceForRatioCalculation(StandardDoasSpecie::BRO, "BrO", "", false, true, false));
    m_references.push_back(ReferenceForRatioCalculation(StandardDoasSpecie::O3, "O3", "", true, true, false));
    m_references.push_back(ReferenceForRatioCalculation(StandardDoasSpecie::RING, "Ring", "", true, true, true));
    m_references.push_back(ReferenceForRatioCalculation(StandardDoasSpecie::RING_LAMBDA4, "Ringxlambda^4", "", true, true, true));
}

void RatioCalculationController::LoadSetup(const std::string& setupFilePath)
{
    try
    {
        // Super basic xml parsing
        std::ifstream file(setupFilePath, std::ios::in);
        std::string line;
        while (std::getline(file, line))
        {
            if (line.find("<Reference>") != std::string::npos)
            {
                auto name = novac::ParseXmlString("Name", line);
                auto path = novac::ParseXmlString("Path", line);
                auto inMajor = novac::ParseXmlString("IncludeInMajor", line);
                auto inMinor = novac::ParseXmlString("IncludeInMinor", line);
                auto doCalculate = novac::ParseXmlString("Calculate", line);

                ReferenceForRatioCalculation* ref = GetReferenceWithName(*this, name);
                if (ref != nullptr)
                {
                    ref->m_path = path;
                    ref->m_includeInMajor = (inMajor == "1");
                    ref->m_includeInMinor = (inMinor == "1");
                    ref->m_automaticallyCalculate = (doCalculate == "1");
                }
            }
            else if (line.find("<SO2Setup>") != std::string::npos)
            {
                auto fromStr = novac::ParseXmlString("From", line);
                auto toStr = novac::ParseXmlString("To", line);
                auto polyOrder = novac::ParseXmlInteger("Poly", line, -1);

                if (!fromStr.empty() && !toStr.empty())
                {
                    double from = std::atof(fromStr.c_str());
                    double to = std::atof(toStr.c_str());
                    if (from < to && from > 0 && to > 0)
                    {
                        m_so2FitRange = novac::WavelengthRange(from, to);
                    }
                }

                if (polyOrder >= 0)
                {
                    m_so2PolynomialOrder = polyOrder;
                }

            }
            else if (line.find("<BrOSetup>") != std::string::npos)
            {
                auto fromStr = novac::ParseXmlString("From", line);
                auto toStr = novac::ParseXmlString("To", line);
                auto polyOrder = novac::ParseXmlInteger("Poly", line, -1);

                if (!fromStr.empty() && !toStr.empty())
                {
                    double from = std::atof(fromStr.c_str());
                    double to = std::atof(toStr.c_str());
                    if (from < to && from > 0 && to > 0)
                    {
                        m_broFitRange = novac::WavelengthRange(from, to);
                    }
                }

                if (polyOrder >= 0)
                {
                    m_broPolynomialOrder = polyOrder;
                }
            }
            else if (line.find("<SelectionSettings>") != std::string::npos)
            {
                auto minDifferenceStr = novac::ParseXmlString("MinColumnDifference", line);
                if (!minDifferenceStr.empty())
                {
                    m_ratioEvaluationSettings.minimumInPlumeColumnDifference = std::atof(minDifferenceStr.c_str());
                }

                auto minCompletenessStr = novac::ParseXmlString("MinPlumeCompleteness", line);
                if (!minCompletenessStr.empty())
                {
                    const double completenessLimit = std::atof(minCompletenessStr.c_str());
                    if (completenessLimit > 0.49 && completenessLimit < 1.01)
                    {
                        m_ratioEvaluationSettings.minimumPlumeCompleteness = completenessLimit;
                    }
                }

                auto minScanAngleStr = novac::ParseXmlString("MinScanAngle", line);
                auto maxScanAngleStr = novac::ParseXmlString("MaxScanAngle", line);
                if (!minScanAngleStr.empty() && !maxScanAngleStr.empty())
                {
                    double from = std::atof(minScanAngleStr.c_str());
                    double to = std::atof(maxScanAngleStr.c_str());
                    if (from < to)
                    {
                        m_ratioEvaluationSettings.minimumScanAngle = from;
                        m_ratioEvaluationSettings.maximumScanAngle = to;
                    }
                }

                auto minSaturationRatioStr = novac::ParseXmlString("MinSaturationRatio", line);
                auto maxSaturationRatioStr = novac::ParseXmlString("MaxSaturationRatio", line);
                if (!minSaturationRatioStr.empty() && !maxSaturationRatioStr.empty())
                {
                    double from = std::atof(minSaturationRatioStr.c_str());
                    double to = std::atof(maxSaturationRatioStr.c_str());
                    if (from < to)
                    {
                        m_ratioEvaluationSettings.minSaturationRatio = from;
                        m_ratioEvaluationSettings.maxSaturationRatio = to;
                    }
                }

                int value = novac::ParseXmlInteger("MinInPlumeSpectra", line, 0);
                if (value > 0)
                {
                    m_ratioEvaluationSettings.minNumberOfSpectraInPlume = value;
                }

                value = novac::ParseXmlInteger("MinOutPlumeSpectra", line, 0);
                if (value > 0)
                {
                    m_ratioEvaluationSettings.numberOfSpectraOutsideOfPlume = value;
                }

                auto flanksStr = novac::ParseXmlInteger("RequireTwoFlanks", line);
                m_ratioEvaluationSettings.requireVisiblePlumeEdges = (flanksStr == 1);
            }
            else if (line.find("<DoasFitType>") != std::string::npos)
            {
                auto fitType = novac::ParseXmlInteger("DoasFitType", line, -1);
                switch (fitType)
                {
                case 0:  m_doasFitType = novac::FIT_TYPE::FIT_HP_DIV; break;
                case 1:  m_doasFitType = novac::FIT_TYPE::FIT_HP_SUB; break;
                default:  m_doasFitType = novac::FIT_TYPE::FIT_POLY; break;
                }
            }
            else if (line.find("<Unit>") != std::string::npos)
            {
                auto unit = (novac::CrossSectionUnit)novac::ParseXmlInteger("Unit", line, (int)novac::CrossSectionUnit::cm2_molecule);
                if (unit == novac::CrossSectionUnit::ppmm || unit == novac::CrossSectionUnit::cm2_molecule)
                {
                    m_crossSectionUnit = unit;
                }
            }
        }
    }
    catch (std::exception&)
    {
    }
}

void RatioCalculationController::SaveSetup(const std::string& setupFilePath)
{
    try
    {
        std::ofstream dst(setupFilePath, std::ios::out);
        dst << "<RatioCalculationDlg>" << std::endl;
        for (int rowIdx = 0; rowIdx < static_cast<int>(m_references.size()); ++rowIdx)
        {
            const auto& reference = m_references[rowIdx];

            dst << "\t<Reference>";
            dst << "<Name>" << reference.m_name << "</Name>";
            dst << "<Path>" << reference.m_path << "</Path>";
            dst << "<IncludeInMajor>" << reference.m_includeInMajor << "</IncludeInMajor>";
            dst << "<IncludeInMinor>" << reference.m_includeInMinor << "</IncludeInMinor>";
            dst << "<Calculate>" << reference.m_automaticallyCalculate << "</Calculate>";
            dst << "</Reference>" << std::endl;
        }

        dst << "\t<DoasFitType>" << (int)m_doasFitType << "</DoasFitType>" << std::endl;

        dst << "\t<SO2Setup>";
        dst << "<From>" << m_so2FitRange.low << "</From>";
        dst << "<To>" << m_so2FitRange.high << "</To>";
        dst << "<Poly>" << m_so2PolynomialOrder << "</Poly>";
        dst << "</SO2Setup>" << std::endl;

        dst << "\t<BrOSetup>";
        dst << "<From>" << m_broFitRange.low << "</From>";
        dst << "<To>" << m_broFitRange.high << "</To>";
        dst << "<Poly>" << m_broPolynomialOrder << "</Poly>";
        dst << "</BrOSetup>" << std::endl;

        dst << "\t<SelectionSettings>";
        dst << "<MinColumnDifference>" << m_ratioEvaluationSettings.minimumInPlumeColumnDifference << "</MinColumnDifference>";
        dst << "<MinPlumeCompleteness>" << m_ratioEvaluationSettings.minimumPlumeCompleteness << "</MinPlumeCompleteness>";
        dst << "<MinInPlumeSpectra>" << m_ratioEvaluationSettings.minNumberOfSpectraInPlume << "</MinInPlumeSpectra>";
        dst << "<MinOutPlumeSpectra>" << m_ratioEvaluationSettings.numberOfSpectraOutsideOfPlume << "</MinOutPlumeSpectra>";
        dst << "<RequireTwoFlanks>" << m_ratioEvaluationSettings.requireVisiblePlumeEdges << "</RequireTwoFlanks>";
        dst << "<MinScanAngle>" << m_ratioEvaluationSettings.minimumScanAngle << "</MinScanAngle>";
        dst << "<MaxScanAngle>" << m_ratioEvaluationSettings.maximumScanAngle << "</MaxScanAngle>";
        dst << "<MinSaturationRatio>" << m_ratioEvaluationSettings.minSaturationRatio << "</MinSaturationRatio>";
        dst << "<MaxSaturationRatio>" << m_ratioEvaluationSettings.maxSaturationRatio << "</MaxSaturationRatio>";
        dst << "</SelectionSettings>" << std::endl;

        dst << "<Unit>" << (int)m_crossSectionUnit << "</Unit>";

        dst << "</RatioCalculationDlg>" << std::endl;
    }
    catch (std::exception&)
    {
    }
}

void RatioCalculationController::SetupPakFileList(const std::vector<std::string>& pakFiles)
{
    m_pakfiles = pakFiles;
    m_currentPakFileIdx = 0;
}

std::vector<std::string> RatioCalculationController::ListPakFiles() const
{
    return m_pakfiles;
}

size_t RatioCalculationController::NumberOfPakFilesInSetup() const
{
    return m_pakfiles.size();
}


void SetupFitWindowReferences(novac::CFitWindow& window, const std::vector<ReferenceForRatioCalculation>& references, const novac::WavelengthRange& wavelengthRange, bool isMajor)
{
    window.ringCalculation = novac::RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING;
    window.includeIntensitySpacePolyominal = true;

    for (const ReferenceForRatioCalculation& ref : references)
    {
        const bool referenceShouldBeIncludedInThisWindow = (ref.m_includeInMajor && isMajor) || (ref.m_includeInMinor && !isMajor);

        if (!referenceShouldBeIncludedInThisWindow)
        {
            continue;
        }

        if (ref.m_path != "")
        {
            window.ref[window.nRef].m_path = ref.m_path;
            window.ref[window.nRef].m_specieName = ref.m_name;
            window.ref[window.nRef].m_shiftOption = novac::SHIFT_TYPE::SHIFT_FIX;
            window.ref[window.nRef].m_shiftValue = 0.0;
            window.ref[window.nRef].m_squeezeOption = novac::SHIFT_TYPE::SHIFT_FIX;
            window.ref[window.nRef].m_squeezeValue = 1.0;

            window.nRef++;
        }
        else if (ref.m_automaticallyCalculate && ref.specie == StandardDoasSpecie::RING)
        {
            window.ringCalculation = novac::RING_CALCULATION_OPTION::CALCULATE_RING;
        }
        else if (ref.m_automaticallyCalculate && ref.specie == StandardDoasSpecie::RING_LAMBDA4)
        {
            window.ringCalculation = novac::RING_CALCULATION_OPTION::CALCULATE_RING_X2;
        }
    }

    // Make sure that all the references can be read.
    if (!novac::ReadReferences(window))
    {
        throw std::invalid_argument("failed to read all references");
    }

    // Use the properties of the first (major) reference for the window
    window.name = window.ref[0].m_specieName;

    if (window.ref[0].m_data->m_waveLength.size() == 0)
    {
        std::stringstream message;
        message << "failed to set the fit range,the reference " << window.ref[0].m_specieName << " does not have a wavelength calibration";
        throw std::invalid_argument(message.str());
    }

    // Setup the channel range where the fit should be done.
    const double fractionalFitLow = window.ref[0].m_data->FindWavelength(wavelengthRange.low);
    const double fractionalFitHigh = window.ref[0].m_data->FindWavelength(wavelengthRange.high);
    if (fractionalFitLow < -0.5 || fractionalFitHigh < -0.5)
    {
        std::stringstream message;
        message << "failed to set the fit range,the reference " << window.ref[0].m_specieName;
        message << " does not cover the fit range: " << wavelengthRange.low << " to " << wavelengthRange.high << " nm";
        throw std::invalid_argument(message.str());
    }
    window.fitLow = static_cast<int>(std::round(fractionalFitLow));
    window.fitHigh = static_cast<int>(std::round(fractionalFitHigh));
}

void RatioCalculationController::VerifyReferenceSetup()
{
    const ReferenceForRatioCalculation* so2Reference = GetReferenceFor(*this, StandardDoasSpecie::SO2);
    if (so2Reference == nullptr || so2Reference->m_path.empty() || !so2Reference->m_includeInMajor)
    {
        throw std::invalid_argument("the SO2 reference must be setup and included in the major (SO2) window.");
    }

    const ReferenceForRatioCalculation* broReference = GetReferenceFor(*this, StandardDoasSpecie::BRO);
    if (broReference == nullptr || broReference->m_path.empty() || !broReference->m_includeInMinor)
    {
        throw std::invalid_argument("the BrO reference must be setup and included in the minor (BrO) window.");
    }
}

std::shared_ptr<RatioCalculationFitSetup> RatioCalculationController::SetupFitWindows()
{
    // Start by verifying the setup
    VerifyReferenceSetup();

    auto result = std::make_shared<RatioCalculationFitSetup>();

    result->so2Window.polyOrder = m_so2PolynomialOrder;
    result->so2Window.fitType = m_doasFitType;
    result->broWindow.polyOrder = m_broPolynomialOrder;
    result->broWindow.fitType = m_doasFitType;

    // Setup the references in the fit windows. Notice that the order of the references is different since SO2 and BrO should always be first in their respective window.
    SetupFitWindowReferences(
        result->so2Window,
        m_references,
        m_so2FitRange,
        true);

    SetupFitWindowReferences(
        result->broWindow,
        std::vector<ReferenceForRatioCalculation>{
            m_references[1],
            m_references[0],
            m_references[2],
            m_references[3],
            m_references[4],
        },
        m_broFitRange,
        false);

    return result;
}

novac::BasicScanEvaluationResult RatioCalculationController::DoInitialEvaluation(novac::IScanSpectrumSource& scan, std::shared_ptr<RatioCalculationFitSetup> ratioFitWindows)
{
    novac::BasicScanEvaluationResult result;
    novac::LogContext context; // TODO: Get from input

    // For each spectrum in the scan, do a DOAS evaluation
    novac::CSpectrum measuredSkySpectrum;
    int readSpectrumReturnCode = scan.GetSky(measuredSkySpectrum);
    if (0 != readSpectrumReturnCode || measuredSkySpectrum.m_length == 0)
    {
        throw std::invalid_argument("cannot perform an evaluation on: '" + scan.GetFileName() + "' no sky spectrum found.");
    }

    novac::CSpectrum measuredDarkSpectrum;
    readSpectrumReturnCode = scan.GetDark(measuredDarkSpectrum);
    if (0 != readSpectrumReturnCode || measuredDarkSpectrum.m_length == 0)
    {
        throw std::invalid_argument("cannot perform an evaluation on: '" + scan.GetFileName() + "' no dark spectrum found.");
    }
    measuredSkySpectrum.Sub(measuredDarkSpectrum);

    novac::CFitWindow localCopyOfWindow = ratioFitWindows->so2Window;

    if (localCopyOfWindow.fitType != novac::FIT_TYPE::FIT_HP_DIV)
    {
        const auto filteredSkySpectrum = novac::DoasFitPreparation::PrepareSkySpectrum(measuredSkySpectrum, localCopyOfWindow.fitType);
        (void)AddAsSky(localCopyOfWindow, filteredSkySpectrum, novac::SHIFT_TYPE::SHIFT_FREE);
    }

    novac::DoasFit doas;
    doas.Setup(localCopyOfWindow);

    // TODO: This could be the basis for a (future) scan evaluation class based on the new DoasFit class...
    scan.ResetCounter();
    novac::CSpectrum measuredSpectrum;
    novac::SpectrometerModel spectrometerModel = GetModelForMeasurement(measuredSkySpectrum.m_info.m_device);
    while (0 == scan.GetNextMeasuredSpectrum(context, measuredSpectrum))
    {
        // Check the intensity and save this, such that we can use this later to verify if the evaluated column was good or not.
        measuredSpectrum.m_info.m_peakIntensity = (float)measuredSpectrum.MaxValue(0, measuredSpectrum.m_length - 2);
        measuredSpectrum.m_info.m_fitIntensity = (float)measuredSpectrum.MaxValue(localCopyOfWindow.fitLow, localCopyOfWindow.fitHigh);

        // Dark-correct and prepare the spectrum for the fit
        measuredSpectrum.Sub(measuredDarkSpectrum);
        const auto filteredMeasuredSpectrum = novac::DoasFitPreparation::PrepareMeasuredSpectrum(measuredSpectrum, measuredSkySpectrum, localCopyOfWindow.fitType);

        // do the actual DOAS fit.
        novac::DoasResult doasResult;
        doas.Run(filteredMeasuredSpectrum.data(), filteredMeasuredSpectrum.size(), doasResult);

        // Convert the DoasResult into an CEvaluationResult
        novac::CEvaluationResult evaluationResult = doasResult;

        // Check if the measurement was good or not
        evaluationResult.CheckGoodnessOfFit(measuredSpectrum.m_info, &spectrometerModel);

        result.AppendResult(evaluationResult, measuredSpectrum.m_info);
    }

    return result;
}

bool RatioCalculationController::HasMoreScansToEvaluate() const
{
    if (m_pakfiles.size() == 0 || (size_t)m_currentPakFileIdx >= m_pakfiles.size())
    {
        return false;
    }
    return true;
}

void RatioCalculationController::ResetResults()
{
    m_currentPakFileIdx = 0;
    m_results.clear();
}

RatioCalculationResult RatioCalculationController::EvaluateNextScan(std::shared_ptr<RatioCalculationFitSetup> ratioFitWindows)
{
    if (!HasMoreScansToEvaluate())
    {
        RatioCalculationResult result;
        result.debugInfo.errorMessage = "No more scans to evaluate";
        return result;
    }

    const auto& pakFileName = m_pakfiles[m_currentPakFileIdx++];
    novac::CScanFileHandler scan(m_log);
    novac::LogContext context;
    scan.CheckScanFile(context, pakFileName);

    const auto initialResult = DoInitialEvaluation(scan, ratioFitWindows);

    return EvaluateScan(scan, initialResult, ratioFitWindows);
}

RatioCalculationResult RatioCalculationController::EvaluateScan(
    novac::IScanSpectrumSource& scan,
    const novac::BasicScanEvaluationResult& initialResult,
    std::shared_ptr<RatioCalculationFitSetup> ratioFitWindows)
{
    RatioCalculationResult result;
    result.filename = scan.GetFileName();
    result.initialEvaluation = initialResult;
    result.startTime = scan.GetScanStartTime();
    result.endTime = scan.GetScanStopTime();
    result.deviceSerial = scan.GetDeviceSerial();
    result.evaluatedAt.SetToNow();

    Configuration::CDarkSettings darkSettings; // default dark-settings

    // Calculate the properties of the plume (low large portion of the plume we see and at what angle).
    novac::CPlumeInScanProperty plumeInScanProperties;
    novac::CalculatePlumeOffset(initialResult, 0, plumeInScanProperties);
    const bool plumeIsVisible = novac::CalculatePlumeCompleteness(initialResult, 0, plumeInScanProperties);
    result.plumeInScanProperties = plumeInScanProperties;
    if (!plumeIsVisible)
    {
        result.debugInfo.errorMessage = "No visible plume";
        m_results.push_back(result);
        return result;
    }

    const auto spectrometerModel = GetModelForMeasurement(initialResult.m_specInfo.front().m_device);

    // Setup the ratio evaluation and run it.
    novac::RatioEvaluation ratioEvaluation{ m_ratioEvaluationSettings, darkSettings };
    ratioEvaluation.SetupFitWindows(ratioFitWindows->so2Window, std::vector<novac::CFitWindow>{ ratioFitWindows->broWindow });
    ratioEvaluation.SetupFirstResult(initialResult, plumeInScanProperties, &spectrometerModel);
    const auto ratios = ratioEvaluation.Run(scan, &result.debugInfo);

    // Extract the result
    if (ratios.size() > 0)
    {
        result.ratio = ratios.front();
    }

    // Update the last result as well.
    m_results.push_back(result);

    return result;
}

const novac::SpectrometerModel RatioCalculationController::GetModelForMeasurement(const std::string& deviceSerial) const
{
    if (m_spectrometerModel != nullptr)
    {
        return *(m_spectrometerModel);
    }

    return novac::CSpectrometerDatabase::GetInstance().GuessModelFromSerial(deviceSerial);
}

ReferenceForRatioCalculation* GetReferenceFor(RatioCalculationController& controller, StandardDoasSpecie specie)
{
    for (auto& ref : controller.m_references)
    {
        if (ref.specie == specie)
        {
            return &ref;
        }
    }

    return nullptr; // not found
}

ReferenceForRatioCalculation* GetReferenceWithName(RatioCalculationController& controller, const std::string& name)
{
    for (auto& ref : controller.m_references)
    {
        if (ref.m_name == name)
        {
            return &ref;
        }
    }

    return nullptr; // not found
}

const RatioCalculationResult* GetFirstSuccessfulResult(const std::vector<RatioCalculationResult>& allResults)
{
    for (int ii = 0; ii < static_cast<int>(allResults.size()); ++ii)
    {
        if (allResults[ii].RatioCalculationSuccessful())
        {
            return &allResults[ii];
        }
    }

    return nullptr;
}

void RatioCalculationController::SaveResultsToCsvFile(const std::string& filename, const std::vector< RatioCalculationResult>& resultsToSave, bool overwrite, std::string columnSeparator)
{
    if (resultsToSave.size() == 0)
    {
        return; // too few results to write to file
    }

    const bool writeHeaderLine = !novac::IsExistingFile(filename);

    std::ios::openmode fileMode = overwrite ? std::ios::out : std::ios::app;

    std::ofstream file(filename, fileMode);

    // write the header
    if (writeHeaderLine)
    {
        file << "Device" << columnSeparator;
        file << "ScanStartedAt" << columnSeparator;
        file << "ScanEndedAt" << columnSeparator;
        file << "EvaluatedAt" << columnSeparator;
        file << "Filename" << columnSeparator;
        file << "RatioSuccessfullyCalculated" << columnSeparator;
        file << "BrO/SO2 Ratio" << columnSeparator << "BrO/SO2 RatioError" << columnSeparator;
        file << "BrODetectionSignificant" << columnSeparator;
        file << "PlumeCompleteness" << columnSeparator;
        file << "PlumeCenter" << columnSeparator;
        file << "InPlumeSpectrum_ExposureNum" << columnSeparator;
        file << "OutOfPlumeSpectrum_ExposureNum" << columnSeparator;

        // Write the specie-names from the first successful result, assuming that this hasn't changed for the other results.
        const auto firstSuccessfulResult = GetFirstSuccessfulResult(resultsToSave);
        if (firstSuccessfulResult != nullptr)
        {
            for (size_t windowIdx = 0; windowIdx < firstSuccessfulResult->debugInfo.doasResults.size(); ++windowIdx)
            {
                const size_t humanFriendlyWindowIdx = windowIdx + 1; // humans start counting at 1
                const auto& doasResult = firstSuccessfulResult->debugInfo.doasResults[windowIdx];

                file << "Window" << humanFriendlyWindowIdx << "_FitLow" << columnSeparator;
                file << "Window" << humanFriendlyWindowIdx << "_FitHigh" << columnSeparator;
                file << "Window" << humanFriendlyWindowIdx << "_Chi2" << columnSeparator;
                file << "Window" << humanFriendlyWindowIdx << "_Delta" << columnSeparator;

                for (const auto& referenceResult : doasResult.referenceResult)
                {
                    file << "Window" << humanFriendlyWindowIdx << "_" << referenceResult.name << "_Column" << columnSeparator;
                    file << "Window" << humanFriendlyWindowIdx << "_" << referenceResult.name << "_ColumnError" << columnSeparator;
                    file << "Window" << humanFriendlyWindowIdx << "_" << referenceResult.name << "_Shift" << columnSeparator;
                    file << "Window" << humanFriendlyWindowIdx << "_" << referenceResult.name << "_ShiftError" << columnSeparator;
                    file << "Window" << humanFriendlyWindowIdx << "_" << referenceResult.name << "_Squeeze" << columnSeparator;
                    file << "Window" << humanFriendlyWindowIdx << "_" << referenceResult.name << "_SqueezeError" << columnSeparator;
                }
            }
        }

        file << std::endl;
    }

    // write the data
    for (const auto& result : resultsToSave)
    {
        file << result.deviceSerial << columnSeparator;
        file << result.startTime << columnSeparator;
        file << result.endTime << columnSeparator;
        file << result.evaluatedAt << columnSeparator;
        file << result.filename << columnSeparator;
        file << result.RatioCalculationSuccessful() << columnSeparator; // prints '1' whenever the error message is empty, i.e. the ratio is calculated without any errors.
        file << result.ratio.ratio << columnSeparator << result.ratio.error << columnSeparator;
        file << result.SignificantMinorSpecieDetection() << columnSeparator;
        file << result.plumeInScanProperties.completeness << columnSeparator;
        file << result.plumeInScanProperties.plumeCenter << columnSeparator;
        file << result.debugInfo.inPlumeSpectrum.m_info.m_numSpec << columnSeparator;
        file << result.debugInfo.outOfPlumeSpectrum.m_info.m_numSpec << columnSeparator;

        for (const auto& doasResult : result.debugInfo.doasResults)
        {
            file << doasResult.fitLow << columnSeparator;
            file << doasResult.fitHigh << columnSeparator;
            file << doasResult.chiSquare << columnSeparator;
            file << doasResult.delta << columnSeparator;

            for (const auto& referenceResult : doasResult.referenceResult)
            {
                file << referenceResult.column << columnSeparator << referenceResult.columnError << columnSeparator;
                file << referenceResult.shift << columnSeparator << referenceResult.shiftError << columnSeparator;
                file << referenceResult.squeeze << columnSeparator << referenceResult.squeezeError << columnSeparator;
            }
        }

        file << std::endl;
    }
}

void RatioCalculationController::SaveSpectraToPakFile(const std::string& outputDirectory, const RatioCalculationResult& resultToSave)
{
    // Save the spectra as 1) sky, 2) dark and 3) inplume. This simulates the usual order of the spectra in novac data.
    // (Compare with PlumeSpectrumSelector::CreatePlumeSpectrumFile)

    std::stringstream spectrumOutputFileName;
    spectrumOutputFileName << outputDirectory << "/PlumeSpectra_" << resultToSave.deviceSerial;
    spectrumOutputFileName << "_" << FormatDate(resultToSave.skySpectrumInfo);
    spectrumOutputFileName << "_" << FormatTimestamp(resultToSave.startTime);
    spectrumOutputFileName << "_0.pak";

    const std::string filename = spectrumOutputFileName.str();

    // Since the measurement is already dark-corrected, create an all-zero dark-spectrum to use
    novac::CSpectrum darkSpectrum;
    darkSpectrum.m_length = resultToSave.debugInfo.inPlumeSpectrum.m_length;

    // Write the data to file.
    novac::CSpectrumIO spectrumWriter;
    spectrumWriter.AddSpectrumToFile(filename, resultToSave.debugInfo.outOfPlumeSpectrum, nullptr, 0, true);
    spectrumWriter.AddSpectrumToFile(filename, darkSpectrum);
    spectrumWriter.AddSpectrumToFile(filename, resultToSave.debugInfo.inPlumeSpectrum);
}

void RatioCalculationController::SaveSpectraToStdFile(const std::string& outputDirectory, const RatioCalculationResult& resultToSave)
{
    std::stringstream outOfPlumeSpectrumFileName;
    outOfPlumeSpectrumFileName << outputDirectory << "/ReferenceSpectrum_" << resultToSave.deviceSerial;
    outOfPlumeSpectrumFileName << "_" << FormatDate(resultToSave.skySpectrumInfo);
    outOfPlumeSpectrumFileName << "_" << FormatTimestamp(resultToSave.startTime);
    outOfPlumeSpectrumFileName << "_0.std";

    std::stringstream inPlumeSpectrumFileName;
    inPlumeSpectrumFileName << outputDirectory << "/InPlumeSpectrum_" << resultToSave.deviceSerial;
    inPlumeSpectrumFileName << "_" << FormatDate(resultToSave.skySpectrumInfo);
    inPlumeSpectrumFileName << "_" << FormatTimestamp(resultToSave.startTime);
    inPlumeSpectrumFileName << "_0.std";

    novac::CSTDFile fileWriter;
    fileWriter.WriteSpectrum(resultToSave.debugInfo.outOfPlumeSpectrum, outOfPlumeSpectrumFileName.str(), 1);
    fileWriter.WriteSpectrum(resultToSave.debugInfo.inPlumeSpectrum, inPlumeSpectrumFileName.str(), 1);
}