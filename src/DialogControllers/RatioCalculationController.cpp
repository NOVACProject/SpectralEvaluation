#include <SpectralEvaluation/DialogControllers/RatioCalculationController.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/Evaluation/DoasFitPreparation.h>
#include <SpectralEvaluation/Evaluation/RatioEvaluation.h>
#include <SpectralEvaluation/Configuration/DarkSettings.h>
#include <SpectralEvaluation/File/XmlUtil.h>

#include <fstream>
#include <sstream>

RatioCalculationController::RatioCalculationController()
    : m_so2FitRange(314.8, 326.8), m_broFitRange(330.6, 352.8)
{
    InitializeToDefault();
}

void RatioCalculationController::InitializeToDefault()
{
    m_pakfiles.clear();

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

    // For each spectrum in the scan, do a DOAS evaluation
    novac::CSpectrum measuredSkySpectrum;
    if (0 != scan.GetSky(measuredSkySpectrum))
    {
        throw std::invalid_argument("cannot perform an evaluation on: '" + scan.GetFileName() + "' no sky spectrum found.");
    }

    novac::CSpectrum measuredDarkSpectrum;
    if (0 != scan.GetDark(measuredDarkSpectrum))
    {
        throw std::invalid_argument("cannot perform an evaluation on: '" + scan.GetFileName() + "' no dark spectrum found.");
    }
    measuredSkySpectrum.Sub(measuredDarkSpectrum);

    const auto filteredSkySpectrum = novac::DoasFitPreparation::PrepareSkySpectrum(measuredSkySpectrum, novac::FIT_TYPE::FIT_POLY);

    novac::CFitWindow localCopyOfWindow = ratioFitWindows->so2Window;
    (void)AddAsSky(localCopyOfWindow, filteredSkySpectrum, novac::SHIFT_TYPE::SHIFT_FREE);

    novac::DoasFit doas;
    doas.Setup(localCopyOfWindow);

    // TODO: This could be the basis for a (future) scan evaluation class based on the new DoasFit class...
    scan.ResetCounter();
    novac::CSpectrum measuredSpectrum;
    while (0 == scan.GetNextMeasuredSpectrum(measuredSpectrum))
    {
        measuredSpectrum.Sub(measuredDarkSpectrum);
        const auto filteredMeasuredSpectrum = novac::DoasFitPreparation::PrepareMeasuredSpectrum(measuredSpectrum, measuredSkySpectrum, novac::FIT_TYPE::FIT_POLY);

        novac::DoasResult doasResult;
        doas.Run(filteredMeasuredSpectrum.data(), filteredMeasuredSpectrum.size(), doasResult);

        novac::CEvaluationResult evaluationResult = doasResult;
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

RatioCalculationResult RatioCalculationController::EvaluateNextScan(std::shared_ptr<RatioCalculationFitSetup> ratioFitWindows)
{
    if (!HasMoreScansToEvaluate())
    {
        RatioCalculationResult result;
        result.errorMessage = "No more scans to evaluate";
        return result;
    }

    const auto& pakFileName = m_pakfiles[m_currentPakFileIdx++];
    novac::CScanFileHandler scan;
    scan.CheckScanFile(pakFileName);

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
        result.errorMessage = "No visible plume";
        m_results.push_back(result);
        return result;
    }

    // Setup the ratio evaluation and run it.
    novac::RatioEvaluation ratioEvaluation{ m_ratioEvaluationSettings, darkSettings };
    ratioEvaluation.SetupFitWindows(ratioFitWindows->so2Window, std::vector<novac::CFitWindow>{ ratioFitWindows->broWindow });
    ratioEvaluation.SetupFirstResult(initialResult, plumeInScanProperties);
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


void RatioCalculationController::SaveResultsToCsvFile(const std::string& filename, std::string columnSeparator) const
{
    if (m_results.size() == 0)
    {
        return; // too few results to write to file
    }

    std::ofstream file(filename);

    // write the header
    {
        file << "Device" << columnSeparator;
        file << "ScanStartedAt" << columnSeparator;
        file << "ScanEndedAt" << columnSeparator;
        file << "EvaluatedAt" << columnSeparator;
        file << "Filename" << columnSeparator;
        file << "Ratio" << columnSeparator << "RatioError" << columnSeparator;
        file << "PlumeCompleteness" << columnSeparator;
        file << "PlumeCenter" << columnSeparator;
        file << "InPlumeSpectrum_ExposureNum" << columnSeparator;
        file << "OutOfPlumeSpectrum_ExposureNum" << columnSeparator;

        // Write the specie-names from the first result, assuming that this hasn't changed for the other results.
        for (size_t windowIdx = 0; windowIdx < m_results.front().debugInfo.doasResults.size(); ++windowIdx)
        {
            const size_t humanFriendlyWindowIdx = windowIdx + 1; // humans start counting at 1
            const auto& doasResult = m_results.front().debugInfo.doasResults[windowIdx];

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

        file << std::endl;
    }

    // write the data
    for (const auto& result : m_results)
    {
        file << result.deviceSerial << columnSeparator;
        file << result.startTime << columnSeparator;
        file << result.endTime << columnSeparator;
        file << result.evaluatedAt << columnSeparator;
        file << result.filename << columnSeparator;
        file << result.ratio.ratio << columnSeparator << result.ratio.error << columnSeparator;
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