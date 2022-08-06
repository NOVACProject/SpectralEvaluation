#include <SpectralEvaluation/DialogControllers/RatioCalculationController.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/Evaluation/RatioEvaluation.h>
#include <SpectralEvaluation/Configuration/DarkSettings.h>

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
    window.fitType = novac::FIT_TYPE::FIT_POLY;
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
    result->broWindow.polyOrder = m_broPolynomialOrder;

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


RatioCalculationResult RatioCalculationController::EvaluateNextScan(std::shared_ptr<RatioCalculationFitSetup> ratioFitWindows)
{
    if (m_pakfiles.size() == 0 || (size_t)m_currentPakFileIdx >= m_pakfiles.size())
    {
        RatioCalculationResult result;
        return result;
    }

    const auto pakFileName = m_pakfiles[m_currentPakFileIdx];
    novac::CScanFileHandler scan;
    scan.CheckScanFile(pakFileName);

    // TODO: Do a first evaluation, such that we have an initial result.
    novac::BasicScanEvaluationResult initialResult;

    // Now run the ratio evaluation
    auto result = EvaluateScan(scan, initialResult, ratioFitWindows);

    return result;
}

RatioCalculationResult RatioCalculationController::EvaluateScan(
    novac::IScanSpectrumSource& scan,
    const novac::BasicScanEvaluationResult& initialResult,
    std::shared_ptr<RatioCalculationFitSetup> ratioFitWindows)
{
    RatioCalculationResult result;
    result.filename = scan.GetFileName();

    Configuration::CDarkSettings darkSettings; // default dark-settings

    // Calculate the properties of the plume (low large portion of the plume we see and at what angle).
    novac::CPlumeInScanProperty plumeInScanProperties;
    novac::CalculatePlumeOffset(initialResult, 0, plumeInScanProperties);
    const bool plumeIsVisible = novac::CalculatePlumeCompleteness(initialResult, 0, plumeInScanProperties);
    if (!plumeIsVisible)
    {
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

    return result;
}

ReferenceForRatioCalculation* GetReferenceFor(RatioCalculationController& controller, StandardDoasSpecie specie) {
    for (auto& ref : controller.m_references)
    {
        if (ref.specie == specie)
        {
            return &ref;
        }
    }

    return nullptr; // not found
}