#include <SpectralEvaluation/DialogControllers/RatioCalculationController.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/Evaluation/RatioEvaluation.h>
#include <SpectralEvaluation/Configuration/DarkSettings.h>


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

void SetupFitWindowReferences(novac::CFitWindow& window, const std::vector<ReferenceForRatioCalculation>& references, std::string name, bool isMajor)
{
    window.fitType = novac::FIT_TYPE::FIT_POLY;
    window.name = name;
    window.polyOrder = 3;
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
            window.ref[window.nRef].m_squeezeValue = 0.0;

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

    // Setup the references in the fit windows. Notice that the order of the references is different since SO2 and BrO should always be first in their respective window.
    SetupFitWindowReferences(
        result->m_so2Window,
        m_references,
        "SO2",
        true);
    SetupFitWindowReferences(
        result->m_broWindow,
        std::vector<ReferenceForRatioCalculation>{
            m_references[1],
            m_references[0],
            m_references[2],
            m_references[3],
            m_references[4],
        },
        "BrO",
        false);

    // TODO: Setup the fit ranges to something reasonable...

    return result;
}

void RatioCalculationController::EvaluateNextScan()
{
    const auto pakFileName = m_pakfiles[0];
    novac::CScanFileHandler scan;
    scan.CheckScanFile(pakFileName);

    // Do a first evaluation such that we know if there are spectra to evaluate...

    // Now run the ratio evaluation
    // EvaluateScan(scan);
}

bool RatioCalculationController::EvaluateScan(novac::IScanSpectrumSource& scan, const novac::BasicScanEvaluationResult& initialResult, std::shared_ptr<RatioCalculationFitSetup> ratioFitWindows)
{
    std::string errorMessage;

    Configuration::CDarkSettings darkSettings; // default dark-settings

    // Calculate the properties of the plume (low large portion of the plume we see and at what angle).
    novac::CPlumeInScanProperty plumeInScanProperties;
    novac::CalculatePlumeOffset(initialResult, 0, plumeInScanProperties);
    const bool plumeIsVisible = novac::CalculatePlumeCompleteness(initialResult, 0, plumeInScanProperties);
    if (!plumeIsVisible)
    {
        return false;
    }


    // Setup the ratio evaluation and run it.
    novac::RatioEvaluation ratioEvaluation{ m_ratioEvaluationSettings, darkSettings };
    ratioEvaluation.SetupFitWindows(ratioFitWindows->m_so2Window, std::vector<novac::CFitWindow>{ ratioFitWindows->m_broWindow });
    ratioEvaluation.SetupFirstResult(initialResult, plumeInScanProperties);
    const auto ratios = ratioEvaluation.Run(scan, &errorMessage);

    return true;
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