#include <SpectralEvaluation/DialogControllers/RatioCalculationController.h>

RatioCalculationController::RatioCalculationController()
{
    InitializeToDefault();
}

void RatioCalculationController::InitializeToDefault()
{
    m_pakfiles.clear();

    m_so2Window.Clear();
    m_broWindow.Clear();

    // Insert the default species
    m_references.clear();
    m_references.push_back(ReferenceForRatioCalculation(StandardDoasSpecie::SO2, "SO2", "", true, true, false));
    m_references.push_back(ReferenceForRatioCalculation(StandardDoasSpecie::BRO, "BrO", "", false, true, false));
    m_references.push_back(ReferenceForRatioCalculation(StandardDoasSpecie::O3, "O3", "", true, true, false));
    m_references.push_back(ReferenceForRatioCalculation(StandardDoasSpecie::RING, "Ring", "", true, true, true));
    m_references.push_back(ReferenceForRatioCalculation(StandardDoasSpecie::RING_LAMBDA4, "Ringxlambda^4", "", true, true, true));


    // TODO: Setup the fit ranges to something reasonable...
    m_so2Window.fitType = novac::FIT_TYPE::FIT_POLY;
    m_so2Window.name = "SO2";
    m_so2Window.polyOrder = 3;
    m_so2Window.ringCalculation = novac::RING_CALCULATION_OPTION::CALCULATE_RING_X2;
    m_so2Window.includeIntensitySpacePolyominal = true;


    m_broWindow.fitType = novac::FIT_TYPE::FIT_POLY;
    m_broWindow.name = "BrO";
    m_broWindow.polyOrder = 3;
    m_broWindow.ringCalculation = novac::RING_CALCULATION_OPTION::CALCULATE_RING_X2;
    m_broWindow.includeIntensitySpacePolyominal = true;
}

