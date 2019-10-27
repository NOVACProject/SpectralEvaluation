#include <numeric>

#include <SpectralEvaluation/Evaluation/PlumeSpectrumSelector.h>

#include <SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h>
#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>

using namespace Evaluation;

void PlumeSpectrumSelector::SelectSpectra(
    const BasicScanEvaluationResult& scanResult,
    const CPlumeInScanProperty& properties,
    int mainSpecieIndex,
    std::vector<int>& referenceSpectra,
    std::vector<int>& inPlumeSpectra)
{
    PlumeSpectrumSelectionSettings defaultSettings;
    return SelectSpectra(defaultSettings, scanResult, properties, mainSpecieIndex, referenceSpectra, inPlumeSpectra);
}


void PlumeSpectrumSelector::SelectSpectra(
    const PlumeSpectrumSelectionSettings& settings,
    const BasicScanEvaluationResult& scanResult,
    const CPlumeInScanProperty& properties,
    int mainSpecieIndex,
    std::vector<int>& referenceSpectra,
    std::vector<int>& inPlumeSpectra)
{
    referenceSpectra.clear();
    inPlumeSpectra.clear();

    const size_t referenceRegionWidth = 10U;

    if (scanResult.m_spec.size() <= referenceRegionWidth)
    {
        return; // Cannot retrieve a region, too few spectra...
    }

    // Step 1, find the in-plume region.
    for (size_t idx = 0; idx < scanResult.m_spec.size(); ++idx)
    {
        if (scanResult.m_specInfo[idx].m_scanAngle >= properties.plumeEdgeLow &&
            scanResult.m_specInfo[idx].m_scanAngle <= properties.plumeEdgeHigh &&
            scanResult.m_spec[idx].m_referenceResult[mainSpecieIndex].m_column >= settings.minInPlumeColumn)
        {
            inPlumeSpectra.push_back((int)idx);
        }
    }

    // Step 2, find the reference region (10 adjacent spectra with lowest avg column value)
    size_t startIdxOfReferenceRegion = 0U;
    double lowestMeanColumnValue = std::numeric_limits<double>::max();
    for (size_t startIdxCandidate = 0U; startIdxCandidate < scanResult.m_spec.size() - referenceRegionWidth; ++startIdxCandidate)
    {
        const double meanColumnValue = AverageColumnValue(scanResult, mainSpecieIndex, startIdxCandidate, startIdxCandidate + referenceRegionWidth);
        if (meanColumnValue < lowestMeanColumnValue)
        {
            lowestMeanColumnValue = meanColumnValue;
            startIdxOfReferenceRegion = startIdxCandidate;
        }
    }

    referenceSpectra = std::vector<int>(referenceRegionWidth);
    std::iota(begin(referenceSpectra), end(referenceSpectra), (int)startIdxOfReferenceRegion);

    return;
}

double PlumeSpectrumSelector::AverageColumnValue(
    const BasicScanEvaluationResult& scanResult,
    int specieIndex,
    size_t startIdx,
    size_t endIdx)
{
    double sum = 0.0;
    for (size_t ii = startIdx; ii < endIdx; ++ii)
    {
        sum += scanResult.m_spec[ii].m_referenceResult[specieIndex].m_column;
    }
    return sum / (double)(endIdx - startIdx);
}

bool PlumeSpectrumSelector::IsSuitableScanForRatioEvaluation(
    const PlumeSpectrumSelectionSettings& settings,
    const BasicScanEvaluationResult& scanResult,
    const CPlumeInScanProperty& properties)
{
    if (scanResult.m_spec.size() < settings.minNumberOfSpectraInPlume + settings.minNumberOfReferenceSpectra)
    {
        return false; // not enough spectra
    }
    if (properties.completeness < 0.7)
    {
        return false;
    }
    return true; // TODO: Add more checks...
}