#include <numeric>

#include <SpectralEvaluation/Evaluation/PlumeSpectrumSelector.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h>
#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>

using namespace Evaluation;

void PlumeSpectrumSelector::SelectSpectra(
    FileHandler::CScanFileHandler& scanFile,
    const BasicScanEvaluationResult& scanResult,
    const CPlumeInScanProperty& properties,
    int mainSpecieIndex,
    std::vector<size_t>& referenceSpectra,
    std::vector<size_t>& inPlumeSpectra)
{
    PlumeSpectrumSelectionSettings defaultSettings;
    return SelectSpectra(defaultSettings, scanFile, scanResult, properties, mainSpecieIndex, referenceSpectra, inPlumeSpectra);
}


void PlumeSpectrumSelector::SelectSpectra(
    const PlumeSpectrumSelectionSettings& settings,
    FileHandler::CScanFileHandler& scanFile,
    const BasicScanEvaluationResult& scanResult,
    const CPlumeInScanProperty& properties,
    int mainSpecieIndex,
    std::vector<size_t>& referenceSpectra,
    std::vector<size_t>& inPlumeSpectra)
{
    referenceSpectra.clear();
    inPlumeSpectra.clear();

    if (scanResult.m_spec.size() <= settings.minNumberOfSpectraOutsideOfPlume)
    {
        return; // Cannot retrieve a region, too few spectra...
    }

    // Find a proprosal for the in-plume region.
    auto inPlumeProposal = FindSpectraInPlume(scanResult, properties, mainSpecieIndex, settings);

    inPlumeProposal = FilterSpectraUsingIntensity(inPlumeProposal, scanFile, settings);

    // Step 2, find the reference region (10 adjacent spectra with lowest avg column value)
    const size_t referenceRegionWidth = 10U;
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

    referenceSpectra = std::vector<size_t>(referenceRegionWidth);
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
    if (scanResult.m_spec.size() < settings.minNumberOfSpectraInPlume + settings.minNumberOfSpectraOutsideOfPlume)
    {
        return false; // not enough spectra
    }
    if (properties.completeness < 0.7)
    {
        return false; // need to see the plume
    }
    return true; // TODO: Add more checks...
}

std::vector<size_t> PlumeSpectrumSelector::FindSpectraInPlume(
    const BasicScanEvaluationResult& scanResult,
    const CPlumeInScanProperty& properties,
    int mainSpecieIndex,
    const PlumeSpectrumSelectionSettings& settings)
{
    std::vector<size_t> indices;

    for (size_t idx = 0; idx < scanResult.m_spec.size(); ++idx)
    {
        if (scanResult.m_specInfo[idx].m_scanAngle >= properties.plumeEdgeLow &&
            scanResult.m_specInfo[idx].m_scanAngle <= properties.plumeEdgeHigh &&
            scanResult.m_spec[idx].m_referenceResult[mainSpecieIndex].m_column >= settings.minInPlumeColumn)
        {
            indices.push_back(idx);
        }
    }

    return indices;
}

std::vector<size_t> PlumeSpectrumSelector::FilterSpectraUsingIntensity(
    const std::vector<size_t>& proposedIndices,
    FileHandler::CScanFileHandler& scanFile,
    const PlumeSpectrumSelectionSettings& settings)
{
    std::vector<size_t> result;
    result.reserve(proposedIndices.size());

    CSpectrum spectrum;
    for (size_t idx : proposedIndices)
    {
        if (scanFile.GetSpectrum(spectrum, (long)idx))
        {
            if (spectrum.MaxValue(0, spectrum.m_length) <= settings.maxSaturationRatio &&
                spectrum.MaxValue(0, spectrum.m_length) >= settings.minSaturationRatio)
            {
                result.push_back(idx);
            }
        }
    }

    return result;
}
