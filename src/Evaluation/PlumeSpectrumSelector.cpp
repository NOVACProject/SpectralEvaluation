#include <SpectralEvaluation/Evaluation/PlumeSpectrumSelector.h>
#include <SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h>
#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>

using namespace Evaluation;

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