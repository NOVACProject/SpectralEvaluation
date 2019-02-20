#include "RatioEvaluation.h"
#include "BasicScanEvaluationResult.h"
#include "EvaluationResult.h"

namespace Evaluation
{
    RatioEvaluation::RatioEvaluation(const RatioEvaluationSettings& settings)
        : m_settings(settings)
    {
    }

    std::vector<Ratio> RatioEvaluation::Run()
    {
        std::vector<Ratio> result;
        // TODO: Implement

        return result;
    }

    bool IsSuitableScanForRatioEvaluation(const RatioEvaluationSettings& settings, const BasicScanEvaluationResult& scanResult)
    {
        if (scanResult.m_spec.size() < settings.minNumberOfSpectraInPlume + settings.minNumberOfReferenceSpectra)
        {
            return false; // not enough spectra
        }

        return true; // TODO: Implement
    }

    void SelectSpectraForRatioEvaluation(const RatioEvaluationSettings& /*settings*/, const BasicScanEvaluationResult& /*scanResult*/, std::vector<int>& /*referenceSpectra*/, std::vector<int>& /*inPlumeSpectra*/)
    {
        // TODO: Implement-me
        return;
    }

}