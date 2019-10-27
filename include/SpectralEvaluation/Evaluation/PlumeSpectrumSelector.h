#pragma once

namespace Evaluation
{
    class BasicScanEvaluationResult;
    class CPlumeInScanProperty;

    // This is a helper class for selecting 'in plume' and 'out of plume' 
    //  (based on some basic selection criteria) spectra with 
    //  the intention of creating reference spectra for performing e.g. 
    //  ratio evaluations or detailed spectral analysis of compounds in the plume.
    public class PlumeSpectrumSelector
    {
        static bool IsSuitableScanForRatioEvaluation(const RatioEvaluationSettings& settings, const BasicScanEvaluationResult& scanResult, const CPlumeInScanProperty& properties);


    };
}
