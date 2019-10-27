#pragma once

class CPlumeInScanProperty;

namespace Evaluation
{
    class BasicScanEvaluationResult;

    struct PlumeSpectrumSelectionSettings
    {
        int minNumberOfSpectraInPlume = 5;
        int minNumberOfReferenceSpectra = 5;
    };


    // This is a helper class for selecting 'in plume' and 'out of plume' 
    //  (based on some basic selection criteria) spectra with 
    //  the intention of creating reference spectra for performing e.g. 
    //  ratio evaluations or detailed spectral analysis of compounds in the plume.
    class PlumeSpectrumSelector
    {
        static bool IsSuitableScanForRatioEvaluation(const PlumeSpectrumSelectionSettings& settings, const BasicScanEvaluationResult& scanResult, const CPlumeInScanProperty& properties);


    };
}
