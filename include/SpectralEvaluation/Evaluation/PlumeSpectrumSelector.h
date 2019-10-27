#pragma once
#include <vector>

class CPlumeInScanProperty;

namespace Evaluation
{
    class BasicScanEvaluationResult;

    // This is a helper class for selecting 'in plume' and 'out of plume' 
    //  (based on some basic selection criteria) spectra with 
    //  the intention of creating reference spectra for performing e.g. 
    //  ratio evaluations or detailed spectral analysis of compounds in the plume.
    class PlumeSpectrumSelector
    {
        struct PlumeSpectrumSelectionSettings
        {
            // The minimum number of spectra which needs to be selected in the plume for the ratio calculation to be successful.
            int minNumberOfSpectraInPlume = 7;

            // The minimum (SO2) column for the selected spectra in the plume.
            double minInPlumeColumn = 1e17;

            // The minimum number of spectra which needs to be averaged outside of the plume for the calculation to be successful.
            int minNumberOfReferenceSpectra = 7;
        };

        // Checks the provided evaluated scan using default settings 
        //  and returns the indices of the spectra which can 
        //  be saved as in-plume and out-of-plume spectra.
        // If this check fails, then the two vectors are empty.
        // TODO: This also needs a reference to the ScanFile such that it can check the intensity of the spectra in various regions!
        void SelectSpectra(const BasicScanEvaluationResult& scanResult, const CPlumeInScanProperty& properties, int mainSpecieIndex, std::vector<int>& referenceSpectra, std::vector<int>& inPlumeSpectra);

    private:

        static bool IsSuitableScanForRatioEvaluation(const PlumeSpectrumSelectionSettings& settings, const BasicScanEvaluationResult& scanResult, const CPlumeInScanProperty& properties);

        static void SelectSpectra(const PlumeSpectrumSelectionSettings& settings, const BasicScanEvaluationResult& scanResult, const CPlumeInScanProperty& properties, int mainSpecieIndex, std::vector<int>& referenceSpectra, std::vector<int>& inPlumeSpectra);

        // Calculates the average column value of the given specie in the index [startIdx, endIdx[
        static double AverageColumnValue(const BasicScanEvaluationResult& scanResult, int specieIndex, size_t startIdx, size_t endIdx);

    };
}
