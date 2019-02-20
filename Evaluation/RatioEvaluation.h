#pragma once

#include "ScanEvaluationBase.h"
#include "FitWindow.h"
#include <vector>

namespace Evaluation
{
    class BasicScanEvaluationResult;

    struct Ratio
    {
        double value; //<- The estimated quotient.
        double error; //<- An estimation of the error in the quotient.
    };

    // TODO: Move to Configuration
    struct RatioEvaluationSettings
    {
        // The minimum number of spectra which needs to be selected in the plume for the ratio calculation to be successful.
        int minNumberOfSpectraInPlume = 10;

        // The minimum (SO2) column for the selected spectra in the plume.
        double minInPlumeColumn = 1e19;

        // The minimum number of spectra which needs to be averaged outside of the plume for the calculation to be successful.
        int minNumberOfReferenceSpectra = 10;
    };

    /** The class RatioEvaluation helps with evaluating the ratios of specific elements (e.g. BrO/SO2-ratio)
        from single scans by evaluating the species in separate fit-windows and calculating a ratio. 
        This needs a number of things:
            1) A series of fit-windows from which to calculate a column amount of a specific specie.
            2) A process for selecting which spectrum to use as sky and which spectrum to use as the spectrum to evaluate.
                typically is the spectrum to evaluate selected from the highest few spectra in the plume and the sky
                selected from a few spectra outside of the plume.
                This in turn needs an already evaluated scan from which it is possible to determine which spectrum
                to use as sky and which to use as spectrum to evaluate.
            3) A way to configure all this from the user!
        */
    class RatioEvaluation : public ScanEvaluationBase
    {
    public:
        RatioEvaluation(const RatioEvaluationSettings& settings);

        void Initialize(const CFitWindow& mainFitWindow, const std::vector<CFitWindow>& referenceFitWindows)
        {
            m_masterFitWindow = mainFitWindow;
            m_referenceFit = referenceFitWindows;
        }

        /** Runs all the ratio evaluations.
            @return A vector with all the calculated quotients. The length of this vector
                equals the number of reference-fit windows passed to 'Initialize' which must have been called before this.
            @return an empty vector if the evaluations fail. */
        std::vector<Ratio> Run();

    private:
        /** The fit window against which the ratio will be calculated (typically SO2). 
            The ratio will be performed against the first specie in the fit window. */
        CFitWindow m_masterFitWindow;

        /** The fit windows which will be used to estimate the minor specie, typically BrO.
            The name of each fit-window must equal the name of the sub-specie to calculate.
            Being a vector makes it possible to evaluate multiple species in the same round. */
        std::vector<CFitWindow> m_referenceFit;

        /** The settings for how the ratio-calculation should be performed. */
        RatioEvaluationSettings m_settings;
    };

    /** @return true if the provided evaluation result is suitable for performing a ratio-evaluation. */
    bool IsSuitableScanForRatioEvaluation(const RatioEvaluationSettings& settings, const BasicScanEvaluationResult& result);

    /** Estimates which spectra should be used for a ratio-evaluation, assuming that one should be performed. 
        @param referenceSpectra Will be filled with the index of the spectra which should be averaged to a reference spectrum.
        @param inPlumeSpectra Will be filled with the index of the spectra which should be averaged to an in-plume spectrum. */
    void SelectSpectraForRatioEvaluation(const RatioEvaluationSettings& settings, const BasicScanEvaluationResult& scanResult, std::vector<int>& referenceSpectra, std::vector<int>& inPlumeSpectra);

}