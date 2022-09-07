#pragma once

#include <SpectralEvaluation/Evaluation/Ratio.h>
#include <SpectralEvaluation/Evaluation/FitWindow.h>
#include <SpectralEvaluation/Evaluation/DoasFit.h>
#include <SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Configuration/RatioEvaluationSettings.h>
#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>
#include <vector>

namespace Configuration
{
    struct CDarkSettings;
}

namespace novac
{
    class IScanSpectrumSource;

    // A basic structure used to be able to extract information from the RatioEvaluation below.
    // This can be used to e.g. present intermediates to the user for a better understanding of the process.
    struct RatioEvaluationDebugInformation
    {
        // A possible error message, empty if no error has occurred.
        std::string errorMessage;

        // The model name of the spectrometer which generated this measurement.
        std::string spectrometerModel;

        // The full dynamic range (of one readout) of the generating spectrometer.
        double spectrometerFullDynamicRange = 4095.0;

        // The zero-based indices of the spectra used to create the in-plume spectrum
        std::vector<int> plumeSpectrumIndices;

        // The zero-based indices of the spectra used to create the out-of-plume spectrum
        std::vector<int> outOfPlumeSpectrumIndices;

        // Listing the indices (in the original scan) which for some reason were not considered for participating in the in-plume or out-of-plume spectra,
        // together with the reason for not including them.
        std::vector<std::pair<int, std::string>> rejectedSpectrumIndices;

        // The final (dark corrected) in-plume spectrum.
        CSpectrum inPlumeSpectrum;

        // The final (dark corrected) out-of-plume spectrum.
        CSpectrum outOfPlumeSpectrum;

        // The results of each of the evaluations performed (one result for the major evaluation and one for the minor).
        std::vector<DoasResult> doasResults;
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
    class RatioEvaluation
    {
    public:
        RatioEvaluation(const Configuration::RatioEvaluationSettings& settings, const Configuration::CDarkSettings& darkSettings);

        // Calculates the ratio of an evaluation in a major window and a minor window.
        // The ratio is calculated using the first reference result in each DoasResult.
        static Ratio CalculateRatio(const DoasResult& majorWindow, const DoasResult& minorWindow);

        /** Sets up this ratio-evaluation with the fit windows to use for evaluation.
            @param fitWindow Is the fit window for the main specie to evaluate (typically SO2)
            @param referenceFitWindows Are the fit widnows for the minor species to evaluate (typically BrO). */
        void SetupFitWindows(const CFitWindow& fitWindow, const std::vector<CFitWindow>& referenceFitWindows)
        {
            m_masterFitWindow = fitWindow;
            m_referenceFit = referenceFitWindows;
        }

        /** Sets up this ratio-evaluation with the result from the first evaluation which has been done using the main fit window.
            @param result Is the result from evaluation using this FitWindow
            @param properties Is the plume-shape properties for this result. */
        void SetupFirstResult(const BasicScanEvaluationResult& result, const CPlumeInScanProperty& properties)
        {
            m_masterResult = result;
            m_masterResultProperties = properties;
        }

        /** Runs all the ratio evaluations.
            @param scan A handle to the .pak file to evaluate
            @param errorMessage If not null then this will be filled with the reason the evaluation failed in case of errors.
            @return A vector with all the calculated quotients. The length of this vector equals
            the number of reference-fit windows passed to 'SetupFitWindows', which must have been called before this,
            or an empty vector if the evaluations fail.  */
        std::vector<Ratio> Run(IScanSpectrumSource& scan, RatioEvaluationDebugInformation* debugInfo = nullptr);

    private:
        /** The fit window against which the ratio will be calculated (typically SO2).
            The ratio will be performed against the first specie in the fit window. */
        CFitWindow m_masterFitWindow;

        /** The results from evaluation from the scan using the m_masterFitWindow */
        BasicScanEvaluationResult m_masterResult;
        CPlumeInScanProperty m_masterResultProperties;

        /** Settings for dark-correction */
        const Configuration::CDarkSettings& m_darkSettings;

        /** The fit windows which will be used to estimate the minor specie, typically BrO.
            The name of each fit-window must equal the name of the sub-specie to calculate.
            Being a vector makes it possible to evaluate multiple species in the same round. */
        std::vector<CFitWindow> m_referenceFit;

        /** The settings for how the ratio-calculation should be performed. */
        Configuration::RatioEvaluationSettings m_settings;

        std::vector<Ratio> Run(IScanSpectrumSource& scan, RatioEvaluationDebugInformation& debugInfo);

        // Performs a dark-correction of the provided spectrum, which is assumed to be an average of several spectra in the scan.
        // @throws std::invalid_argument if the dark correction failed.
        void DarkCorrectSpectrum(IScanSpectrumSource& scan, CSpectrum& spectrum) const;

        // Return true if any of the fit windows here requires the calculation of a Ring spectrum
        bool AnyFitWindowRequiresRingSpectrum() const;

        // Validates the settings. @throws std::invalid_argument if the setup is not coherent.s
        void ValidateSetup() const;
    };

    /** Estimates which spectra should be used for a ratio-evaluation, assuming that one should be performed.
        @param referenceSpectra Will be filled with the index of the spectra which should be averaged to a reference spectrum.
        @param inPlumeSpectra Will be filled with the index of the spectra which should be averaged to an in-plume spectrum. */
    void SelectSpectraForRatioEvaluation(const Configuration::RatioEvaluationSettings& settings, const BasicScanEvaluationResult& scanResult, const CPlumeInScanProperty& properties, std::vector<int>& referenceSpectra, std::vector<int>& inPlumeSpectra);

}