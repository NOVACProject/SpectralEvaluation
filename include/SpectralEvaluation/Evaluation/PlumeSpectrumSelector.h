#pragma once
#include <vector>

#include <SpectralEvaluation/Spectra/Spectrum.h>

namespace novac
{
    class CPlumeInScanProperty;
    class CSpectrum;
    class BasicScanEvaluationResult;
    class CScanFileHandler;

    /** Struct used to store the selected in-plume and out-of-plume spectra
        from one single scan, with the intention that these can be used later for 
        evaluating the ratio between gases in the plume. */
    struct PlumeSpectra
    {
        // The out-of-plume reference spectrum. Not dark-corrected.
        std::unique_ptr<CSpectrum> referenceSpectrum;

        // The in-plume spectrum. Not dark-corrected.
        std::unique_ptr<CSpectrum> inPlumeSpectrum;

        // The dark spectrum.
        std::unique_ptr<CSpectrum> darkSpectrum;

        // Listing the indices (in the original scan) which were used to create the reference spectrum.
        std::vector<size_t> referenceSpectrumIndices;

        // Listing the indices (in the original scan) which were used to create the in plume spectrum.
        std::vector<size_t> inPlumeSpectrumIndices;

        // The CSpectrumInfo of the sky-spectrum of the scan. Gives info on e.g. device and start time.
        CSpectrumInfo skySpectrumInfo;
    };

    // This is a helper class for selecting 'in plume' and 'out of plume' spectra
    //  (based on some basic selection criteria) with the intention of creating 
    //  reference spectra for performing e.g. ratio evaluations or detailed 
    //  spectral analysis of compounds in the plume.
    class PlumeSpectrumSelector
    {
    public:
        struct PlumeSpectrumSelectionSettings
        {
            int minNumberOfSpectraInPlume = 4;

            int maxNumberOfSpectraInPlume = 10;

            int numberOfSpectraOutsideOfPlume = 10;

            // Lowest allowed angle to include, in degrees from zenith. Used to exclude spectra too close to the horizon.
            double minimumScanAngle = -75.0;

            // Highest allowed angle to include, in degrees from zenith. Used to exclude spectra too close to the horizon.
            double maximumScanAngle = +75.0;

            // The maximum saturation ratio (intensity / maximum intensity of spectrometer)
            double maxSaturationRatio = 0.88;

            // The minimum saturation ratio (intensity / maximum intensity of spectrometer)
            double minSaturationRatio = 0.12;
        };

        /**
         * @brief Selects in-plume and out-of-plume spectra from the given scan file with given evaluation result and selection criteria.
         * @param originalScanFile The .pak file where the spectra are found.
         * @param scanResult The result of evaluating the provided scan file.
         * @param selectionProperties Criteria for how to select the scans.
         * @param mainSpecieIndex The index of the main species (typically SO2) in the scanResult.
         * @return A created PlumeSpectra struct, or nullptr if none could be created.
        */
        std::unique_ptr<PlumeSpectra> CreatePlumeSpectra(
            CScanFileHandler& originalScanFile,
            const BasicScanEvaluationResult& scanResult,
            const CPlumeInScanProperty& selectionProperties,
            int mainSpecieIndex = 0);

        /**
         * @brief Selects in-plume and out-of-plume spectra from the given scan file with given evaluation result and selection criteria
         * and saves these to file.
         * @param originalScanFile The .pak file where the spectra are found.
         * @param scanResult The result of evaluating the provided scan file.
         * @param selectionProperties Criteria for how to select the scans.
         * @param mainSpecieIndex The index of the main species (typically SO2) in the scanResult.
         * @param outputDirectory The destination directory where the output file should be saved.
        */
        void CreatePlumeSpectrumFile(
            CScanFileHandler& originalScanFile,
            const BasicScanEvaluationResult& scanResult,
            const CPlumeInScanProperty& selectionProperties,
            int mainSpecieIndex,
            const std::string& outputDirectory);

    private:
        double m_maximumSpectrometerIntensity = 4095.0;

        int m_mainSpecieIndex = 0;

        PlumeSpectrumSelectionSettings m_settings;

        bool IsSuitableScanForRatioEvaluation(
            const CSpectrum& skySpectrum,
            const CSpectrum& darkSpectrum,
            const BasicScanEvaluationResult& scanResult,
            const CPlumeInScanProperty& properties);

        // Checks the provided evaluated scan using default settings 
        //  and returns the indices of the spectra which can 
        //  be saved as in-plume and out-of-plume spectra.
        // This will use the provided CPlumeInScanProperty to find out
        //   if the scan sees the plume at all and where the edges of 
        //   the plume are located.
        // If this check fails, then the two vectors are empty.
        void SelectSpectra(
            CScanFileHandler& scanFile,
            const CSpectrum& darkSpectrum,
            const BasicScanEvaluationResult& scanResult,
            const CPlumeInScanProperty& properties,
            std::vector<size_t>& referenceSpectra,
            std::vector<size_t>& inPlumeSpectra);

        std::vector<size_t> FindSpectraInPlume(
            const BasicScanEvaluationResult& scanResult,
            const CPlumeInScanProperty& properties);

        std::vector<size_t> FindSpectraOutOfPlume(
            CScanFileHandler& scanFile,
            const CSpectrum& darkSpectrum,
            const BasicScanEvaluationResult& scanResult,
            const std::vector<size_t>& inPlumeProposal);

        std::vector<size_t> FilterSpectraUsingIntensity(
            const std::vector<size_t>& proposedIndices,
            CScanFileHandler& scanFile,
            const CSpectrum& darkSpectrum);

        bool SpectrumFulfillsIntensityRequirement(
            const CSpectrum& spectrum,
            const CSpectrum& darkSpectrum);
    };
}
