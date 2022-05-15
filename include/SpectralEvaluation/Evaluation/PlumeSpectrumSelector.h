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

    /** This is a helper class for selecting 'in plume' and 'out of plume' spectra
        (based on some basic selection criteria) with the intention of creating
        reference spectra for performing e.g. ratio evaluations or detailed
        spectral analysis of compounds in the plume. */
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
         * @param plumeProperties Describes how much of the plume is visible in the scan, and at what angle.
         * @param mainSpecieIndex The index of the main species (typically SO2) in the scanResult.
         * @return A created PlumeSpectra struct, or nullptr if none could be created.
        */
        std::unique_ptr<PlumeSpectra> CreatePlumeSpectra(
            CScanFileHandler& originalScanFile,
            const BasicScanEvaluationResult& scanResult,
            const CPlumeInScanProperty& plumeProperties,
            int mainSpecieIndex = 0);

        /**
         * @brief Selects in-plume and out-of-plume spectra from the given scan file with given evaluation result and selection criteria
         * and saves these to file.
         * @param originalScanFile The .pak file where the spectra are found.
         * @param scanResult The result of evaluating the provided scan file.
         * @param plumeProperties Describes how much of the plume is visible in the scan, and at what angle.
         * @param mainSpecieIndex The index of the main species (typically SO2) in the scanResult.
         * @param outputDirectory The destination directory where the output file should be saved.
        */
        void CreatePlumeSpectrumFile(
            CScanFileHandler& originalScanFile,
            const BasicScanEvaluationResult& scanResult,
            const CPlumeInScanProperty& plumeProperties,
            int mainSpecieIndex,
            const std::string& outputDirectory);

    private:

        int m_mainSpecieIndex = 0;

        PlumeSpectrumSelectionSettings m_settings;

        /**
         * @brief Returns true if the provided scanResult represents a scan which is suitable for performing ratio evaluation
         * @param skySpectrum The sky spectrum of the scan.
         * @param darkSpectrum The dark spectrum of the scan.
         * @param scanResult The evaluation result.
         * @param properties The properties of the scan, must have center and completeness filled in.
         * @param spectrometerModel The model of the spectrometer which collected the measurement.
         * @return True if the scan is suitable.
        */
        bool IsSuitableScanForRatioEvaluation(
            const CSpectrum& skySpectrum,
            const CSpectrum& darkSpectrum,
            const BasicScanEvaluationResult& scanResult,
            const CPlumeInScanProperty& properties,
            const SpectrometerModel& spectrometerModel);

        /**
         * @brief Checks the provided evaluated scan using default settings and returns the indices of the spectra which can
         * be saved as in-plume and out-of-plume spectra.This will use the provided CPlumeInScanProperty to find out
         *   if the scan sees the plume at all and where the edges of the plume are located.
         * @param scanFile The .pak file to check.
         * @param darkSpectrum The dark spectrum to be used for EVERY spectrum in the scan.
         * @param scanResult The evaluation result of the scan
         * @param properties The properties of the plume, must be filled in with completeness, plumeHalfLow and plumeHalfHigh.
         * @param spectrometerModel Model of the spectrometer which collected this scan.
         * @param referenceSpectra Will on return be filled with the indices of the spectra to be used as out-of-plume spectra.
         * @param inPlumeSpectra Will on return be filled with the indices of the spectra to be used as in-plume spectra.
        */
        void SelectSpectra(
            CScanFileHandler& scanFile,
            const CSpectrum& darkSpectrum,
            const BasicScanEvaluationResult& scanResult,
            const CPlumeInScanProperty& properties,
            const SpectrometerModel& spectrometerModel,
            std::vector<size_t>& referenceSpectra,
            std::vector<size_t>& inPlumeSpectra);

        std::vector<size_t> FindSpectraInPlume(
            const BasicScanEvaluationResult& scanResult,
            const CPlumeInScanProperty& properties);

        std::vector<size_t> FindSpectraOutOfPlume(
            CScanFileHandler& scanFile,
            const CSpectrum& darkSpectrum,
            const BasicScanEvaluationResult& scanResult,
            const SpectrometerModel& spectrometerModel,
            const std::vector<size_t>& inPlumeProposal);

        std::vector<size_t> FilterSpectraUsingIntensity(
            const std::vector<size_t>& proposedIndices,
            CScanFileHandler& scanFile,
            const CSpectrum& darkSpectrum,
            const SpectrometerModel& spectrometerModel);

        /**
         * @brief Return true if the provided spectrum has a good intensity for being included in a ratio evaluation.
         * @param spectrum The spectrum to verify, not dark-corrected
         * @param darkSpectrum The dark-spectrum corresponding to the provded spectrum.
         * @param spectrometerModel The model of the spectrometer which collected the measurement.
         * @return True if the provided spectrum has suitable properties.
        */
        bool SpectrumFulfillsIntensityRequirement(
            const CSpectrum& spectrum,
            const CSpectrum& darkSpectrum,
            const SpectrometerModel& spectrometerModel);
    };
}
