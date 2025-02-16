#pragma once

#include <vector>
#include <memory>
#include <SpectralEvaluation/Log.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>

namespace Configuration
{
struct RatioEvaluationSettings;
}

namespace novac
{
class CPlumeInScanProperty;
class CSpectrum;
class BasicScanEvaluationResult;
class IScanSpectrumSource;

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
    std::vector<int> referenceSpectrumIndices;

    // Listing the indices (in the original scan) which were used to create the in plume spectrum.
    std::vector<int> inPlumeSpectrumIndices;

    // Listing the indices (in the original scan) which for some reason were not considered for participating in the in-plume or out-of-plume spectra,
    // together with the reason for not including them.
    std::vector<std::pair<int, std::string>> rejectedSpectrumIndices;

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

    PlumeSpectrumSelector(novac::ILogger& log) : m_log(log) { }

    /**
     * @brief Selects in-plume and out-of-plume spectra from the given scan file with given evaluation result and selection criteria.
     * @param originalScanFile The .pak file where the spectra are found.
     * @param originalScanResult The result of evaluating the provided scan file.
     * @param plumeProperties Describes how much of the plume is visible in the scan, and at what angle.
     * @param mainSpecieIndex The index of the main species (typically SO2) in the originalScanResult.
     * @param errorMessage If present, will be filled with the reason PlumeSpectra is null if it could not be created.
     * @return A created PlumeSpectra struct, or nullptr if none could be created.
    */
    std::unique_ptr<PlumeSpectra> CreatePlumeSpectra(
        novac::LogContext context,
        IScanSpectrumSource& originalScanFile,
        const BasicScanEvaluationResult& originalScanResult,
        const CPlumeInScanProperty& plumeProperties,
        const Configuration::RatioEvaluationSettings& settings,
        const SpectrometerModel& spectrometerModel,
        int mainSpecieIndex = 0);

    /**
     * @brief Selects in-plume and out-of-plume spectra from the given scan file with given evaluation result and selection criteria
     * and saves these to file.
     * @param originalScanFile The .pak file where the spectra are found.
     * @param originalScanResult The result of evaluating the provided scan file.
     * @param plumeProperties Describes how much of the plume is visible in the scan, and at what angle.
     * @param mainSpecieIndex The index of the main species (typically SO2) in the originalScanResult.
     * @param outputDirectory The destination directory where the output file should be saved.
    */
    void CreatePlumeSpectrumFile(
        novac::LogContext context,
        IScanSpectrumSource& originalScanFile,
        const BasicScanEvaluationResult& originalScanResult,
        const CPlumeInScanProperty& plumeProperties,
        const Configuration::RatioEvaluationSettings& settings,
        const SpectrometerModel& spectrometerModel,
        int mainSpecieIndex,
        const std::string& outputDirectory);

private:

    novac::ILogger& m_log;

    int m_mainSpecieIndex = 0;

    // Helper structure, collecting together data which we need for selecting spectra in the scan.
    struct InitialEvaluationData
    {
        int indexInScan = 0;
        double scanAngle = 0.0;
        double offsetCorrectedColumn = 0.0;
        double peakSaturation = 0.0;
        double peakSaturationAfterDarkCorrection = 0.0;
    };

    /**
     * @brief Returns true if the provided originalScanResult represents a scan which is suitable for performing ratio evaluation
     * @param originalScanResult The evaluation result.
     * @param properties The properties of the scan, must have center and completeness filled in.
     * @return True if the scan is suitable.
    */
    bool IsSuitableScanForRatioEvaluation(
        const BasicScanEvaluationResult& originalScanResult,
        const CPlumeInScanProperty& properties,
        const Configuration::RatioEvaluationSettings& settings,
        std::string* errorMessage = nullptr);

    /**
     * @brief Checks the provided evaluated scan using default settings and returns the indices of the spectra which can
     * be saved as in-plume and out-of-plume spectra.This will use the provided CPlumeInScanProperty to find out
     *   if the scan sees the plume at all and where the edges of the plume are located.
     * @param scanFile The .pak file to check.
     * @param darkSpectrum The dark spectrum to be used for EVERY spectrum in the scan.
     * @param originalScanResult The evaluation result of the scan
     * @param properties The properties of the plume, must be filled in with completeness, plumeHalfLow and plumeHalfHigh.
     * @param spectrometerModel Model of the spectrometer which collected this scan.
     * @param referenceSpectra Will on return be filled with the indices of the spectra to be used as out-of-plume spectra.
     * @param inPlumeSpectra Will on return be filled with the indices of the spectra to be used as in-plume spectra.
     * @return true if the selection was successful.
    */
    bool SelectSpectra(
        const std::vector< InitialEvaluationData>& evaluationData,
        const CPlumeInScanProperty& properties,
        const Configuration::RatioEvaluationSettings& settings,
        std::vector<int>& referenceSpectra,
        std::vector<int>& inPlumeSpectra,
        std::string* errorMessage = nullptr);

    std::vector<int> FindSpectraInPlume(
        const std::vector< InitialEvaluationData>& evaluationData,
        const CPlumeInScanProperty& properties,
        const Configuration::RatioEvaluationSettings& settings);

    std::vector<int> FindSpectraOutOfPlume(
        const std::vector< InitialEvaluationData>& evaluationData,
        const Configuration::RatioEvaluationSettings& settings,
        const std::vector<int>& inPlumeProposal,
        double minimumColumnInPlume);

    /**
     * @brief Return the peak intensity of the provided spectrum before and after dark-correction.
     * @param spectrum The spectrum to verify, not dark-corrected
     * @param darkSpectrum The dark-spectrum corresponding to the provded spectrum.
     * @param spectrometerModel The model of the spectrometer which collected the measurement.
     * @param peakSaturation Will on return equal the highest intensity of the spectrum
     * @param peakSaturationAfterDarkCorrection Will on return equal the highest intensity of the spectrum
     * after the dark spectrum has been subtracted. */
    void GetIntensityOfSpectrum(
        const CSpectrum& spectrum,
        const CSpectrum& darkSpectrum,
        const SpectrometerModel& spectrometerModel,
        double& peakSaturation,
        double& peakSaturationAfterDarkCorrection) const;

    void ExtractEvaluationData(
        novac::LogContext context,
        IScanSpectrumSource& originalScanFile,
        const CSpectrum& darkSpectrum,
        const BasicScanEvaluationResult& originalScanResult,
        const Configuration::RatioEvaluationSettings& settings,
        const SpectrometerModel& spectrometermodel,
        std::vector< InitialEvaluationData>& evaluationData,
        std::vector<std::pair<int, std::string>>& rejectedIndices) const;

    // Helper method for finding data points.
    double ColumnAtScanIndex(const std::vector< InitialEvaluationData> evaluationData, int scanIndex) const;
};
}
