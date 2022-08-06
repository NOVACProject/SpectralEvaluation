#pragma once

#include <memory>
#include <string>
#include <vector>
#include <SpectralEvaluation/Evaluation/RatioEvaluation.h>
#include <SpectralEvaluation/Spectra/WavelengthRange.h>

// Building a set of standard DOAS species
enum class StandardDoasSpecie
{
    NotSet,
    SO2,
    BRO,
    O3,
    RING,
    RING_LAMBDA4,
    O4,
    CH20,
    NO2
};

// The structure ReferenceForRatioCalculation represents a reference file which
// may be included in the DOAS fit for performing a ratio calculation.
// This is intended to be an easier way of allowing the user to control the user interface
struct ReferenceForRatioCalculation
{
    ReferenceForRatioCalculation() = default;

    ReferenceForRatioCalculation(StandardDoasSpecie specieEnum, const std::string& name, const std::string& path, bool includeInMajor, bool includeInMinor, bool automaticallyCalculate) :
        specie(specieEnum), m_name(name), m_path(path), m_includeInMajor(includeInMajor), m_includeInMinor(includeInMinor), m_automaticallyCalculate(automaticallyCalculate)
    {
    }

    // The type of specie.
    StandardDoasSpecie specie;

    // The (user given) name of the reference. Automatically set to a reasonable reference name.
    std::string m_name = "";

    // The full file path to the location of the reference file.
    std::string m_path = "";

    // Set to true if this reference should be included in the DOAS fit calculation of the Major window (SO2)
    bool m_includeInMajor = true;

    // Set to true if this reference should be included in the DOAS fit calculation of the Minor window (BrO)
    bool m_includeInMinor = true;

    // Set to true to automatically calculate this reference from the measured data (only possible for Ring and RingxLambda4)
    bool m_automaticallyCalculate = false;

    // @return true if this reference can be automatically calculated from the measured data (only possible for Ring and RingxLambda4)
    bool CanBeAutomaticallyCalculated() const { return specie == StandardDoasSpecie::RING || specie == StandardDoasSpecie::RING_LAMBDA4; }
};

// The required input for the ratio calculation.
struct RatioCalculationFitSetup
{
    // The evaluation settings for the major species (typically SO2).
    // The first reference here is used for the ratio calculation.
    novac::CFitWindow so2Window;

    // The evaluation settings for the minor species (typically BrO)
    // The first reference here is used for the ratio calculation.
    novac::CFitWindow broWindow;
};

// The result of the ratio calculation, including debug information
struct RatioCalculationResult
{
    // The final output
    novac::Ratio ratio;

    std::string filename;

    // The initial evaluation, faster evaluations on each spectrum in the scan.
    novac::BasicScanEvaluationResult initialEvaluation;

    // The calculated properties of the plume (from initialEvaluation). Includes completeness and offset.
    novac::CPlumeInScanProperty plumeInScanProperties;

    // Detailed information the evaluation
    novac::RatioEvaluationDebugInformation debugInfo;
};

class RatioCalculationController
{
public:
    RatioCalculationController();

    // Sets up all the values, and the fit windows, to their default values
    void InitializeToDefault();

    // Sets up this controller with the provided list of .pak files, each containing one scan.
    void SetupPakFileList(const std::vector<std::string>& pakFiles);

    // Retrieves the current list of .pak files.
    std::vector<std::string> ListPakFiles() const;

    // Returns the number of .pak files in the current setup (same as ListPakFiles().size())
    size_t NumberOfPakFilesInSetup() const;

    // The list of references which are to be included in the two evaluations.
    std::vector<ReferenceForRatioCalculation> m_references;

    // The wavelength range for the major species window
    novac::WavelengthRange m_so2FitRange;

    // The polynomial order to use in the SO2 fit
    int m_so2PolynomialOrder = 3;

    // The wavelength range for the minor species window
    novac::WavelengthRange m_broFitRange;

    // The polynomial order to use in the BrO fit
    int m_broPolynomialOrder = 3;

    // Settings for how to select the in-plume and out-of-plume spectra.
    novac::RatioEvaluationSettings m_ratioEvaluationSettings;

    // The last result, as calculated by EvaluateNextScan or EvaluateScan.
    RatioCalculationResult m_lastResult;

    // Sets up the SO2 and BrO fit windows using m_references, m_so2FitRange and m_broFitRange
    // @throws std::invalid_argument if the setup is incorrect
    std::shared_ptr<RatioCalculationFitSetup> SetupFitWindows();

    // @return true if there are more scans available to evaluate in the list of .pak-files.
    bool HasMoreScansToEvaluate() const;

    // Performs the evaluation of the next scan (in the list m_pakFiles)
    // If a prior call to HasMoreScansToEvaluate() would return false, then this will return an empty result.
    // This require that SetupFitWindows() has been called since it will use the contents of the fit windows.
    RatioCalculationResult EvaluateNextScan(std::shared_ptr<RatioCalculationFitSetup> ratioFitWindows);

    // Performs the initial evaluation (typically of SO2) for the provided scan.
    // Produces an initial result, when can then be used to see if the scan is suitable for ratio calculations, and if so to select in-plume and out-of-plume spectra.
    // @throws std::invalid_argument if the setup is such that the scan cannot be evaluated.
    novac::BasicScanEvaluationResult DoInitialEvaluation(novac::IScanSpectrumSource& scan, std::shared_ptr<RatioCalculationFitSetup> ratioFitWindows);

    // Performs the evaluation of the given scan.
    // This require that SetupFitWindows() has been called since it will use the contents of the fit windows.
    // @return a vector containing the calculated ratio if the evaluation succeeded.
    RatioCalculationResult EvaluateScan(novac::IScanSpectrumSource& scan, const novac::BasicScanEvaluationResult& initialResult, std::shared_ptr<RatioCalculationFitSetup> ratioFitWindows);

private:

    // The list of .pak files to calculate a ratio for
    std::vector<std::string> m_pakfiles;

    // The index (into m_pakFiles) which we're currently at.
    int m_currentPakFileIdx = 0;

    // Verifies that the references have been setup as expected. Throws invalid_argument if not.
    void VerifyReferenceSetup();
};

// Helper method, picks out the reference with the given specie and returns a pointer to it.
// Returns nullptr if there is no reference with the given specie.
ReferenceForRatioCalculation* GetReferenceFor(RatioCalculationController& controller, StandardDoasSpecie specie);
