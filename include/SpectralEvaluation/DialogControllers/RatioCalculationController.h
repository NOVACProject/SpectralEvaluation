#pragma once

#include <string>
#include <vector>
#include <unordered_map>
#include <SpectralEvaluation/Evaluation/FitWindow.h>
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

class RatioCalculationController
{
public:
    RatioCalculationController();

    // Sets up all the values, and the fit windows, to their default values
    void InitializeToDefault();

    // The list of .pak files to calculate a ratio for
    std::vector<std::string> m_pakfiles;

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

    // Sets up m_so2Window and m_broWindow using m_references, m_so2FitRange and m_broFitRange
    void SetupFitWindows();

    // Performs the evaluation of the next scan (in the list m_pakFiles)
    // This require that SetupFitWindows() has been called since it will use the contents of the fit windows.
    void EvaluateNextScan();

private:

    // The evaluation settings for the major species (typically SO2).
    // The first reference here is used for the ratio calculation.
    // Notice that the list of references here is empty until running (TODO: function name) when they are filled in from m_references.
    novac::CFitWindow m_so2Window;

    // The evaluation settings for the minor species (typically BrO)
    // The first reference here is used for the ratio calculation.
    // Notice that the list of references here is empty until running (TODO: function name) when they are filled in from m_references.
    novac::CFitWindow m_broWindow;

};