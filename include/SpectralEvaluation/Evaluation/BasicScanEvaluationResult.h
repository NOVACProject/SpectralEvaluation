#pragma once

#include <memory>
#include <vector>
#include <SpectralEvaluation/Molecule.h>
#include <SpectralEvaluation/NovacEnums.h>
#include <SpectralEvaluation/Evaluation/EvaluationResult.h>
#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>
#include <SpectralEvaluation/Spectra/SpectrumInfo.h>

namespace novac
{
/** An instance of the class BasicScanEvaluationResult contains a minimal amount of information
    regarding the evaluation of a full scan using one fit-window.
    This forms the base-class for the classes CScanResult in NovacProgram and CScanResult in NovacPPP. */
class BasicScanEvaluationResult
{
public:

#pragma region Properties on the spectra in the scan itself

    /** The results of evaluating the spectra.
        There is one evaluation result for each spectrum in the scan.
        These are ordered such that m_spec[i] contains the result for spectrum #i in the scan. */
    std::vector<CEvaluationResult> m_spec;

    /** General information about the collected spectra.
        There is one entry here for each spectrum in the scan, and element #i must match
         element #i in the vector m_spec. */
    std::vector<CSpectrumInfo> m_specInfo;

    /** information about the sky-spectrum used */
    CSpectrumInfo m_skySpecInfo;

    /** information about the dark-spectrum used, if any */
    CSpectrumInfo m_darkSpecInfo;

    /** information about the offset-spectrum used, if any */
    CSpectrumInfo m_offsetSpecInfo;

    /** information about the dark-current spectrum used, if any */
    CSpectrumInfo m_darkCurSpecInfo;

    /** A path where the .pak file associated with this result resides.
        This makes it possible to re-read the scan and analyze it again. */
    std::string m_path = "";

    /** A list of which spectra were corrupted and could not be evaluated.
        There is one entry here for each spectrum which wasn't evaluated. */
    std::vector<unsigned int> m_corruptedSpectra;

    /** The type of the instrument used for this scan */
    NovacInstrumentType m_instrumentType = NovacInstrumentType::Gothenburg;

    /** Flag to signal if this is a wind measurement, a scan, or something else. */
    MeasurementMode m_measurementMode = MeasurementMode::Unknown;

#pragma endregion Properties on the spectra in the scan itself

#pragma region Calculated properties on the scan

    /** This contains the parameters of the plume that is seen in this scan,
        such as the completeness or the centre angle of the plume. */
    novac::CPlumeInScanProperty m_plumeProperties;

#pragma endregion Calculated properties on the scan

    /** Appends the provided result to the list of calculated results */
    int AppendResult(const novac::CEvaluationResult& evalRes, const novac::CSpectrumInfo& specInfo);

    /** Removes the spectrum number 'specIndex' from the list of calcualted results */
    int RemoveResult(unsigned int specIndex);

    /** Intializes the memory arrays to have, initially, space for 'size' spectra. */
    void InitializeArrays(long size);

    /** Returns the index of the specie with the provided name.
        E.g. checks the scan result for which index corresponds to SO2.
        The comparison is done ignoring case.
        @return -1 if the specie could not be found */
    int GetSpecieIndex(const std::string& specieName) const;

    /** Returns the index of the specie.
        @return -1 if the specie could not be found */
    int GetSpecieIndex(Molecule molecule) const;

    #pragma region Getting properties of the scan

    unsigned long NumberOfEvaluatedSpectra() const { return m_specNum; }

    // Attempts to retrieve the location of the instrument from the SpectrumInfo.
    // If no GPS data is found then the retuned location has both latitude and longitude = 0.0.
    CGPSData GetLocation() const;

    #pragma endregion Getting properties of the scan

protected:

    /** The number of evaluated spectra. */
    unsigned long m_specNum = 0;

};

/** @return all the evaluated columns for the specie with the provided index.
    @return an empty vector if result is empty specieIndex is invalid. */
std::vector<double> GetColumns(const BasicScanEvaluationResult& result, int specieIndex);

/** @return the errors of all the evaluated columns for the specie with the provided index.
    @return an empty vector if result is empty specieIndex is invalid. */
std::vector<double> GetColumnErrors(const BasicScanEvaluationResult& result, int specieIndex);

/** Tries to find a plume in the provided scan result.
    If the plume is found then the calculated properties are returned otherwise this returns nullptr.
    This will fill in the scan offset, the plume completeness and the plume positions. */
std::unique_ptr<novac::CPlumeInScanProperty> CalculatePlumeProperties(const BasicScanEvaluationResult& result, const Molecule& specie, std::string& message);

/** Checks the kind of measurement mode of the provided scan */
novac::MeasurementMode CheckMeasurementMode(const BasicScanEvaluationResult& result);

bool IsFluxMeasurement(const BasicScanEvaluationResult& result);

bool IsWindMeasurement(const BasicScanEvaluationResult& result);

bool IsWindMeasurement_Gothenburg(const BasicScanEvaluationResult& result);

bool IsWindMeasurement_Heidelberg(const BasicScanEvaluationResult& result);

bool IsStratosphereMeasurement(const BasicScanEvaluationResult& result);

bool IsDirectSunMeasurement(const BasicScanEvaluationResult& result);

bool IsLunarMeasurement(const BasicScanEvaluationResult& result);

bool IsCompositionMeasurement(const BasicScanEvaluationResult& result);

}