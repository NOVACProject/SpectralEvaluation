#pragma once

#include <vector>
#include <random>

// -----------------------------------------------------------------------------------------------------------------------------
// - This header contains a helper struct used to perform the wavelength calibration of a spectrometer using a ransac approach -
// -----------------------------------------------------------------------------------------------------------------------------


namespace novac
{

/// <summary>
/// A Correspondence is an essential part in the ransac calibration routine, it
/// represents a connection between a point in the measured spectrum and a point
/// in the theoretical fraunhofer spectrum.
/// </summary>
struct Correspondence
{
    Correspondence()
        : measuredIdx(0), theoreticalIdx(0), error(0.0)
    { }

    Correspondence(size_t measured, size_t theoretical, double error)
        : measuredIdx(measured), theoreticalIdx(theoretical), error(error)
    { }

    /// <summary>
    /// The index of the keypoint in the measured spectrum
    /// </summary>
    size_t measuredIdx = 0;

    /// <summary>
    /// The value of the keypoint in the measured spectrum (pixel, in the case of wavelength calibration)
    /// </summary>
    double measuredValue = 0.0;

    /// <summary>
    /// The index of the keypoint in the theoretical spectrum
    /// </summary>
    size_t theoreticalIdx = 0;

    /// <summary>
    /// The value of the keypoint in the theoretical spectrum (wavelength, in the case of wavelength calibration)
    /// </summary>
    double theoreticalValue = 0.0;

    /// <summary>
    /// An error measure between the keypoints in the two spectra, lower is better
    /// </summary>
    double error = 0.0;
};

/** Returns true if the provided set of correspondences all have unique measuredIdx */
bool AllMeasuredPointsAreUnique(const std::vector<Correspondence>& correspondences);

/** Returns true if the provided set of correspondences all have unique theoreticalIdx */
bool AllTheoreticalPointsAreUnique(const std::vector<Correspondence>& correspondences);

/** Lists all possible combinations of the provided vector of Correspondences such that
    each combination contains 'number' elements and each combination have all unique measureIdx and all unique theoreticalIdx. */
std::vector<std::vector<Correspondence>> ListAllPossibleKeypointCombinations(size_t number, const std::vector<Correspondence>& allCorrespondences);

/** Selects one combination of the provided vector of Correspondences such that it contains
    'number' elements and has all unique measureIdx and all unique theoreticalIdx.
    The provided randomGenerator will be used to select elements from 'allCorrespondences'.
    @throws std::invalid_argument if 'number' is smaller than the number of elements in 'allCorrespondences' */
void SelectMaybeInliers(size_t number, const std::vector<Correspondence>& allCorrespondences, std::mt19937& randomGenerator, std::vector<Correspondence>& result);

/** Finds the index of the first Correspondence which has theoreticalIdx equal to theoreticalIdxOfEntryToFind
    Returns -1 if no such Correspondence can be found. */
int FindCorrespondenceWithTheoreticalIdx(const std::vector<Correspondence>& allEntries, size_t theoreticalIdxOfEntryToFind);

/** Calculates and returns the measured pixel-distance between the correspondence with the smallest 
    measuredIdx and the correspondence with the largest measuredIdx. */
double MeasurePixelSpan(const std::vector<Correspondence>& correspondences);

/** Sorts the provided set of correspondences into a vector where the element at position 'N' lists 
    all the correspondences for the measured keypoint 'N' */
std::vector<std::vector<Correspondence>> ArrangeByMeasuredKeypoint(const std::vector<Correspondence>& allCorrespondences);

/** Returns a vector of equal length to 'selectedValues' where element ii is true if the correspondence selectedValues[ii]
    is present in allValues */
std::vector<bool> ListInliers(const std::vector<Correspondence>& selectedValues, const std::vector<Correspondence>& allValues);


}
