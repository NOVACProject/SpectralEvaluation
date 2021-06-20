#include <SpectralEvaluation/Calibration/WavelengthCalibrationByRansac.h>

#include <iostream>
#include <omp.h>
#include <random>
#include <set>
#include <stdexcept>
#include <SpectralEvaluation/Math/PolynomialFit.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Spectra/SpectrumUtils.h>
#include <SpectralEvaluation/VectorUtils.h>

#include <SpectralEvaluation/Evaluation/BasicMath.h>

#undef min
#undef max

namespace novac
{

// ---------------------------------- RansacWavelengthCalibrationResult ----------------------------------
RansacWavelengthCalibrationResult::RansacWavelengthCalibrationResult() :
    bestFittingModelCoefficients(41),
    modelPolynomialOrder(3),
    highestNumberOfInliers(0U),
    smallestError(std::numeric_limits<double>::max()),
    numberOfPossibleCorrelations(0U)
{
}

RansacWavelengthCalibrationResult::RansacWavelengthCalibrationResult(size_t polynomialOrder) :
    bestFittingModelCoefficients(polynomialOrder + 1),
    modelPolynomialOrder(polynomialOrder),
    highestNumberOfInliers(0U),
    smallestError(std::numeric_limits<double>::max()),
    numberOfPossibleCorrelations(0U)
{
}

RansacWavelengthCalibrationResult::RansacWavelengthCalibrationResult(const RansacWavelengthCalibrationResult& other) :
    bestFittingModelCoefficients(begin(other.bestFittingModelCoefficients), end(other.bestFittingModelCoefficients)),
    modelPolynomialOrder(other.modelPolynomialOrder),
    highestNumberOfInliers(other.highestNumberOfInliers),
    correspondenceIsInlier(other.correspondenceIsInlier),
    smallestError(other.smallestError),
    numberOfPossibleCorrelations(other.numberOfPossibleCorrelations)
{
}

RansacWavelengthCalibrationResult& RansacWavelengthCalibrationResult::operator=(const RansacWavelengthCalibrationResult& other)
{
    this->modelPolynomialOrder = other.modelPolynomialOrder;
    this->bestFittingModelCoefficients = std::vector<double>(begin(other.bestFittingModelCoefficients), end(other.bestFittingModelCoefficients));
    this->highestNumberOfInliers = other.highestNumberOfInliers;
    this->correspondenceIsInlier = std::vector<bool>(begin(other.correspondenceIsInlier), end(other.correspondenceIsInlier));
    this->smallestError = other.smallestError;
    this->numberOfPossibleCorrelations = other.numberOfPossibleCorrelations;

    return *this;
}

// ---------------------------------- Free functions used to select the correspondences ----------------------------------

double MeasureCorrespondenceError(
    const CSpectrum& measuredSpectrum,
    double pixelInMeasuredSpectrum,
    const CSpectrum& fraunhoferSpectrum,
    double pixelInFraunhoferSpectrum,
    const CorrespondenceSelectionSettings& settings)
{
    // extract a small region around 'idxOfMeasuredSpectrum' in the measured spectrum
    const size_t measuredLocationStart = (size_t)std::max(0.0, std::min((double)(measuredSpectrum.m_length - settings.pixelRegionSizeForCorrespondenceErrorMeasurement), pixelInMeasuredSpectrum - settings.pixelRegionSizeForCorrespondenceErrorMeasurement * 0.5));
    std::vector<double> measuredRegion{ measuredSpectrum.m_data + measuredLocationStart, measuredSpectrum.m_data + measuredLocationStart + settings.pixelRegionSizeForCorrespondenceErrorMeasurement };
    RemoveMean(measuredRegion);

    // likewise, extract a small region around 'fraunhoferLcation' in the fraunhofer spectrum
    const size_t fraunhoferLocationStart = (size_t)std::max(0.0, std::min((double)(fraunhoferSpectrum.m_length - settings.pixelRegionSizeForCorrespondenceErrorMeasurement), pixelInFraunhoferSpectrum - settings.pixelRegionSizeForCorrespondenceErrorMeasurement * 0.5));
    std::vector<double> fraunhoferRegion{ fraunhoferSpectrum.m_data + fraunhoferLocationStart, fraunhoferSpectrum.m_data + fraunhoferLocationStart + settings.pixelRegionSizeForCorrespondenceErrorMeasurement };
    RemoveMean(fraunhoferRegion);

    return SumOfSquaredDifferences(measuredRegion, fraunhoferRegion);
}

std::vector<novac::Correspondence> ListPossibleCorrespondences(
    const std::vector<novac::SpectrumDataPoint>& measuredKeypoints,
    const CSpectrum& measuredSpectrum,
    const std::vector<novac::SpectrumDataPoint>& fraunhoferKeypoints,
    const CSpectrum& fraunhoferSpectrum,
    const CorrespondenceSelectionSettings& correspondenceSettings)
{
    // to reduce the impact of the noise present, we do a low-pass filtering of the measured spectrum first
    std::unique_ptr<CSpectrum> lowPassFilteredSpectrum = std::make_unique<CSpectrum>(measuredSpectrum);
    ::CBasicMath math;
    math.LowPassBinomial(lowPassFilteredSpectrum->m_data, lowPassFilteredSpectrum->m_length, 5);

    std::vector<novac::Correspondence> possibleCorrespondences;
    std::vector<double> errors; // keep track of the errors, we will use this later to filter the correspondences
    for (size_t ii = 0; ii < measuredKeypoints.size(); ii++)
    {
        if (measuredKeypoints[ii].pixel >= correspondenceSettings.measuredPixelStart && measuredKeypoints[ii].pixel <= correspondenceSettings.measuredPixelStop)
        {
            for (size_t jj = 0; jj < fraunhoferKeypoints.size(); jj++)
            {
                if (measuredKeypoints[ii].type == fraunhoferKeypoints[jj].type &&
                    std::abs(measuredKeypoints[ii].pixel - fraunhoferKeypoints[jj].pixel) <= correspondenceSettings.maximumPixelDistanceForPossibleCorrespondence)
                {
                    novac::Correspondence corr;
                    corr.measuredIdx = ii;
                    corr.measuredValue = measuredKeypoints[ii].pixel;
                    corr.theoreticalIdx = jj;
                    corr.theoreticalValue = fraunhoferKeypoints[jj].wavelength;
                    corr.error = novac::MeasureCorrespondenceError(*lowPassFilteredSpectrum, measuredKeypoints[ii].pixel, fraunhoferSpectrum, fraunhoferKeypoints[jj].pixel, correspondenceSettings);
                    possibleCorrespondences.push_back(corr);

                    errors.push_back(corr.error);
                }
            }
        }
    }

    if (possibleCorrespondences.size() == 0)
    {
        return possibleCorrespondences;
    }

    // Filter out correspondences by keeping only the ones which have an error < median(error)
    std::sort(begin(errors), end(errors));
    size_t idx = std::min(errors.size() - 1, (size_t)(correspondenceSettings.percentageOfCorrespondencesToSelect * errors.size()));
    const double errorLimit = errors[idx];

    std::vector<novac::Correspondence> filteredCorrespondences;
    filteredCorrespondences.reserve(possibleCorrespondences.size());
    for (const novac::Correspondence& corr : possibleCorrespondences)
    {
        if (corr.error <= errorLimit)
        {
            filteredCorrespondences.push_back(corr);
        }
    }

    return filteredCorrespondences;
}

// ---------------------------------- Free functions used to run the calibration ----------------------------------
void SelectMaybeInliers(size_t number, const std::vector<Correspondence>& allCorrespondences, std::mt19937& randomGenerator, std::vector<Correspondence>& result)
{
    if (number == allCorrespondences.size())
    {
        result = allCorrespondences;
        return;
    }
    else if (number > allCorrespondences.size())
    {
        throw std::invalid_argument("Cannot select a larger number of inliers than the entire data set.");
    }

    result.resize(number);
    std::uniform_int_distribution<size_t> uniform_dist(0, allCorrespondences.size() - 1);

    // "Randomly" select correspondences with the following restrictions
    //  1) A given measured or theoretical point must not be present twice in the selected result
    std::set<size_t> selectedIndices;
    std::set<size_t> selecedMeasuredPoints;
    std::set<size_t> selecedTheoreticalPoints;
    while (selectedIndices.size() < number)
    {
        const size_t suggestedIdx = uniform_dist(randomGenerator);

        if (selecedMeasuredPoints.find(allCorrespondences[suggestedIdx].measuredIdx) == selecedMeasuredPoints.end() &&
            selecedTheoreticalPoints.find(allCorrespondences[suggestedIdx].theoreticalIdx) == selecedTheoreticalPoints.end())
        {
            selecedMeasuredPoints.insert(allCorrespondences[suggestedIdx].measuredIdx);
            selecedTheoreticalPoints.insert(allCorrespondences[suggestedIdx].theoreticalIdx);
            selectedIndices.insert(suggestedIdx);
        }
    }


    size_t resultIdx = 0;
    for (size_t idx : selectedIndices)
    {
        result[resultIdx] = allCorrespondences[idx];
        ++resultIdx;
    }
}

double PolynomialValueAt(const std::vector<double>& coefficients, double x)
{
    if (coefficients.size() == 1)
    {
        return coefficients[0];
    }

    auto it = coefficients.rbegin();
    double result = *it;
    ++it;
    for (; it != coefficients.rend(); ++it)
    {
        result = x * result + *it;
    }
    return result;
}

inline double CalculateValueAt(const std::vector<double>& coefficients, double x)
{
    double value = coefficients.back();
    for (size_t ii = 1; ii < coefficients.size(); ++ii)
    {
        value = x * value + coefficients[coefficients.size() - ii - 1];
    }
    return value;
}

int FindCorrespondenceWithTheoreticalIdx(const std::vector<Correspondence>& allEntries, size_t theoreticalIdxOfEntryToFind)
{
    for (int ii = 0; ii < static_cast<int>(allEntries.size()); ++ii)
    {
        if (allEntries[ii].theoreticalIdx == theoreticalIdxOfEntryToFind)
        {
            return ii;
        }
    }

    return -1; // not found
}

/// <summary>
/// Counts the number of inliers in the provided model.
/// </summary>
/// <param name="polynomialCoefficients"></param>
/// <param name="possibleCorrespondences"></param>
/// <param name="toleranceInWavelength"></param>
/// <param name="inlier">Will on return be filled with the found inliers. inlier.size() == return value</param>
/// <param name="averageError">Will on return be filled with an error estimate describing how well the inliers fit to the model</param>
/// <param name="isMonotonic">Will on return be set to 'true' if the polynomial is jugded to be monotonically increasing with pixel</param>
/// <returns>The number of inliers found</returns>
size_t CountInliers(
    const std::vector<double>& polynomialCoefficients,
    const std::vector<std::vector<Correspondence>>& possibleCorrespondences,
    double toleranceInWavelength,
    std::vector<Correspondence>& inlier,
    double& averageError,
    bool& isMonotonic)
{
    // Order the correspondences by the measured keypoint they belong to and only select one correspondence per keypoint
    averageError = 0.0;
    inlier.clear();
    inlier.reserve(100); // guess for the upper bound
    std::vector<double> distances;
    distances.reserve(100); // guess for the upper bound.
    std::vector<std::pair<double, double>> pixelToWavelengthMappings; // the mapping evaluated at the measured keypoints
    pixelToWavelengthMappings.resize(possibleCorrespondences.size());

    for (size_t measuredKeypointIdx = 0; measuredKeypointIdx < possibleCorrespondences.size(); ++measuredKeypointIdx)
    {
        if (possibleCorrespondences[measuredKeypointIdx].size() == 0)
        {
            continue;
        }

        // These correspondeces all have the same measured value, use that to only calculate the wavelength value once
        const double pixelValue = possibleCorrespondences[measuredKeypointIdx][0].measuredValue;
        const double predictedWavelength = CalculateValueAt(polynomialCoefficients, pixelValue);
        pixelToWavelengthMappings[measuredKeypointIdx].first = pixelValue;
        pixelToWavelengthMappings[measuredKeypointIdx].second = predictedWavelength;

        // Find the best possible correspondence for this measured keypoint.
        Correspondence bestcorrespondence;
        double smallestDistance = std::numeric_limits<double>::max();
        for (const Correspondence& corr : possibleCorrespondences[measuredKeypointIdx])
        {
            // Calculate the distance, in nm air, between the wavelength of this keypoint in the fraunhofer spectrum and the predicted wavelength of the measured keypoint
            const double wavelength = corr.theoreticalValue;
            const double distance = std::abs(predictedWavelength - wavelength); // in nm air

            if (distance < toleranceInWavelength && distance < smallestDistance)
            {
                bestcorrespondence = corr;
                smallestDistance = distance;
            }
        }

        if (smallestDistance < toleranceInWavelength)
        {
            // There are correspondences with the same theoreticalIdx in the incoming list
            //  but the selected inliers must be unique wrt both the measuredIdx and theoreticalIdx.
            //  Hence, check if we have already inserted a correspondence with this same theoreticalIdx in the result list.
            const int idxOfDuplicateEntry = FindCorrespondenceWithTheoreticalIdx(inlier, bestcorrespondence.theoreticalIdx);

            if (idxOfDuplicateEntry < 0)
            {
                // unique
                inlier.push_back(bestcorrespondence);
                distances.push_back(smallestDistance);
                averageError += smallestDistance;
            }
            else
            {
                // this theoreticalIdx already exists, compare the distances and keep the entry with the lowest distance
                if (smallestDistance < distances[idxOfDuplicateEntry])
                {
                    // Replace the existing entry with this
                    distances[idxOfDuplicateEntry] = smallestDistance;
                    inlier[idxOfDuplicateEntry] = bestcorrespondence;
                }
                // else, do nothing...
            }
        }
    }

    averageError /= (double)inlier.size();

    // Determine if the mapping is monotonically increasing
    std::sort(begin(pixelToWavelengthMappings), end(pixelToWavelengthMappings), [](const std::pair<double, double>& a, const std::pair<double, double>& b) { return a.first < b.first; });
    isMonotonic = true;
    for (size_t ii = 1; ii < pixelToWavelengthMappings.size(); ++ii)
    {
        if (pixelToWavelengthMappings[ii].second < pixelToWavelengthMappings[ii - 1].second)
        {
            isMonotonic = false;
            break;
        }
    }

    return inlier.size();
}

// ---------------------------------- RansacWavelengthCalibrationSetup ----------------------------------
RansacWavelengthCalibrationSetup::RansacWavelengthCalibrationSetup(RansacWavelengthCalibrationSettings calibrationSettings) :
    settings(calibrationSettings)
{
}

/// <summary>
/// Sorts the provided set of correspondences into a vector where the element at position 'N' lists all the correspondences for the measured keypoint 'N'
/// </summary>
std::vector<std::vector<Correspondence>> ArrangeByMeasuredKeypoint(const std::vector<Correspondence>& allCorrespondences)
{
    std::vector<std::vector<Correspondence>> result;
    result.reserve(100); // initial guess

    for (const Correspondence& c : allCorrespondences)
    {
        if (c.measuredIdx >= result.size())
        {
            result.resize(1 + c.measuredIdx);
        }
        result[c.measuredIdx].push_back(c);
    }

    return result;
}

std::vector<bool> ListInliers(const std::vector<Correspondence>& selectedValues, const std::vector<Correspondence>& allValues)
{
    std::vector<bool> result(allValues.size(), false);

    // TODO: This is not efficient!
    for (const auto& selectedValue : selectedValues)
    {
        for (size_t ii = 0; ii < allValues.size(); ++ii)
        {
            if (allValues[ii].measuredIdx == selectedValue.measuredIdx && allValues[ii].theoreticalIdx == selectedValue.theoreticalIdx)
            {
                result[ii] = true;
                break;
            }
        }
    }

    return result;
}


RansacWavelengthCalibrationResult RansacWavelengthCalibrationSetup::RunRansacCalibrations(
    const std::vector<Correspondence>& possibleCorrespondences,
    const std::vector<std::vector<Correspondence>>& possibleCorrespondencesOrderedByMeasuredKeypoint,
    int numberOfIterations) const
{
    std::random_device r;
    std::mt19937 rnd{ r() };
    PolynomialFit polyFit{ static_cast<int>(settings.modelPolynomialOrder) };

    RansacWavelengthCalibrationResult result(settings.modelPolynomialOrder);
    result.numberOfPossibleCorrelations = possibleCorrespondences.size();
    result.highestNumberOfInliers = 0U;
    result.smallestError = std::numeric_limits<double>::max();

    // variables in the loop, such that we don't have to recreate them every time
    const size_t ransacSampleSize = (settings.sampleSize > 0) ? settings.sampleSize : (settings.modelPolynomialOrder + 1);
    std::vector<double> selectedPixelValues(ransacSampleSize);
    std::vector<double> selectedWavelengths(ransacSampleSize);
    std::vector<double> suggestionForPolynomial;
    std::vector<Correspondence> selectedCorrespondences;

    for (int iteration = 0; iteration < numberOfIterations; ++iteration)
    {
        SelectMaybeInliers(ransacSampleSize, possibleCorrespondences, rnd, selectedCorrespondences);

        // Create a new (better?) model from these selected correspondences
        for (size_t ii = 0; ii < ransacSampleSize; ++ii)
        {
            selectedPixelValues[ii] = selectedCorrespondences[ii].measuredValue;
            selectedWavelengths[ii] = selectedCorrespondences[ii].theoreticalValue;
        }
        if (!polyFit.FitPolynomial(selectedPixelValues, selectedWavelengths, suggestionForPolynomial))
        {
            std::cout << "Polynomial fit failed" << std::endl;
            continue;
        }

        {
            const double firstValue = PolynomialValueAt(suggestionForPolynomial, 0);
            const double midpointValue = PolynomialValueAt(suggestionForPolynomial, settings.detectorSize * 0.5);
            const double lastValue = PolynomialValueAt(suggestionForPolynomial, (double)(settings.detectorSize - 1));

            if (firstValue > lastValue || firstValue > midpointValue || midpointValue > lastValue)
            {
                // This is not strictly a test for monotonically increasing function, just sampling. But it's fast and that is the main point here.
                // std::cout << "Found polynomial is not monotonically increasing, skipping" << std::endl;
                continue;
            }
        }

        // Evaluate if this suggested polynomial fits better than the guess we already have by counting how many of the 
        //  possible correspondences fits with the provided model.
        std::vector<Correspondence> inlierCorrespondences;
        double meanErrorOfModel = 0.0;
        bool isMonotonicallyIncreasing = true;
        size_t numberOfInliers = CountInliers(suggestionForPolynomial, possibleCorrespondencesOrderedByMeasuredKeypoint, settings.inlierLimitInWavelength, inlierCorrespondences, meanErrorOfModel, isMonotonicallyIncreasing);

        if (iteration == 0 ||
            (numberOfInliers > result.highestNumberOfInliers && isMonotonicallyIncreasing) ||
            (numberOfInliers == result.highestNumberOfInliers && isMonotonicallyIncreasing && meanErrorOfModel < result.smallestError))
        {
            if (settings.refine && numberOfInliers > settings.modelPolynomialOrder)
            {
                std::vector<double> pixelValues;
                std::vector<double> wavelengths;
                pixelValues.reserve(numberOfInliers);
                wavelengths.reserve(numberOfInliers);
                for (const auto& corr : inlierCorrespondences)
                {
                    pixelValues.push_back(corr.measuredValue);
                    wavelengths.push_back(corr.theoreticalValue);
                }

                if (!polyFit.FitPolynomial(pixelValues, wavelengths, result.bestFittingModelCoefficients))
                {
                    std::cout << "Polynomial fit failed" << std::endl;
                    continue;
                }

                // recount the inliers
                numberOfInliers = CountInliers(suggestionForPolynomial, possibleCorrespondencesOrderedByMeasuredKeypoint, settings.inlierLimitInWavelength, inlierCorrespondences, meanErrorOfModel, isMonotonicallyIncreasing);

                // std::cout << "Updating model at iteration " << iteration << ", new number of inliers is: " << numberOfInliers << ", mean error: " << meanErrorOfModel << std::endl;

                result.highestNumberOfInliers = numberOfInliers;
                result.correspondenceIsInlier = ListInliers(inlierCorrespondences, possibleCorrespondences);
                result.smallestError = meanErrorOfModel;
            }
            else
            {
                // std::cout << "Updating model at iteration " << iteration << ", new number of inliers is: " << numberOfInliers << std::endl;

                result.bestFittingModelCoefficients = std::move(suggestionForPolynomial);

                result.highestNumberOfInliers = numberOfInliers;
                result.correspondenceIsInlier = ListInliers(inlierCorrespondences, possibleCorrespondences);
                result.smallestError = meanErrorOfModel;
            }
        }
    }

    return result;
}

int GetBestResult(const std::vector< RansacWavelengthCalibrationResult>& partialResults)
{
    int bestIdx = 0;
    for (int ii = 1; ii < (int)partialResults.size(); ++ii)
    {
        if (partialResults[ii].highestNumberOfInliers > partialResults[bestIdx].highestNumberOfInliers ||
            (partialResults[ii].highestNumberOfInliers == partialResults[bestIdx].highestNumberOfInliers && partialResults[ii].smallestError < partialResults[bestIdx].smallestError))
        {
            bestIdx = ii;
        }
    }

    return bestIdx;
}


constexpr int OPENMP_DEFAULT_NUMBER_OF_THREADS = 4;
constexpr size_t OPENMP_MAX_NUMBER_OF_THREADS = 64LLU;

RansacWavelengthCalibrationResult RansacWavelengthCalibrationSetup::DoWavelengthCalibration(
    const std::vector<Correspondence>& possibleCorrespondences)
{
    const auto possibleCorrespondencesOrderedByMeasuredKeypoint = ArrangeByMeasuredKeypoint(possibleCorrespondences);
    std::vector< RansacWavelengthCalibrationResult> partialResults;

    const int numberOfThreads = (settings.numberOfThreads > 0) ? static_cast<int>(std::min(settings.numberOfThreads, OPENMP_MAX_NUMBER_OF_THREADS)) : OPENMP_DEFAULT_NUMBER_OF_THREADS;

    partialResults.resize(numberOfThreads);
    omp_set_num_threads(numberOfThreads);

#pragma omp parallel for
    for (int threadIdx = 0; threadIdx < numberOfThreads; ++threadIdx)
    {
        const int numberOfIterationsInThisThread = settings.numberOfRansacIterations / numberOfThreads;
        partialResults[threadIdx] = RunRansacCalibrations(possibleCorrespondences, possibleCorrespondencesOrderedByMeasuredKeypoint, numberOfIterationsInThisThread);
    }

    // Combine the results of the threads results into one final result
    const int bestResultIdx = GetBestResult(partialResults);
    RansacWavelengthCalibrationResult finalResult = partialResults[bestResultIdx];

    return finalResult;
}

}
