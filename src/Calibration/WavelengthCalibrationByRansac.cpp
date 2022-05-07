#include <SpectralEvaluation/Calibration/WavelengthCalibrationByRansac.h>

#include <algorithm>
#include <iostream>
#include <omp.h>
#include <memory>
#include <random>
#include <stdexcept>
#include <SpectralEvaluation/Calibration/Correspondence.h>
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
        largestPixelSpan(other.largestPixelSpan),
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
        this->largestPixelSpan = other.largestPixelSpan;
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

        // likewise, extract a small region around 'fraunhoferLocation' in the fraunhofer spectrum
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
        ::CBasicMath math;

        // to reduce the impact of the noise present, we do a low-pass filtering of the measured spectrum.
        // to reduce the impact from the slopes of the measured spectrum, we also do a high pass filtering
        std::unique_ptr<CSpectrum> filteredMeasuredSpectrum = std::make_unique<CSpectrum>(measuredSpectrum);
        math.HighPassBinomial(filteredMeasuredSpectrum->m_data, filteredMeasuredSpectrum->m_length, 500);
        math.LowPassBinomial(filteredMeasuredSpectrum->m_data, filteredMeasuredSpectrum->m_length, 5);

        // also perform a high pass filtering of the Fraunhofer spectrum, to make them compareable.
        std::unique_ptr<CSpectrum> filteredFraunhoferSpectrum = std::make_unique<CSpectrum>(fraunhoferSpectrum);
        math.HighPassBinomial(filteredFraunhoferSpectrum->m_data, filteredFraunhoferSpectrum->m_length, 500);

        std::vector<novac::Correspondence> possibleCorrespondences;
        possibleCorrespondences.reserve(3 * measuredKeypoints.size());
        for (size_t ii = 0; ii < measuredKeypoints.size(); ii++)
        {
            if (measuredKeypoints[ii].pixel >= correspondenceSettings.measuredPixelStart &&
                measuredKeypoints[ii].pixel <= correspondenceSettings.measuredPixelStop)
            {
                std::vector<novac::Correspondence> possibleCorrespondencesForThisMeasuredKeypoint;
                std::vector<double> errors; // keep track of the errors, we will use this later to filter the correspondences

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
                        corr.error = novac::MeasureCorrespondenceError(*filteredMeasuredSpectrum, measuredKeypoints[ii].pixel, *filteredFraunhoferSpectrum, fraunhoferKeypoints[jj].pixel, correspondenceSettings);
                        possibleCorrespondencesForThisMeasuredKeypoint.push_back(corr);

                        errors.push_back(corr.error);
                    }
                }

                // Filter out correspondences by, for each measured keypoint, keeping only the correspondence with the lowest error.
                if (possibleCorrespondencesForThisMeasuredKeypoint.size() > 0)
                {
                    if (possibleCorrespondencesForThisMeasuredKeypoint.size() > 1)
                    {
                        std::sort(
                            begin(possibleCorrespondencesForThisMeasuredKeypoint),
                            end(possibleCorrespondencesForThisMeasuredKeypoint),
                            [](const novac::Correspondence& c1, const novac::Correspondence& c2) { return c1.error < c2.error; });
                    }
                    possibleCorrespondences.push_back(possibleCorrespondencesForThisMeasuredKeypoint.front());
                }
            }
        }

        return possibleCorrespondences;
    }

    // ---------------------------------- Free functions used to run the calibration ----------------------------------


    bool FitPolynomial(novac::PolynomialFit& polyFit, const std::vector<Correspondence>& selectedCorrespondences, std::vector<double>& resultingPolynomial)
    {
        std::vector<double> selectedPixelValues(selectedCorrespondences.size());
        std::vector<double> selectedWavelengths(selectedCorrespondences.size());

        for (size_t ii = 0; ii < selectedCorrespondences.size(); ++ii)
        {
            selectedPixelValues[ii] = selectedCorrespondences[ii].measuredValue;
            selectedWavelengths[ii] = selectedCorrespondences[ii].theoreticalValue;
        }

        return polyFit.FitPolynomial(selectedPixelValues, selectedWavelengths, resultingPolynomial);
    }

    bool IsPossiblePixelToWavelengthCalibrationPolynomial(const std::vector<double>& candidatePolynomial, size_t numberOfPixels)
    {
        // Attempt to check if the polynomial is monotonically increasing by checking that the derivative is always positive
        //  this is done by checking that it is positive at the first value and does not change sign anywhere.
        std::vector<double> polynomialDerivative;
        polynomialDerivative.resize(candidatePolynomial.size() - 1);
        for (size_t orderIdx = 1; orderIdx < candidatePolynomial.size(); ++orderIdx)
        {
            polynomialDerivative[orderIdx - 1] = candidatePolynomial[orderIdx] * orderIdx;
        }

        if (PolynomialValueAt(polynomialDerivative, 0.0) < 0.0)
        {
            // Derivative is negative at the first point
            return false;
        }

        std::vector<std::complex<double>> polynomialRoots;
        if (FindRoots(polynomialDerivative, polynomialRoots))
        {
            bool hasRootInInterval = false;
            for (const auto& root : polynomialRoots)
            {
                if (root.real() > 0.0 && root.real() < numberOfPixels)
                {
                    hasRootInInterval = true;
                    break;
                }
            }

            if (hasRootInInterval)
            {
                return false;
            }
        }
        else
        {
            // root-checking failed, sample the polynomial instead
            const double firstValue = PolynomialValueAt(candidatePolynomial, 0);
            const double midpointValue = PolynomialValueAt(candidatePolynomial, numberOfPixels * 0.5);
            const double lastValue = PolynomialValueAt(candidatePolynomial, (double)(numberOfPixels - 1));

            if (firstValue > lastValue || firstValue > midpointValue || midpointValue > lastValue)
            {
                // This is not strictly a test for monotonically increasing function, just sampling. But it's fast and that is the main point here.
                // std::cout << "Found polynomial is not monotonically increasing, skipping" << std::endl;
                return false;
            }
        }

        // No problem found, seems legit.
        return true;
    }

    bool IsPossiblePixelToWavelengthCalibration(const std::vector<double>& pixelToWavelengthMapping)
    {
        if (pixelToWavelengthMapping.size() < 2)
        {
            return false;
        }
        for (size_t ii = 1; ii < pixelToWavelengthMapping.size(); ++ii)
        {
            if (pixelToWavelengthMapping[ii] < pixelToWavelengthMapping[ii - 1])
            {
                return false;
            }
        }

        return true;
    }


    // ---------------------------------- RansacWavelengthCalibrationSetup ----------------------------------
    RansacWavelengthCalibrationSetup::RansacWavelengthCalibrationSetup(RansacWavelengthCalibrationSettings calibrationSettings) :
        settings(calibrationSettings)
    {
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
        if (settings.initialModelCoefficients.size() == settings.modelPolynomialOrder + 1)
        {
            // We have an initial guess for the pixel-to-wavelength mapping. Measure how good it is.
            std::vector<Correspondence> inlierCorrespondences;
            bool isMonotonicallyIncreasing = true;
            result.highestNumberOfInliers = CountInliers(
                                                settings.initialModelCoefficients,
                                                possibleCorrespondencesOrderedByMeasuredKeypoint,
                                                settings.inlierLimitInWavelength,
                                                inlierCorrespondences,
                                                result.smallestError,
                                                isMonotonicallyIncreasing);

            if (!isMonotonicallyIncreasing)
            {
                result.highestNumberOfInliers = 0U;
                result.smallestError = std::numeric_limits<double>::max();
            }
        }
        else
        {
            result.highestNumberOfInliers = 0U;
            result.smallestError = std::numeric_limits<double>::max();
        }

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
            if (!FitPolynomial(polyFit, selectedCorrespondences, suggestionForPolynomial))
            {
                std::cout << "Polynomial fit failed" << std::endl;
                continue;
            }

            if (!IsPossiblePixelToWavelengthCalibrationPolynomial(suggestionForPolynomial, settings.detectorSize))
            {
                continue;
            }

            // Evaluate if this suggested polynomial fits better than the guess we already have by counting how many of the 
            //  possible correspondences fits with the provided model.
            std::vector<Correspondence> inlierCorrespondences;
            double meanErrorOfModel = 0.0;
            bool isMonotonicallyIncreasing = true;
            size_t numberOfInliers = CountInliers(suggestionForPolynomial, possibleCorrespondencesOrderedByMeasuredKeypoint, settings.inlierLimitInWavelength, inlierCorrespondences, meanErrorOfModel, isMonotonicallyIncreasing);
            const double inlierPixelSpan = GetMeasuredValueSpan(inlierCorrespondences);

            if ((numberOfInliers > result.highestNumberOfInliers && isMonotonicallyIncreasing) ||
                (numberOfInliers == result.highestNumberOfInliers && isMonotonicallyIncreasing && inlierPixelSpan > result.largestPixelSpan) ||
                (numberOfInliers == result.highestNumberOfInliers && isMonotonicallyIncreasing && std::abs(inlierPixelSpan - result.largestPixelSpan) < 0.1 && meanErrorOfModel < result.smallestError))
            {
                if (settings.refine && numberOfInliers > settings.modelPolynomialOrder + 1)
                {
                    if (!FitPolynomial(polyFit, inlierCorrespondences, suggestionForPolynomial))
                    {
                        std::cout << "Polynomial fit failed" << std::endl;
                        continue;
                    }

                    // recount the inliers
                    numberOfInliers = CountInliers(suggestionForPolynomial, possibleCorrespondencesOrderedByMeasuredKeypoint, settings.inlierLimitInWavelength, inlierCorrespondences, meanErrorOfModel, isMonotonicallyIncreasing);
                    result.largestPixelSpan = GetMeasuredValueSpan(inlierCorrespondences);
                }

                result.bestFittingModelCoefficients = std::move(suggestionForPolynomial);
                result.highestNumberOfInliers = numberOfInliers;
                result.correspondenceIsInlier = ListInliers(inlierCorrespondences, possibleCorrespondences);
                result.smallestError = meanErrorOfModel;
                result.largestPixelSpan = inlierPixelSpan;
            }
        }

        return result;
    }

    RansacWavelengthCalibrationResult RansacWavelengthCalibrationSetup::RunDeterministicCalibration(
        const std::vector<Correspondence>& possibleCorrespondences,
        const std::vector<std::vector<Correspondence>>& possibleCorrespondencesOrderedByMeasuredKeypoint) const
    {
        PolynomialFit polyFit{ static_cast<int>(settings.modelPolynomialOrder) };

        RansacWavelengthCalibrationResult result(settings.modelPolynomialOrder);
        result.numberOfPossibleCorrelations = possibleCorrespondences.size();
        result.highestNumberOfInliers = 0U;
        result.smallestError = std::numeric_limits<double>::max();

        // variables in the loop, such that we don't have to recreate them every time
        const size_t ransacSampleSize = (settings.sampleSize > 0) ? settings.sampleSize : (settings.modelPolynomialOrder + 1);
        std::vector<double> suggestionForPolynomial;
        std::vector<Correspondence> selectedCorrespondences;

        const auto allCorrespondenceCombinations = ListAllPossibleKeypointCombinations(ransacSampleSize, possibleCorrespondences);

        for (size_t iteration = 0; iteration < allCorrespondenceCombinations.size(); ++iteration)
        {
            selectedCorrespondences = allCorrespondenceCombinations[iteration];

            // TODO: Refactor and duplicated code in RunRansacCalibrations

            // Create a new (better?) model from these selected correspondences
            if (!FitPolynomial(polyFit, selectedCorrespondences, suggestionForPolynomial))
            {
                std::cout << "Polynomial fit failed" << std::endl;
                continue;
            }

            if (!IsPossiblePixelToWavelengthCalibrationPolynomial(suggestionForPolynomial, settings.detectorSize))
            {
                continue;
            }

            // Evaluate if this suggested polynomial fits better than the guess we already have by counting how many of the 
            //  possible correspondences fits with the provided model.
            std::vector<Correspondence> inlierCorrespondences;
            double meanErrorOfModel = 0.0;
            bool isMonotonicallyIncreasing = true;
            size_t numberOfInliers = CountInliers(suggestionForPolynomial, possibleCorrespondencesOrderedByMeasuredKeypoint, settings.inlierLimitInWavelength, inlierCorrespondences, meanErrorOfModel, isMonotonicallyIncreasing);
            double inlierPixelSpan = GetMeasuredValueSpan(inlierCorrespondences);

            if (iteration == 0 ||
                (numberOfInliers > result.highestNumberOfInliers && isMonotonicallyIncreasing) ||
                (numberOfInliers == result.highestNumberOfInliers && isMonotonicallyIncreasing && inlierPixelSpan > result.largestPixelSpan) ||
                (numberOfInliers == result.highestNumberOfInliers && isMonotonicallyIncreasing && std::abs(inlierPixelSpan - result.largestPixelSpan) < 0.1 && meanErrorOfModel < result.smallestError))
            {
                if (settings.refine && numberOfInliers > settings.modelPolynomialOrder + 1)
                {
                    if (!FitPolynomial(polyFit, inlierCorrespondences, suggestionForPolynomial))
                    {
                        std::cout << "Polynomial fit failed" << std::endl;
                        continue;
                    }

                    // recount the inliers
                    numberOfInliers = CountInliers(suggestionForPolynomial, possibleCorrespondencesOrderedByMeasuredKeypoint, settings.inlierLimitInWavelength, inlierCorrespondences, meanErrorOfModel, isMonotonicallyIncreasing);
                    inlierPixelSpan = GetMeasuredValueSpan(inlierCorrespondences);
                }

                result.bestFittingModelCoefficients = std::move(suggestionForPolynomial);
                result.highestNumberOfInliers = numberOfInliers;
                result.correspondenceIsInlier = ListInliers(inlierCorrespondences, possibleCorrespondences);
                result.smallestError = meanErrorOfModel;
                result.largestPixelSpan = inlierPixelSpan;
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

        if (possibleCorrespondences.size() < 10 * settings.modelPolynomialOrder)
        {
            // No use in random sampling, simply pick all the possible combinations.
            return RunDeterministicCalibration(possibleCorrespondences, possibleCorrespondencesOrderedByMeasuredKeypoint);
        }
        else
        {
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

}
