#include <SpectralEvaluation/Calibration/WavelengthCalibrationByRansac.h>
#include <random>
#include <set>
#include <stdexcept>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Spectra/SpectrumUtils.h>
#include <SpectralEvaluation/VectorUtils.h>

// TODO: Try to remove these
#include <SpectralEvaluation/Fit/StandardFit.h>
#include <SpectralEvaluation/Fit/StandardMetricFunction.h>
#include <SpectralEvaluation/Fit/PolynomialFunction.h>
#include <SpectralEvaluation/Fit/CubicSplineFunction.h>


namespace novac
{

// ---------------------------------- RansacWavelengthCalibrationResult ----------------------------------
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
    const ::CSpectrum& measuredSpectrum,
    double pixelInMeasuredSpectrum,
    const ::CSpectrum& fraunhoferSpectrum,
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
    const RansacWavelengthCalibrationSettings& ransacSettings,
    const CorrespondenceSelectionSettings& correspondenceSettings)
{
    std::vector<novac::Correspondence> possibleCorrespondences;
    std::vector<double> errors; // keep track of the errors, we will use this later to filter the correspondences
    for (size_t ii = 0; ii < measuredKeypoints.size(); ii++)
    {
        if (measuredKeypoints[ii].pixel >= correspondenceSettings.measuredPixelStart && measuredKeypoints[ii].pixel <= correspondenceSettings.measuredPixelStop)
        {
            for (size_t jj = 0; jj < fraunhoferKeypoints.size(); jj++)
            {
                if (measuredKeypoints[ii].type == fraunhoferKeypoints[jj].type &&
                    std::abs(measuredKeypoints[ii].pixel - fraunhoferKeypoints[jj].pixel) <= ransacSettings.maximumPixelDistanceForPossibleCorrespondence)
                {
                    novac::Correspondence corr;
                    corr.measuredIdx = ii;
                    corr.measuredValue = measuredKeypoints[ii].pixel;
                    corr.theoreticalIdx = jj;
                    corr.theoreticalValue = fraunhoferKeypoints[jj].wavelength;
                    corr.error = novac::MeasureCorrespondenceError(measuredSpectrum, measuredKeypoints[ii].pixel, fraunhoferSpectrum, fraunhoferKeypoints[jj].pixel, correspondenceSettings);
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
std::vector<Correspondence> SelectMaybeInliers(size_t number, const std::vector<Correspondence>& allCorrespondences, std::mt19937& randomGenerator)
{
    if (number == allCorrespondences.size())
    {
        return allCorrespondences;
    }
    else if (number > allCorrespondences.size())
    {
        throw std::invalid_argument("Cannot select a larger number of inliers than the entire data set.");
    }

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

    std::vector<Correspondence> result;
    result.reserve(number);
    for (size_t idx : selectedIndices)
    {
        result.push_back(allCorrespondences[idx]);
    }

    return result;
}

// TODO: MOVE AND MERGE WITH FITFUNCTION IN INSTRUMENTLINESHAPE.CPP
bool FitPolynomial(MathFit::CVector& xData, MathFit::CVector& yData, size_t order, std::vector<double>& polynomialCoefficients)
{
    try
    {
        MathFit::CCubicSplineFunction cubicSplineRepresentation{ xData, yData };

        MathFit::CPolynomialFunction functionToFit{ static_cast<int>(order) };
        MathFit::CStandardMetricFunction diff(cubicSplineRepresentation, functionToFit);

        MathFit::CStandardFit fit(diff);
        fit.SetFitRange(xData);
        fit.SetMaxFitSteps(500);
        fit.PrepareMinimize();

        if (!fit.Minimize())
        {
            return false;
        }
        fit.FinishMinimize();

        // Now extract the coefficients
        {
            const auto& modelVector = functionToFit.GetCoefficients();
            polynomialCoefficients.resize(modelVector.GetSize());
            for (int orderIdx = 0; orderIdx < modelVector.GetSize(); ++orderIdx)
            {
                polynomialCoefficients[orderIdx] = modelVector.GetAt(orderIdx);
            }
        }

        return true;
    }
    catch (MathFit::CFitException& e)
    {
        std::cout << "Fit failed: " << e.mMessage << std::endl;

        return false;
    }
}
bool FitPolynomial(std::vector<double>& xData, std::vector<double>& yData, size_t order, std::vector<double>& polynomialCoefficients)
{
    if (xData.size() != yData.size())
    {
        return false;
    }

    const bool autoReleaseData = false;
    MathFit::CVector xVector{ &xData[0], (int)xData.size(), 1, autoReleaseData };
    MathFit::CVector yVector{ &yData[0], (int)yData.size(), 1, autoReleaseData };

    return FitPolynomial(xVector, yVector, order, polynomialCoefficients);
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

// Hard coded polynomial calculation for order = 3
inline double CalculateValueAt(const std::vector<double>& coefficients, double x)
{
    return coefficients[0] + x * (coefficients[1] + x * (coefficients[2] + x * coefficients[3]));
}

size_t CountInliers(const std::vector<double>& polynomialCoefficients, const std::vector<Correspondence> possibleCorrespondences, double toleranceInWavelength, std::vector<bool>& isInlier, double& averageError)
{
    size_t count = 0;
    isInlier.resize(possibleCorrespondences.size());
    std::fill(begin(isInlier), end(isInlier), false); // set all values to false

    struct CorrespondenceWithDistance
    {
        Correspondence correspondence;
        double distance = 0.0;
        size_t index = 0U;
    };

    // First make a prelimiary list of all correspondences which fulfill the distance criterion.
    size_t maxMeasuredIdx = 0U; // This is used in the loop below...
    size_t maxTheoreticalIdx = 0U; // This is used in the loop below...
    std::vector<CorrespondenceWithDistance> selectedCorrespondences;
    selectedCorrespondences.reserve(possibleCorrespondences.size() / 10); // guess for the final number of correspondences.
    for (size_t ii = 0; ii < possibleCorrespondences.size(); ++ii)
    {
        maxMeasuredIdx = std::max(maxMeasuredIdx, possibleCorrespondences[ii].measuredIdx);
        maxTheoreticalIdx = std::max(maxTheoreticalIdx, possibleCorrespondences[ii].theoreticalIdx);

        const double pixelValue = possibleCorrespondences[ii].measuredValue;
        const double wavelength = possibleCorrespondences[ii].theoreticalValue;
        const double predictedWavelength = CalculateValueAt(polynomialCoefficients, pixelValue); // TODO: Generalize
        const double distance = std::abs(predictedWavelength - wavelength);

        if (distance < toleranceInWavelength)
        {
            CorrespondenceWithDistance result;
            result.correspondence = possibleCorrespondences[ii];
            result.distance = distance;
            result.index = ii;
            selectedCorrespondences.push_back(result);
        }
    }

    // Now we have a preliminary list, which may contain duplicates in both measured and theoretical index.
    //  Sort them by increasing distance, and select the ones with the smallest distances first.
    std::sort(begin(selectedCorrespondences), end(selectedCorrespondences), [](const CorrespondenceWithDistance& c1, const CorrespondenceWithDistance& c2) { return c1.distance < c2.distance; });
    std::vector<bool> selectedMeasuredIdx(maxMeasuredIdx + 1, false);
    std::vector<bool> selectedTheoreticalIdx(maxTheoreticalIdx + 1, false);
    averageError = 0.0;
    for each (const CorrespondenceWithDistance & corr in selectedCorrespondences)
    {
        if (selectedMeasuredIdx[corr.correspondence.measuredIdx] == false &&
            selectedTheoreticalIdx[corr.correspondence.theoreticalIdx] == false)
        {
            // These points have not been handled already. Insert.
            isInlier[corr.index] = true;
            ++count;
            averageError += corr.distance;
            selectedMeasuredIdx[corr.correspondence.measuredIdx] = true;
            selectedTheoreticalIdx[corr.correspondence.theoreticalIdx] = true;
        }
    }

    averageError /= (double)count;

    return count;
}

// ---------------------------------- RansacWavelengthCalibrationSetup ----------------------------------
RansacWavelengthCalibrationSetup::RansacWavelengthCalibrationSetup(RansacWavelengthCalibrationSettings calibrationSettings) :
    settings(calibrationSettings)
{
}

RansacWavelengthCalibrationResult RansacWavelengthCalibrationSetup::DoWavelengthCalibration(
    const std::vector<Correspondence>& possibleCorrespondences)
{
    RansacWavelengthCalibrationResult result(settings.modelPolynomialOrder);
    result.numberOfPossibleCorrelations = possibleCorrespondences.size();
    result.highestNumberOfInliers = 0U;
    result.smallestError = std::numeric_limits<double>::max();
    std::vector<bool> correspondenceIsInlier;

    // Seed with a real random value, if available
    std::random_device r;
    std::mt19937 rnd{ r() };

    for (int iteration = 0; iteration < settings.numberOfRansacIterations; ++iteration)
    {
        std::vector<Correspondence> selectedCorrespondences = SelectMaybeInliers(settings.sampleSize, possibleCorrespondences, rnd);

        // Create a new (better?) model from these selected correspondences
        std::vector<double> selectedPixelValues(settings.sampleSize);
        std::vector<double> selectedWavelengths(settings.sampleSize);
        for (size_t ii = 0; ii < settings.sampleSize; ++ii)
        {
            selectedPixelValues[ii] = selectedCorrespondences[ii].measuredValue;
            selectedWavelengths[ii] = selectedCorrespondences[ii].theoreticalValue;
        }
        std::vector<double> suggestionForPolynomial;
        if (!FitPolynomial(selectedPixelValues, selectedWavelengths, settings.modelPolynomialOrder, suggestionForPolynomial))
        {
            std::cout << "Polynomial fit failed" << std::endl;
            continue;
        }

        // Evaluate if this suggested polynomial fits better than the guess we already have by counting how many of the 
        //  possible correspondences fits with the provided model.
        std::vector<bool> isInlier;
        double meanErrorOfModel = 0.0;
        size_t numberOfInliers = CountInliers(suggestionForPolynomial, possibleCorrespondences, settings.inlierLimitInWavelength, isInlier, meanErrorOfModel);
        if (iteration == 0 || numberOfInliers > result.highestNumberOfInliers || (numberOfInliers == result.highestNumberOfInliers && meanErrorOfModel < result.smallestError))
        {
            if (settings.refine && numberOfInliers > settings.modelPolynomialOrder)
            {
                std::vector<double> pixelValues;
                std::vector<double> wavelengths;
                pixelValues.reserve(numberOfInliers);
                wavelengths.reserve(numberOfInliers);
                for (size_t ii = 0; ii < possibleCorrespondences.size(); ++ii)
                {
                    if (isInlier[ii])
                    {
                        pixelValues.push_back(possibleCorrespondences[ii].measuredValue);
                        wavelengths.push_back(possibleCorrespondences[ii].theoreticalValue);
                    }
                }
                if (!FitPolynomial(pixelValues, wavelengths, settings.modelPolynomialOrder, result.bestFittingModelCoefficients))
                {
                    std::cout << "Polynomial fit failed" << std::endl;
                    continue;
                }

                // recount the inliers
                numberOfInliers = CountInliers(suggestionForPolynomial, possibleCorrespondences, settings.inlierLimitInWavelength, isInlier, meanErrorOfModel);

                std::cout << "Updating model at iteration " << iteration << ", new number of inliers is: " << numberOfInliers << ", mean error: " << meanErrorOfModel << std::endl;

                result.highestNumberOfInliers = numberOfInliers;
                result.correspondenceIsInlier = isInlier;
                result.smallestError = meanErrorOfModel;
                correspondenceIsInlier = isInlier;
            }
            else
            {
                std::cout << "Updating model at iteration " << iteration << ", new number of inliers is: " << numberOfInliers << std::endl;

                result.bestFittingModelCoefficients = std::move(suggestionForPolynomial);

                result.highestNumberOfInliers = numberOfInliers;
                result.correspondenceIsInlier = isInlier;
                result.smallestError = meanErrorOfModel;
                correspondenceIsInlier = isInlier;
            }
        }
    }

    return result;
}

}
