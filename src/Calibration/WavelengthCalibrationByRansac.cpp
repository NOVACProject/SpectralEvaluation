#include <SpectralEvaluation/Calibration/WavelengthCalibrationByRansac.h>
#include <random>
#include <set>
#include <stdexcept>

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
    smallestError(other.smallestError),
    numberOfPossibleCorrelations(other.numberOfPossibleCorrelations)
{
}

RansacWavelengthCalibrationResult& RansacWavelengthCalibrationResult::operator=(const RansacWavelengthCalibrationResult& other)
{
    this->modelPolynomialOrder = other.modelPolynomialOrder;
    this->bestFittingModelCoefficients = std::vector<double>(begin(other.bestFittingModelCoefficients), end(other.bestFittingModelCoefficients));
    this->highestNumberOfInliers = other.highestNumberOfInliers;
    this->smallestError = other.smallestError;
    this->numberOfPossibleCorrelations = other.numberOfPossibleCorrelations;

    return *this;
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
    for ( ; it != coefficients.rend(); ++it)
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

                std::cout << "Updating model at iteration " << iteration << ", new number of inliers is: " << numberOfInliers << ", mean eror: " << meanErrorOfModel << std::endl;

                result.highestNumberOfInliers = numberOfInliers;
                result.smallestError = meanErrorOfModel;
                correspondenceIsInlier = isInlier;
            }
            else
            {
                std::cout << "Updating model at iteration " << iteration << ", new number of inliers is: " << numberOfInliers << std::endl;

                result.bestFittingModelCoefficients = std::move(suggestionForPolynomial);

                result.highestNumberOfInliers = numberOfInliers;
                result.smallestError = meanErrorOfModel;
                correspondenceIsInlier = isInlier;
            }
        }
    }

    // Save the result to the debug output file
    /*
    if (debugOutputFile.length() > 3)
    {
        rapidjson::Document doc;
        doc.SetObject();
        rapidjson::Document::AllocatorType& allocator = doc.GetAllocator();

        // The Keypoints in the Fraunhofer spectrum
        {
            rapidjson::Value pixelArray(rapidjson::kArrayType);
            rapidjson::Value wavelengthArray(rapidjson::kArrayType);
            rapidjson::Value intensityArray(rapidjson::kArrayType);
            rapidjson::Value typeArray(rapidjson::kArrayType);
            rapidjson::Value calcWavlengthArray(rapidjson::kArrayType);
            for (const SpectrumDataPoint& pt : fraunhoferKeypoints)
            {
                pixelArray.PushBack(pt.pixel, allocator);
                wavelengthArray.PushBack(pt.wavelength, allocator);
                intensityArray.PushBack(pt.intensity, allocator);
                typeArray.PushBack(pt.type, allocator);
                calcWavlengthArray.PushBack(result.bestFittingModel.GetValue(pt.pixel), allocator);
            }

            rapidjson::Value fraunhoferKeypointsDoc(rapidjson::kObjectType);
            fraunhoferKeypointsDoc.AddMember("description", "This is the set of peaks/valleys found in the convolved Fraunhofer spectrum", allocator);
            fraunhoferKeypointsDoc.AddMember("pixel", pixelArray, allocator);
            fraunhoferKeypointsDoc.AddMember("wavelength", wavelengthArray, allocator);
            fraunhoferKeypointsDoc.AddMember("intensity", intensityArray, allocator);
            fraunhoferKeypointsDoc.AddMember("type", typeArray, allocator);
            fraunhoferKeypointsDoc.AddMember("calcWavelength", calcWavlengthArray, allocator);

            doc.AddMember("fraunhoferKeypoints", fraunhoferKeypointsDoc, allocator);
        }

        // The Keypoints in the Measured spectrum
        {
            rapidjson::Value pixelArray(rapidjson::kArrayType);
            rapidjson::Value intensityArray(rapidjson::kArrayType);
            rapidjson::Value typeArray(rapidjson::kArrayType);
            rapidjson::Value calcWavlengthArray(rapidjson::kArrayType);
            for (const SpectrumDataPoint& pt : measuredKeypoints)
            {
                pixelArray.PushBack(pt.pixel, allocator);
                intensityArray.PushBack(pt.intensity, allocator);
                typeArray.PushBack(pt.type, allocator);
                calcWavlengthArray.PushBack(result.bestFittingModel.GetValue(pt.pixel), allocator);
            }

            rapidjson::Value measKeypointsDoc(rapidjson::kObjectType);
            measKeypointsDoc.AddMember("description", "This is the set of peaks/valleys found in the measured spectrum", allocator);
            measKeypointsDoc.AddMember("pixel", pixelArray, allocator);
            measKeypointsDoc.AddMember("intensity", intensityArray, allocator);
            measKeypointsDoc.AddMember("type", typeArray, allocator);
            measKeypointsDoc.AddMember("calcWavelength", calcWavlengthArray, allocator);

            doc.AddMember("measuredKeypoints", measKeypointsDoc, allocator);
        }

        // The Correspondences
        {
            rapidjson::Value pixelArray(rapidjson::kArrayType);
            rapidjson::Value wavelengthArray(rapidjson::kArrayType);
            rapidjson::Value calcWavlengthArray(rapidjson::kArrayType);
            rapidjson::Value measIntensityArray(rapidjson::kArrayType);
            rapidjson::Value fraunhoferIntensityArray(rapidjson::kArrayType);
            rapidjson::Value selectedArray(rapidjson::kArrayType);
            rapidjson::Value scoreArray(rapidjson::kArrayType);
            for (size_t ii = 0; ii < possibleCorrespondences.size(); ++ii)
            {
                const SpectrumDataPoint& pixelPt = measuredKeypoints[possibleCorrespondences[ii].measuredIdx];
                const SpectrumDataPoint& wavelPt = fraunhoferKeypoints[possibleCorrespondences[ii].theoreticalIdx];
                const bool isSelected = correspondenceIsInlier[ii];

                pixelArray.PushBack(pixelPt.pixel, allocator);
                wavelengthArray.PushBack(wavelPt.wavelength, allocator);
                calcWavlengthArray.PushBack(result.bestFittingModel.GetValue(pixelPt.pixel), allocator);
                measIntensityArray.PushBack(pixelPt.intensity, allocator);
                fraunhoferIntensityArray.PushBack(wavelPt.intensity, allocator);
                selectedArray.PushBack(isSelected, allocator);
                scoreArray.PushBack(possibleCorrespondences[ii].error, allocator);
            }

            rapidjson::Value correspondencesDoc(rapidjson::kObjectType);
            correspondencesDoc.AddMember("description", "This is the set of correspondences between the measured and the fraunhofer keypoints. Pixel is from the measured, wavelength is from the fraunhofer spectrum and calcWavelength is from the model", allocator);
            correspondencesDoc.AddMember("pixel", pixelArray, allocator);
            correspondencesDoc.AddMember("wavelength", wavelengthArray, allocator);
            correspondencesDoc.AddMember("calcWavelength", calcWavlengthArray, allocator);
            correspondencesDoc.AddMember("measIntensity", measIntensityArray, allocator);
            correspondencesDoc.AddMember("fraunhoferIntensity", fraunhoferIntensityArray, allocator);
            correspondencesDoc.AddMember("selected", selectedArray, allocator);
            correspondencesDoc.AddMember("score", scoreArray, allocator);

            doc.AddMember("correspondences", correspondencesDoc, allocator);
        }

        // The initial calibration
        if (initialWavelengthCalibrationForDebug != nullptr)
        {
            // initialWavelengthCalibrationForDebug
            rapidjson::Value calibrationArray(rapidjson::kArrayType);
            for (double value : *initialWavelengthCalibrationForDebug)
            {
                calibrationArray.PushBack(value, allocator);
            }
            doc.AddMember("initialCalibration", calibrationArray, allocator);
        }

        // The final model
        {
            rapidjson::Value polynomial;
            polynomial.SetObject();

            polynomial.AddMember("order", (int)result.modelPolynomialOrder, allocator);

            for (int orderIdx = 0; orderIdx <= (int)result.modelPolynomialOrder; ++orderIdx)
            {
                char name[64];
                sprintf(name, "c%d", orderIdx);

                rapidjson::Value coefficientName(name, allocator);

                polynomial.AddMember(coefficientName, result.bestFittingModel.GetCoefficient(orderIdx), allocator);
            }
            doc.AddMember("polynomial", polynomial, allocator);
        }

        rapidjson::StringBuffer sb;
        rapidjson::PrettyWriter<rapidjson::StringBuffer> writer(sb);
        doc.Accept(writer);    // Accept() traverses the DOM and generates Handler events.

        std::ofstream dst{ debugOutputFile };
        dst << sb.GetString() << newLine;
    }

    */

    return result;
}

}
