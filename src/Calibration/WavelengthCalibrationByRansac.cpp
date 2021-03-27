#include <SpectralEvaluation/Calibration/WavelengthCalibrationByRansac.h>
#include <random>
#include <set>
#include <stdexcept>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Spectra/SpectrumUtils.h>
#include <SpectralEvaluation/VectorUtils.h>

#undef min
#undef max

// TODO: Try to remove these
#include <SpectralEvaluation/Fit/StandardFit.h>
#include <SpectralEvaluation/Fit/StandardMetricFunction.h>
#include <SpectralEvaluation/Fit/PolynomialFunction.h>
#include <SpectralEvaluation/Fit/CubicSplineFunction.h>


namespace novac
{

// ---------------------------------- PolynomialFit ----------------------------------

class CVectorDataFunction : public MathFit::IParamFunction
{
private:
    MathFit::CVector xVector;
    MathFit::CVector yVector;
public:
    CVectorDataFunction(const MathFit::CVector& x, const MathFit::CVector& y)
        : xVector(x), yVector(y)
    {
        if (xVector.GetSize() != yVector.GetSize())
        {
            throw std::invalid_argument("Cannot setup a vector data function with not equal length of the x- and y- vectors");
        }
    }

    virtual MathFit::TFitData GetLinearBasisFunction(MathFit::TFitData /*fXValue*/, int /*iParamID*/, bool /*bFixedID*/) override
    {
        throw std::domain_error("CVectorDataFunction::GetLinearBasisFunction is not implemented");
    }

    virtual MathFit::TFitData GetValue(MathFit::TFitData /*fXValue*/) override
    {
        throw std::domain_error("CVectorDataFunction::GetValue is not implemented");
    }

    virtual MathFit::CVector& GetValues(MathFit::CVector& vXValues, MathFit::CVector& vYTargetVector) override
    {
        if (vXValues.GetSize() != this->xVector.GetSize())
        {
            throw std::invalid_argument("The CVectorDataFunction is constructed to only use vectors of equal length.");
        }

        const int size = this->xVector.GetSize();
        for (int i = 0; i < size; i++)
        {
            vYTargetVector.SetAt(i, this->yVector.GetAt(i));
        }

        return vYTargetVector;
    }

    virtual MathFit::TFitData GetSlope(MathFit::TFitData /*fXValue*/) override
    {
        throw std::domain_error("CVectorDataFunction::GetSlope is not implemented");
    }
};

// The PolynomialFit is a helper class used to somewhat optimize the fitting of a polynomial to data again and again using the MathFit code from DOASIS
class PolynomialFit
{
private:
    MathFit::CPolynomialFunction functionToFit;

    const int polynomialOrder;

    // Copied from https://en.wikipedia.org/wiki/LU_decomposition
    /* INPUT: A - array of pointers to rows of a square matrix having dimension N
     *        Tol - small tolerance number to detect failure when the matrix is near degenerate
     * OUTPUT: Matrix A is changed, it contains a copy of both matrices L-E and U as A=(L-E)+U such that P*A=L*U.
     *        The permutation matrix is not stored as a matrix, but in an integer vector P of size N+1
     *        containing column indexes where the permutation matrix has "1". The last element P[N]=S+N,
     *        where S is the number of row exchanges needed for determinant computation, det(P)=(-1)^S
     */
    int LUPDecompose(double** A, int N, double Tol, int* P) {

        int i, j, k, imax;
        double maxA, * ptr, absA;

        for (i = 0; i <= N; i++)
            P[i] = i; //Unit permutation matrix, P[N] initialized with N

        for (i = 0; i < N; i++) {
            maxA = 0.0;
            imax = i;

            for (k = i; k < N; k++)
                if ((absA = fabs(A[k][i])) > maxA) {
                    maxA = absA;
                    imax = k;
                }

            if (maxA < Tol) return 0; //failure, matrix is degenerate

            if (imax != i) {
                //pivoting P
                j = P[i];
                P[i] = P[imax];
                P[imax] = j;

                //pivoting rows of A
                ptr = A[i];
                A[i] = A[imax];
                A[imax] = ptr;

                //counting pivots starting from N (for determinant)
                P[N]++;
            }

            for (j = i + 1; j < N; j++) {
                A[j][i] /= A[i][i];

                for (k = i + 1; k < N; k++)
                    A[j][k] -= A[j][i] * A[i][k];
            }
        }

        return 1;  //decomposition done 
    }

    /* INPUT: A,P filled in LUPDecompose; b - rhs vector; N - dimension
     * OUTPUT: x - solution vector of A*x=b
     */
    void LUPSolve(double** A, int* P, double* b, int N, double* x) {

        for (int i = 0; i < N; i++) {
            x[i] = b[P[i]];

            for (int k = 0; k < i; k++)
                x[i] -= A[i][k] * x[k];
        }

        for (int i = N - 1; i >= 0; i--) {
            for (int k = i + 1; k < N; k++)
                x[i] -= A[i][k] * x[k];

            x[i] /= A[i][i];
        }
    }

public:
    PolynomialFit(int order)
        : functionToFit(order), polynomialOrder(order)
    {
    }

    // Specialization of FitPolynomial which fits a cubic polynomial to four data points.
    bool FitCubicPolynomial(std::vector<double>& xData, std::vector<double>& yData, std::vector<double>& polynomialCoefficients)
    {
        polynomialCoefficients.resize(4);

        double A[16]{}; // the Vandermode Matrix
        for (size_t rowIdx = 0; rowIdx < 4; ++rowIdx)
        {
            A[rowIdx * 4] = 1.0;
            A[rowIdx * 4 + 1] = xData[rowIdx];
            A[rowIdx * 4 + 2] = std::pow(xData[rowIdx], 2.0);
            A[rowIdx * 4 + 3] = std::pow(xData[rowIdx], 3.0);
        }
        double B[4]{ yData[0],  yData[1],  yData[2],  yData[3] };
        // Multiply by transpose (TODO: this can be done manually in the setup above)
        double AtA1[4]{};
        double AtA2[4]{};
        double AtA3[4]{};
        double AtA4[4]{};
        double* AtA[4] {AtA1, AtA2, AtA3, AtA4}; // Transpose(A) * A
        double Atb[4]{}; // Transpose(A) * B
        for (size_t rowIdx = 0; rowIdx < 4; ++rowIdx)
        {
            Atb[rowIdx] =
                A[rowIdx + 0] * yData[0] +
                A[rowIdx + 4] * yData[1] +
                A[rowIdx + 8] * yData[2] +
                A[rowIdx + 12] * yData[3];

            for (size_t colIdx = 0; colIdx < 4; ++colIdx)
            {
                AtA[rowIdx][colIdx] =
                    A[rowIdx + 0] * A[colIdx + 0] +
                    A[rowIdx + 4] * A[colIdx + 4] +
                    A[rowIdx + 8] * A[colIdx + 8] +
                    A[rowIdx + 12] * A[colIdx + 12];
            }
        }

        // Solve the normal equations AtA * x = AtB
        int permutations[8]{};
        LUPDecompose(AtA, 4, 1e-9, permutations);
        LUPSolve(AtA, permutations, Atb, 4, polynomialCoefficients.data());

        return true;
    }


    bool FitPolynomial(std::vector<double>& xData, std::vector<double>& yData, std::vector<double>& polynomialCoefficients)
    {
        if (xData.size() != yData.size())
        {
            return false;
        }
        if (xData.size() == 4 && polynomialOrder == 3)
        {
            return FitCubicPolynomial(xData, yData, polynomialCoefficients);
        }

        const bool autoReleaseData = false;
        MathFit::CVector xVector{ &xData[0], (int)xData.size(), 1, autoReleaseData };
        MathFit::CVector yVector{ &yData[0], (int)yData.size(), 1, autoReleaseData };

        try
        {
            // This is a home-made replacement for the spline. In order to speed up the fitting somewhat.
            CVectorDataFunction cubicSplineRepresentation{ xVector, yVector };

            MathFit::CStandardMetricFunction diff(cubicSplineRepresentation, this->functionToFit);

            MathFit::CLeastSquareFit fit(diff);
            fit.SetFitRange(xVector);
            fit.SetMaxFitSteps(500);
            fit.PrepareMinimize();

            (void)fit.Minimize();
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

};

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

// Hard coded polynomial calculation for order = 3
inline double CalculateValueAt(const std::vector<double>& coefficients, double x)
{
    return coefficients[0] + x * (coefficients[1] + x * (coefficients[2] + x * coefficients[3]));
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
/// <param name="averageError"></param>
/// <returns></returns>
size_t CountInliers(
    const std::vector<double>& polynomialCoefficients,
    const std::vector<std::vector<Correspondence>>& possibleCorrespondences,
    double toleranceInWavelength,
    std::vector<Correspondence>& inlier,
    double& averageError)
{
    // Version 2. Order the correspondences by the measured keypoint they belong to and only select one correspondence per keypoint
    averageError = 0.0;
    inlier.clear();
    inlier.reserve(100); // guess for the upper bound
    std::vector<double> distances;
    distances.reserve(100); // guess for the upper bound.
    for (size_t keypointIdx = 0; keypointIdx < possibleCorrespondences.size(); ++keypointIdx)
    {
        if (possibleCorrespondences[keypointIdx].size() == 0)
        {
            continue;
        }

        // Find the best possible correspondence for this measured keypoint.
        Correspondence bestcorrespondence;
        double smallestDistance = std::numeric_limits<double>::max();

        const std::vector<Correspondence>& possibleCorrespondencesForThisKeypoint = possibleCorrespondences[keypointIdx];
        for (const Correspondence& corr : possibleCorrespondencesForThisKeypoint)
        {
            // Calculate the distance, in nm air, between the wavelength of this keypoint in the fraunhofer spectrum and the predicted wavelength of the measured keypoint
            const double pixelValue = corr.measuredValue;
            const double wavelength = corr.theoreticalValue;
            const double predictedWavelength = CalculateValueAt(polynomialCoefficients, pixelValue);
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

RansacWavelengthCalibrationResult RansacWavelengthCalibrationSetup::DoWavelengthCalibration(
    const std::vector<Correspondence>& possibleCorrespondences)
{
    RansacWavelengthCalibrationResult result(settings.modelPolynomialOrder);
    result.numberOfPossibleCorrelations = possibleCorrespondences.size();
    result.highestNumberOfInliers = 0U;
    result.smallestError = std::numeric_limits<double>::max();

    const auto possibleCorrespondencesOrderedByMeasuredKeypoint = ArrangeByMeasuredKeypoint(possibleCorrespondences);

    std::vector<Correspondence> selectedCorrespondences;
    PolynomialFit polyFit{ static_cast<int>(settings.modelPolynomialOrder) };

    // Seed with a real random value, if available
    std::random_device r;
    std::mt19937 rnd{ r() };

    for (int iteration = 0; iteration < settings.numberOfRansacIterations; ++iteration)
    {
        SelectMaybeInliers(settings.sampleSize, possibleCorrespondences, rnd, selectedCorrespondences);

        // Create a new (better?) model from these selected correspondences
        std::vector<double> selectedPixelValues(settings.sampleSize);
        std::vector<double> selectedWavelengths(settings.sampleSize);
        for (size_t ii = 0; ii < settings.sampleSize; ++ii)
        {
            selectedPixelValues[ii] = selectedCorrespondences[ii].measuredValue;
            selectedWavelengths[ii] = selectedCorrespondences[ii].theoreticalValue;
        }
        std::vector<double> suggestionForPolynomial;
        if (!polyFit.FitPolynomial(selectedPixelValues, selectedWavelengths, suggestionForPolynomial))
        {
            std::cout << "Polynomial fit failed" << std::endl;
            continue;
        }

        // Evaluate if this suggested polynomial fits better than the guess we already have by counting how many of the 
        //  possible correspondences fits with the provided model.
        std::vector<Correspondence> inlierCorrespondences;
        double meanErrorOfModel = 0.0;
        size_t numberOfInliers = CountInliers(suggestionForPolynomial, possibleCorrespondencesOrderedByMeasuredKeypoint, settings.inlierLimitInWavelength, inlierCorrespondences, meanErrorOfModel);
        if (iteration == 0 || numberOfInliers > result.highestNumberOfInliers || (numberOfInliers == result.highestNumberOfInliers && meanErrorOfModel < result.smallestError))
        {
            if (settings.refine && numberOfInliers > settings.modelPolynomialOrder)
            {
                std::vector<double> pixelValues;
                std::vector<double> wavelengths;
                pixelValues.reserve(numberOfInliers);
                wavelengths.reserve(numberOfInliers);
                for (const auto & corr : inlierCorrespondences)
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
                numberOfInliers = CountInliers(suggestionForPolynomial, possibleCorrespondencesOrderedByMeasuredKeypoint, settings.inlierLimitInWavelength, inlierCorrespondences, meanErrorOfModel);

                std::cout << "Updating model at iteration " << iteration << ", new number of inliers is: " << numberOfInliers << ", mean error: " << meanErrorOfModel << std::endl;

                result.highestNumberOfInliers = numberOfInliers;
                result.correspondenceIsInlier = ListInliers(inlierCorrespondences, possibleCorrespondences);
                result.smallestError = meanErrorOfModel;
            }
            else
            {
                std::cout << "Updating model at iteration " << iteration << ", new number of inliers is: " << numberOfInliers << std::endl;

                result.bestFittingModelCoefficients = std::move(suggestionForPolynomial);

                result.highestNumberOfInliers = numberOfInliers;
                result.correspondenceIsInlier = ListInliers(inlierCorrespondences, possibleCorrespondences);
                result.smallestError = meanErrorOfModel;
            }
        }
    }

    return result;
}

}
