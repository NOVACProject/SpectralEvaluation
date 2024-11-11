#include <SpectralEvaluation/VectorUtils.h>
#include <SpectralEvaluation/Math/PolynomialFit.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <assert.h>
#include <limits>

double Max(const std::vector<double>& values, size_t& idx)
{
    idx = 0;

    if (values.size() == 0)
    {
        return 0.0;
    }

    double m = values[0];
    for (size_t ii = 1; ii < values.size(); ++ii)
    {
        if (values[ii] > m)
        {
            m = values[ii];
            idx = ii;
        }
    }

    return m;
}

double Max(const std::vector<double>& values)
{
    size_t idx = 0U;
    return Max(values, idx);
}

double Max(std::vector<double>::const_iterator start, std::vector<double>::const_iterator end)
{
    if (start > end)
    {
        return 0.0;
    }
    if (start == end)
    {
        return *start;
    }

    double m = *start;
    for (auto it = start; it < end; ++it)
    {
        m = std::max(m, *it);
    }

    return m;
}


double MaxAbs(const std::vector<double>& values, size_t& idx)
{
    idx = 0;

    if (values.size() == 0)
    {
        return 0.0;
    }

    double m = std::abs(values[0]);
    for (size_t ii = 1; ii < values.size(); ++ii)
    {
        if (std::abs(values[ii]) > m)
        {
            m = std::abs(values[ii]);
            idx = ii;
        }
    }

    return m;
}

double MaxAbs(const std::vector<double>& values)
{
    size_t idx = 0U;
    return MaxAbs(values, idx);
}

double Min(const std::vector<double>& values, size_t& idx)
{
    idx = 0;

    if (values.size() == 0)
    {
        return 0.0;
    }

    double m = values[0];
    for (size_t ii = 1; ii < values.size(); ++ii)
    {
        if (values[ii] < m)
        {
            m = values[ii];
            idx = ii;
        }
    }

    return m;
}

double Min(const std::vector<double>& values)
{
    size_t idx = 0U;
    return Min(values, idx);
}


double Min(std::vector<double>::const_iterator start, std::vector<double>::const_iterator end)
{
    if (start > end)
    {
        return 0.0;
    }
    if (start == end)
    {
        return *start;
    }

    double m = *start;
    for (auto it = start; it < end; ++it)
    {
        m = std::min(m, *it);
    }

    return m;
}

std::pair<double, double> MinMax(const std::vector<double>& values)
{
    std::pair<size_t, size_t> idx;
    return MinMax(values, idx);
}

std::pair<double, double> MinMax(const std::vector<double>& values, std::pair<size_t, size_t>& idx)
{
    if (values.size() == 0)
    {
        idx.first = 0U;
        idx.second = 0U;
        return std::pair<double, double>(0.0, 0.0);
    }

    double minValue = values[0];
    double maxValue = values[1];
    size_t minIdx = 0U;
    size_t maxIdx = 0U;
    for (size_t ii = 1; ii < values.size(); ++ii)
    {
        if (values[ii] < minValue)
        {
            minValue = values[ii];
            minIdx = ii;
        }
        else if (values[ii] > maxValue)
        {
            maxValue = values[ii];
            maxIdx = ii;
        }
    }

    idx.first = minIdx;
    idx.second = maxIdx;
    return std::pair<double, double>(minValue, maxValue);
}

std::pair<double, double> MinMax(std::vector<double>::const_iterator start, std::vector<double>::const_iterator end, std::pair<size_t, size_t>& idx)
{
    if (start > end)
    {
        return std::pair<double, double>(0.0, 0.0);
    }
    if (start == end)
    {
        return std::pair<double, double>(*start, *start);
    }

    double minValue = *start;
    double maxValue = *start;
    size_t minIdx = 0U;
    size_t maxIdx = 0U;
    for (auto it = start; it < end; ++it)
    {
        if (*it < minValue)
        {
            minValue = *it;
            minIdx = it - start;
        }
        else if (*it > maxValue)
        {
            maxValue = *it;
            maxIdx = it - start;
        }
    }

    idx.first = minIdx;
    idx.second = maxIdx;
    return std::pair<double, double>(minValue, maxValue);
}

double Sum(std::vector<double>::const_iterator start, std::vector<double>::const_iterator end)
{
    double sum = 0.0;
    for (auto it = start; it < end; ++it)
    {
        sum += *it;
    }

    return sum;
}

double Sum(const std::vector<double>& values)
{
    return Sum(begin(values), end(values));
}

double SumAbs(const std::vector<double>& values)
{
    double sum = 0.0;
    for (size_t ii = 0; ii < values.size(); ++ii)
    {
        sum += std::abs(values[ii]);
    }

    return sum;
}

double SumOfSquaredDifferences(const std::vector<double>& a, const std::vector<double>& b)
{
    if (a.size() != b.size())
    {
        return -1.0;
    }

    double sum = 0.0;
    for (size_t ii = 0; ii < a.size(); ++ii)
    {
        const double diff = a[ii] - b[ii];
        sum += diff * diff;
    }

    return sum;
}

void Mult(std::vector<double>& values, double factor)
{
    for (double& v : values)
    {
        v *= factor;
    }
}

void Div(std::vector<double>& first, const std::vector<double>& second)
{
    if (first.size() != second.size())
    {
        throw std::invalid_argument("Invalid argument to 'Div'. The two vectors must have equal length.");
    }

    for (size_t i = 0; i < first.size(); i++)
    {
        if (second[i] != 0)
        {
            first[i] /= second[i];
        }
        else
        {
            first[i] = 0;
        }
    }
}

void Add(std::vector<double>& values, double factor)
{
    for (double& v : values)
    {
        v += factor;
    }
}

void Invert(std::vector<double>& values)
{
    for (double& v : values)
    {
        v = 1.0 / v;
    }
}

void Reverse(std::vector<double>& values)
{
    if (values.size() <= 1)
    {
        return;
    }

    const size_t midpoint = values.size() / 2;
    for (size_t ii = 0; ii < midpoint; ++ii)
    {
        const double tempValue = values[ii];
        values[ii] = values[values.size() - ii - 1];
        values[values.size() - ii - 1] = tempValue;
    }
}

void Mult(const std::vector<double>& firstVector, std::vector<double>& secondVectorAndResult)
{
    if (firstVector.size() != secondVectorAndResult.size())
    {
        throw std::invalid_argument("In multiplication, the first and the second vector must have equal length");
    }

    for (size_t ii = 0; ii < firstVector.size(); ++ii)
    {
        secondVectorAndResult[ii] = firstVector[ii] * secondVectorAndResult[ii];
    }
}

void Exp(std::vector<double>& values)
{
    for (double& v : values)
    {
        v = std::exp(v);
    }
}

void Log(std::vector<double>& values)
{
    for (double& v : values)
    {
        v = (v <= 0) ? 0.0 : std::log(v);
    }
}

void LowPassBinomial(std::vector<double>& values, int nIterations)
{
    if (values.size() < 3)
    {
        return; // nothing to do.
    }

    // The binomial filtering is done by, for each, calculating a new value
    // as (0.25 * value to the left + 0.5 * this value + 0.25 * value to the right)
    // In order to do this efficiently, we keep the three values (left, middle, right) in memory
    // and calculate the outputs based on this. Notice that this algo will in each iteration overwrite
    // the contents of the vector, hence it is important that we load the values of the vector into memory
    // before overwriting them.

    for (int iteration = 0; iteration < nIterations; ++iteration)
    {
        double leftValue = values[0];
        double thisValue = values[0];
        double rightValue = values[1];

        for (size_t ii = 0; ii < values.size() - 2; ++ii)
        {
            const double output = 0.5 * thisValue + 0.25 * (leftValue + rightValue);

            // update the three values before saving in the calculated value
            leftValue = values[ii];
            thisValue = values[ii + 1];
            rightValue = values[ii + 2];

            // save the calculated value
            values[ii] = output;
        }

        // The last two values are calculated separately (in order to not have to check the indices in the loop)
        values[values.size() - 2] = 0.5 * thisValue + 0.25 * (leftValue + rightValue);
        values[values.size() - 1] = 0.5 * rightValue + 0.25 * (thisValue + rightValue);
    }

    return;
}

void HighPassBinomial(std::vector<double>& values, int nIterations)
{
    // create copy of original data
    std::vector<double> fBuffer(begin(values), end(values));

    // create low pass filtered data
    LowPassBinomial(fBuffer, nIterations);

    // remove low pass part from data
    for (size_t i = 0; i < values.size(); i++)
    {
        if (fBuffer[i] != 0.0)
            values[i] /= fBuffer[i];
        else
            values[i] = 0;
    }

    return;
}

double Average(std::vector<double>::const_iterator start, std::vector<double>::const_iterator end)
{
    if (start > end)
    {
        return 0.0;
    }
    if (start == end)
    {
        return *start;
    }
    double sum = 0.0;
    double numberOfValues = 0.0;
    for (auto it = start; it < end; ++it)
    {
        sum += *it;
        numberOfValues += 1;
    }
    return sum / numberOfValues;
}

double Average(const std::vector<double>& values)
{
    return Average(begin(values), end(values));
}

double Variance(const std::vector<double>& values)
{
    if (values.size() <= 1)
    {
        return 0.0;
    }

    // First get the mean-value
    double mean = Average(values);

    double sum = 0;
    for (size_t k = 0; k < values.size(); ++k)
    {
        sum += (values[k] - mean) * (values[k] - mean);
    }
    sum = sum / static_cast<double>(values.size());
    return sum;
}

double Stdev(const std::vector<double>& values)
{
    if (values.size() <= 1)
    {
        return 0.0;
    }

    return std::sqrt(Variance(values));
}

double MinOfAbsolutes(const std::vector<double>& values)
{
    if (values.size() == 0)
    {
        return 0.0;
    }
    if (values.size() == 1)
    {
        return values[0];
    }

    double minValue = std::abs(values[0]);
    for (size_t ii = 1; ii < values.size(); ++ii)
    {
        minValue = std::min(minValue, std::abs(values[ii]));
    }
    return minValue;
}

double WeightedAverage(const std::vector<double>& values, const std::vector<double>& errors)
{
    if (values.size() != errors.size())
    {
        return 0.0;
    }
    if (values.size() == 0)
    {
        return 0.0;
    }
    if (values.size() == 1)
    {
        return values.front();
    }

    // in order to avoid some catastrophic cancellation here, extract the common order of magnitude for both the values and the errors and handle them separately
    const double errorsFactor = MinOfAbsolutes(errors);

    double sumOfValues = 0.0;
    double sumOfErrors = 0.0;
    for (size_t ii = 0; ii < values.size(); ++ii)
    {
        const double errorSquared = errors[ii] * errors[ii] / (errorsFactor * errorsFactor);

        sumOfValues += values[ii] / errorSquared;
        sumOfErrors += 1.0 / errorSquared;
    }

    return (sumOfValues / sumOfErrors);
}

void RemoveMean(std::vector<double>& values)
{
    double mean = Average(values);
    for (double& value : values)
    {
        value -= mean;
    }
}

void RemoveSlope(std::vector<double>& values)
{
    novac::PolynomialFit polyFit{ 1 };

    std::vector<double> indices;
    indices.resize(values.size());
    for (int ii = 0; ii < static_cast<int>(values.size()); ++ii)
    {
        indices[ii] = static_cast<double>(ii);
    }

    std::vector<double> polyomialCoefficients;
    polyFit.FitPolynomial(indices, values, polyomialCoefficients);

    for (int ii = 0; ii < static_cast<int>(values.size()); ++ii)
    {
        const double baseline = novac::PolynomialValueAt(polyomialCoefficients, static_cast<double>(ii));
        values[ii] -= baseline;
    }
}

double Median(std::vector<double>& values)
{
    if (values.size() == 0)
    {
        return 0.0;
    }
    else if (values.size() == 1)
    {
        return values[0];
    }
    std::sort(begin(values), end(values));

    const size_t midpointIndex = (values.size() >> 1);
    if (values.size() % 2 == 0)
    {
        return 0.5 * (values[midpointIndex] + values[midpointIndex + 1]);
    }
    else
    {
        return values[midpointIndex];
    }
}

double Area(const std::vector<double>& values, double xStep)
{
    if (values.size() < 2)
    {
        return 0.0;
    }

    double sum = 0.0;
    for (size_t ii = 1; ii < values.size(); ++ii)
    {
        sum += (values[ii] + values[ii - 1]);
    }

    sum *= 0.5 * xStep;

    return sum;
}


void FindNLowest(const std::vector<double>& input, size_t N, std::vector<double>& result)
{
    if (0 == N || 0 == input.size())
    {
        result.clear();
        return;
    }
    N = std::min(N, input.size()); // never sort beyond the end of the data

    std::vector<double> temp(begin(input), end(input));
    std::nth_element(begin(temp), begin(temp) + N, end(temp));
    result = std::vector<double>(begin(temp), begin(temp) + N);
    std::sort(begin(result), end(result));
}

void Normalize(std::vector<double>& values)
{
    const double minValue = Min(values);
    const double maxValue = Max(values);

    if (std::abs(minValue) > std::numeric_limits<double>::epsilon() || std::abs(maxValue - 1.0) > std::numeric_limits<double>::epsilon())
    {
        for (size_t ii = 0; ii < values.size(); ++ii)
        {
            values[ii] = (values[ii] - minValue) / (maxValue - minValue);
        }
    }
}

void Normalize(const std::vector<double>& input, std::vector<double>& output)
{
    output.resize(input.size());

    const double minValue = Min(input);
    const double maxValue = Max(input);

    for (size_t ii = 0; ii < input.size(); ++ii)
    {
        output[ii] = (input[ii] - minValue) / (maxValue - minValue);
    }
}

void NormalizeArea(const std::vector<double>& input, std::vector<double>& output)
{
    output.resize(input.size());

    if (input.size() == 0)
    {
        return;
    }

    const double minValue = Min(input);
    const double sumOfValues = Sum(input) - minValue * input.size();

    for (size_t ii = 0; ii < input.size(); ++ii)
    {
        output[ii] = (input[ii] - minValue) / sumOfValues;
    }

    assert(fabs(Sum(output) - 1.0) < 0.1);
}

void NormalizeArea(const std::vector<double>& input, double xStep, std::vector<double>& output)
{
    output.resize(input.size());

    const double minValue = Min(input);
    const double maxValue = Area(input, xStep);

    for (size_t ii = 0; ii < input.size(); ++ii)
    {
        output[ii] = (input[ii] - minValue) / (maxValue - minValue);
    }

    assert(fabs(Area(output, xStep) - 1.0) < 0.1);
}

double FindValue(const std::vector<double>& values, double valueToFind, size_t startIdx, size_t stopIdx)
{
    if (stopIdx <= startIdx || startIdx >= values.size())
    {
        return -1.0;
    }

    if (valueToFind == values[startIdx])
    {
        return (double)startIdx;
    }

    stopIdx = std::min(stopIdx, values.size() - 1);

    for (size_t idx = startIdx + 1; idx <= stopIdx; ++idx)
    {
        const double lastValue = values[idx - 1];
        const double thisValue = values[idx];

        if ((thisValue >= valueToFind && lastValue < valueToFind))
        {
            const double alpha = (thisValue - valueToFind) / (thisValue - lastValue);
            return (double)idx - alpha;
        }
        else if ((thisValue <= valueToFind && lastValue > valueToFind))
        {
            const double alpha = (lastValue - valueToFind) / (lastValue - thisValue);
            return (double)idx - 1 + alpha;
        }
    }

    // The vector doesn't contain the value to find.
    return -1.0;
}

double GetAt(const std::vector<double>& values, double idx)
{
    if (idx < 0.0)
    {
        return 0.0;
    }
    if (idx > (double)(values.size() - 1))
    {
        return 0.0;
    }

    // linear interpolation between floor(idx) and ceil(idx)
    double x1 = values.at((int)std::floor(idx));
    double x2 = values.at((int)std::ceil(idx));

    double alpha = idx - std::floor(idx);

    return x1 * (1 - alpha) + x2 * alpha;
}

double Centroid(const std::vector<double>& values)
{
    if (values.size() <= 1)
    {
        return 0.0;
    }

    double sumOfWeights = 0.0;
    double weightedSum = 0.0;
    for (size_t ii = 0; ii < values.size(); ++ii)
    {
        weightedSum += (double)(ii)*values[ii];
        sumOfWeights += values[ii];
    }
    const double centerOfMass = weightedSum / sumOfWeights;

    return centerOfMass;
}

double Resolution(const std::vector<double>& wavelGrid)
{
    const double minValue = wavelGrid.front();
    const double maxValue = wavelGrid.back();
    assert(maxValue > minValue);

    return (maxValue - minValue) / (double)(wavelGrid.size() - 1);
}

std::vector<double> GenerateVector(double minValue, double maxValue, size_t length)
{
    std::vector<double> result;
    result.resize(length);
    for (size_t ii = 0; ii < length; ++ii)
    {
        result[ii] = minValue + (maxValue - minValue) * static_cast<double>(ii) / static_cast<double>(length - 1);
    }
    return result;
}

bool Contains(const std::vector<size_t>& data, size_t value)
{
    for (size_t v : data)
    {
        if (v == value)
        {
            return true;
        }
    }
    return false;
}
