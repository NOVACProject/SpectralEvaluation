#include <SpectralEvaluation/VectorUtils.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <assert.h>

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

double Sum(const std::vector<double>& values)
{
    double sum = 0.0;
    for (size_t ii = 0; ii < values.size(); ++ii)
    {
        sum += values[ii];
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

void Invert(std::vector<double>& values)
{
    for (double& v : values)
    {
        v = 1.0 / v;
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

double Average(const std::vector<double>& values)
{
    if (values.size() == 0)
    {
        return 0.0;
    }
    if (values.size() == 1)
    {
        return values.front();
    }
    const double sum = Sum(values);
    return (sum / (double)values.size());
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

    stopIdx = std::min(stopIdx, values.size());

    for (size_t idx = startIdx + 1; idx < stopIdx; ++idx)
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

