#include "VectorUtils.h"
#include <algorithm>
#include <cmath>

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

double Sum(const std::vector<double>& values)
{
    double sum = 0.0;
    for(size_t ii = 0; ii < values.size(); ++ii)
    {
        sum += values[ii];
    }
    
    return sum;
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

    const double minValue = Min(input);
    const double maxValue = Sum(input);

    for (size_t ii = 0; ii < input.size(); ++ii)
    {
        output[ii] = (input[ii] - minValue) / (maxValue - minValue);
    }
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
    if (idx >(double)(values.size() - 1))
    {
        return 0.0;
    }

    // linear interpolation between floor(idx) and ceil(idx)
    double x1 = values.at((int)std::floor(idx));
    double x2 = values.at((int)std::ceil(idx));

    double alpha = idx - std::floor(idx);

    return x1 * (1 - alpha) + x2 * alpha;
}

