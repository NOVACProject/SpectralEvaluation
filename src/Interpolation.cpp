#include <SpectralEvaluation/Interpolation.h>
#include <cmath>

namespace novac
{

template<class T>
bool LinearInterpolationT(const std::vector<T>& values, double index, double& result)
{
    const size_t ii = (size_t)std::floor(index);

    if (ii >= values.size() - 2)
    {
        return false; // not found.
    }

    const double alpha = index - (double)ii;
    result = values[ii] * (1.0 - alpha) + values[ii + 1] * alpha;
    return true;
}

template<class T>
bool LinearInterpolationT(const T* values, size_t length, double index, double& result)
{
    if (values == nullptr || length < 3)
    {
        return false; // no data
    }

    const size_t ii = (size_t)std::floor(index);

    if (ii >= length - 2)
    {
        return false; // not found.
    }

    const double alpha = index - (double)ii;
    result = values[ii] * (1.0 - alpha) + values[ii + 1] * alpha;
    return true;
}


bool LinearInterpolation(const std::vector<double>& values, double index, double& result)
{
    return LinearInterpolationT<double>(values, index, result);
}
bool LinearInterpolation(const std::vector<float>& values, double index, double& result)
{
    return LinearInterpolationT<float>(values, index, result);
}

bool LinearInterpolation(const float* values, size_t length, double index, double& result)
{
    return LinearInterpolationT<float>(values, length, index, result);
}
bool LinearInterpolation(const double* values, size_t length, double index, double& result)
{
    return LinearInterpolationT<double>(values, length, index, result);
}


double Stdev(const std::vector<double>& values, double average)
{
    if (values.size() <= 1)
    {
        return 0.0;
    }
    double diff = (values[0] - average) * (values[0] - average);
    for (size_t ii = 1; ii < values.size(); ++ii)
    {
        diff += (values[ii] - average) * (values[ii] - average);
    }
    return std::sqrt(diff / (values.size() - 1));
}

// Calculates the standard deviation of a vector of eight values.
double Stdev8(const std::vector<double>& values, double average)
{
    double diff0 = (values[0] - average) * (values[0] - average);
    double diff1 = (values[1] - average) * (values[1] - average);
    double diff2 = (values[2] - average) * (values[2] - average);
    double diff3 = (values[3] - average) * (values[3] - average);

    diff0 += (values[4] - average) * (values[4] - average);
    diff1 += (values[5] - average) * (values[5] - average);
    diff2 += (values[6] - average) * (values[6] - average);
    diff3 += (values[7] - average) * (values[7] - average);

    diff0 += diff1;
    diff2 += diff3;

    return std::sqrt((diff0 + diff2) / 7.0);
}

bool TriLinearInterpolation(const std::vector<double>& values, const std::vector<double>& index, double& result)
{
    double ignoredValue;
    return TriLinearInterpolation(values, index, result, ignoredValue);
}

bool TriLinearInterpolation(const std::vector<double>& values, const std::vector<double>& index, double& result, double& estimatedVariation)
{
    if (values.size() != 8 ||
        index.size() != 3 ||
        index[0] < 0.0 || index[0] > 1.0 ||
        index[1] < 0.0 || index[1] > 1.0 ||
        index[2] < 0.0 || index[2] > 1.0)
    {
        return false;
    }

    const double idxX = index[0];
    const double idxY = index[1];
    const double idxZ = index[2];

    const double c00 = values[0] * (1.0 - idxZ) + values[1] * idxZ;
    const double c01 = values[2] * (1.0 - idxZ) + values[3] * idxZ;
    const double c10 = values[4] * (1.0 - idxZ) + values[5] * idxZ;
    const double c11 = values[6] * (1.0 - idxZ) + values[7] * idxZ;

    const double c0 = c00 * (1.0 - idxY) + c01 * idxY;
    const double c1 = c10 * (1.0 - idxY) + c11 * idxY;

    result = c0 * (1.0 - idxX) + c1 * idxX;
    estimatedVariation = Stdev8(values, result);

    return true;
}

template<class T>
double GetFractionalIndex(const std::vector<T>& values, double valueToFind)
{
    for (size_t ii = 1; ii < values.size(); ++ii)
    {
        if (values[ii - 1] <= valueToFind && values[ii] >= valueToFind)
        {
            return (ii - 1) + (valueToFind - values[ii - 1]) / (values[ii] - values[ii - 1]);
        }
        else if (values[ii - 1] >= valueToFind && values[ii] <= valueToFind)
        {
            return (ii - 1) + 1 - (valueToFind - values[ii]) / (values[ii - 1] - values[ii]);
        }
    }

    // not found
    return NAN;
}

double GetFractionalIndex(const std::vector<float>& values, double valueToFind)
{
    return GetFractionalIndex<float>(values, valueToFind);
}

double GetFractionalIndex(const std::vector<double>& values, double valueToFind)
{
    return GetFractionalIndex<double>(values, valueToFind);
}

}
