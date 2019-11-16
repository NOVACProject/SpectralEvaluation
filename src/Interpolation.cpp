#include <SpectralEvaluation/Interpolation.h>
#include <cmath>

bool LinearInterpolation(const std::vector<float>& values, double index, double& result)
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
    estimatedVariation = c1 - c0; // TODO: this isn't really correct is it?

    return true;
}

double GetFractionalIndex(const std::vector<float>& values, double valueToFind)
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
    return NAN; // not found.
}
