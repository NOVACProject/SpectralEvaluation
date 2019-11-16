#include <SpectralEvaluation/Interpolation.h>
#include <cmath>

double LinearInterpolation(const std::vector<float>& values, double index)
{
    const size_t ii = (size_t)std::floor(index);

    if (ii >= values.size() - 2)
    {
        return NAN; // not found.
    }

    const double alpha = index - (double)ii;

    return values[ii] * (1.0 - alpha) + values[ii + 1] * alpha;
}

double GetFractionalIndex(const std::vector<float>& values, float valueToFind)
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
