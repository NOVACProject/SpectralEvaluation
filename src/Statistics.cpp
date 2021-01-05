#include <SpectralEvaluation/Statistics.h>
#include <cmath>

namespace novac
{
void IncrementalMeanAndVariance::Update(double newValue)
{
    ++numberOfValues;
    const double delta = newValue - mean;
    mean += delta / (double)numberOfValues;
    const double delta2 = newValue - mean;
    M2 += delta * delta2;
}

double IncrementalMeanAndVariance::Mean() const
{
    return mean;
}

double IncrementalMeanAndVariance::Variance() const
{
    if (numberOfValues < 2)
    {
        return 0.0;
    }
    return M2 / (double)(numberOfValues - 1);
}

double IncrementalMeanAndVariance::Stdev() const
{
    return std::sqrt(Variance());
}

std::uint64_t IncrementalMeanAndVariance::Count() const
{
    return numberOfValues;
}

}