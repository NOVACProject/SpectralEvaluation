#pragma once

#include <cstdint>

// ---------------------------------------------------------------------------------------------------------------
// ----------------- This file contains misc. utility methods for simple statistics calculations -----------------
// ---------------------------------------------------------------------------------------------------------------

namespace novac
{

/** Helper class for calculating mean and variance using Welford's online algorithm. 
    See e.g. https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance#Welford's_online_algorithm */
class IncrementalMeanAndVariance
{
public:

    /** Adds one more value to this statistics */
    void Update(double newValue);

    /** Returns the mean value calculated up until this point. */
    double Mean() const;

    /** Returns the sample variance calculated up until this point.
        This will return zero if the number of samples is less than two. */
    double Variance() const;

    /** Returns the sample standard deviation calculated up until this point.
        This is the square-root of the variance */
    double Stdev() const;

    /** The number of samples updated */
    std::uint64_t Count() const;

private:
    /** The number of values added so far. */
    std::uint64_t numberOfValues = 0;

    /** The current average value */
    double mean = 0.0;

    /** The second moment */
    double M2 = 0.0;
};
}