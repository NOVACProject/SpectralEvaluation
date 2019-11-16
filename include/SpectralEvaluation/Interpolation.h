#pragma once

// ---------------------------------------------------------------------------------------------------------------
// ----------------- This file contains misc. utility methods for simple numerical interpolation -----------------
// ---------------------------------------------------------------------------------------------------------------

#include <vector>

/** Finds the value of the provided vector at the given index by linear interpolation.
    @returns the interpolated value, or NAN if the index is negative or larger than (values.size() - 1). */
double LinearInterpolation(const std::vector<float>& values, double index);

/** Finds the provided value in the given vector of values. 
    The vector of values is assumed to be sorted in either ascending or descending order.
    @return the index or NAN if the value could not be found. */
double GetFractionalIndex(const std::vector<float>& values, float valueToFind);
