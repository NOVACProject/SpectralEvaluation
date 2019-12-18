#pragma once

// ---------------------------------------------------------------------------------------------------------------
// ----------------- This file contains misc. utility methods for simple numerical interpolation -----------------
// ---------------------------------------------------------------------------------------------------------------

#include <vector>

/** Finds the value of the provided vector at the given index by linear interpolation.
    @returns true if the interpolated could be found.
    @returns false if the given index is negative or larger than (values.size() - 1). */
bool LinearInterpolation(const std::vector<float>& values, double index, double& result);

/** Performs a tri-linear interpolation to find the value with the index {index[0], index[1], index[2]}
    in the three dimensional data matrix. The 'values' vector is a flattened three dimensional data
    matrix consisting of eight values, one value for each corner of the unit-cube. 
    This assumes that values is a cube with eight values and index is an index vector with three values
        which are all in the range [0, 1].
    @return true if the interpolation could be performed successfully.
    @return false if any of the above assumptions are not fulfilled. */
bool TriLinearInterpolation(const std::vector<double>& values, const std::vector<double>& index, double& result);
bool TriLinearInterpolation(const std::vector<double>& values, const std::vector<double>& index, double& result, double& estimatedVariation);

/** Finds the provided value in the given vector of values. 
    The vector of values is assumed to be sorted in either ascending or descending order.
    @return the index or NAN if the value could not be found. */
double GetFractionalIndex(const std::vector<float>& values, double valueToFind);
