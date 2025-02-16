#pragma once

// ---------------------------------------------------------------------------------------------------------------
// ----------------- This file contains misc. utility methods for simple numerical interpolation -----------------
// ---------------------------------------------------------------------------------------------------------------

#include <vector>
#include <stddef.h>

namespace novac
{

/** Finds the value of the provided vector at the given index by linear interpolation.
    @returns true if the interpolated could be found.
    @returns false if the given index is negative or larger than (values.size() - 1). */
bool LinearInterpolation(const std::vector<float>& values, double index, double& result);
bool LinearInterpolation(const std::vector<double>& values, double index, double& result);

/** Finds the value of the provided data array at the given index by linear interpolation.
    @returns true if the interpolated could be found.
    @returns false if the given index is negative or larger than (length - 1). */
bool LinearInterpolation(const float* values, size_t length, double index, double& result);
bool LinearInterpolation(const double* values, size_t length, double index, double& result);

/** Performs a tri-linear interpolation to find the value with the index {index[0], index[1], index[2]}
    in the three dimensional data matrix. The 'values' vector is a flattened three dimensional data
    matrix consisting of eight values, one value for each corner of the unit-cube.
    This assumes that values is a cube with eight values and index is an index vector with three values
        which are all in the range [0, 1].
    @return true if the interpolation could be performed successfully.
    @return false if any of the above assumptions are not fulfilled. */
bool TriLinearInterpolation(const std::vector<double>& values, const std::vector<double>& index, double& result);
bool TriLinearInterpolation(const std::vector<double>& values, const std::vector<double>& index, double& result, double& estimatedVariation);

/** Locates the provided value in the given vector of values.
    The vector of values is assumed to be sorted in either ascending or descending order,
        the result is undefined if the vector is not sorted.
    @return the index or NAN if the value could not be found. */
double GetFractionalIndex(const std::vector<float>& values, double valueToFind);

/** Finds the provided value in the given vector of values.
    The vector of values is assumed to be sorted in either ascending or descending order.
    @return the index or NAN if the value could not be found. */
double GetFractionalIndex(const std::vector<double>& values, double valueToFind);

/** Performs a resampling of the data set provided by (x, y) onto the new x-axis grid (newX).
    The result will be returned in the provided 'result' vector which will be resampled to same size as 'newX'.
    If 'newX' contains any value outsize of the range (x[0], x[length-1]) then 'result' will be zero at those values.
    @throws std::invalid_argument if the length of x does not equal the length of y.
    @throws std::invalid_argument if x is not monotonically increasing. */
void Resample(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& newX, std::vector<double>& result);

}
