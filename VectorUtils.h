#pragma once

// ---------------------------------------------------------------------------------------------------------------
// --------------------- This file contains misc. utility methods for working on std::vector ---------------------
// ---------------------------------------------------------------------------------------------------------------

#include <vector>

/** @return the maximum value in the provided vector. 
    @param idx will be filled with the index where the maximum value is found.
    If values.size() == 0 then 0.0 is returned and idx is set to zero. */
double Max(const std::vector<double>& values, size_t& idx);

/** @return the maximum value in the provided vector.
    If values.size() == 0 then 0.0 is returned. */
double Max(const std::vector<double>& values);

/** @return the minimum value in the provided vector.
    @param idx will be filled with the index where the minimum value is found.
    If values.size() == 0 then 0.0 is returned and idx is set to zero. */
double Min(const std::vector<double>& values, size_t& idx);

/** @return the minimum value in the provided vector.
    If values.size() == 0 then 0.0 is returned. */
double Min(const std::vector<double>& values);

/** @return the sum of all the given values 
    If values.size() == 0 then 0.0 is returned. */
double Sum(const std::vector<double>& values);

/** @return the Average of all the given values
    If values.size() == 0 then 0.0 is returned. */
double Average(const std::vector<double>& values);

/** Finds the 'N' lowest values in the input vector and fills them into the result vector
    If N == 0 or input.size() == 0 then result will be an empty vector upon return. */
void FindNLowest(const std::vector<double>& input, size_t N, std::vector<double>& result);

/** Normalizes a vector of values such that the highest value will be 1.0 and the lowest 0.0.
    If input.size() == 0 then output.size() will also be zero. */
void Normalize(const std::vector<double>& input, std::vector<double>& output);

/** Normalizes a vector of values such that the lowest value will be 0.0 and the
    sum of all the valeus will be 1.0.
    If input.size() == 0 then output.size() will also be zero. */
void NormalizeArea(const std::vector<double>& input, std::vector<double>& output);

/** Finds the (fractiona) index where the values in the provided vector crosses the y='valueToFind' line 
    in the index-range [startIdx, stopIdx[. 
     Returns -1.0 if the value cannot be found OR startIdx >= stopIdx OR startIdx >= values.size() */
double FindValue(const std::vector<double>& values, double valueToFind, size_t startIdx, size_t stopIdx);

/** Returns the value of the given vector at the fractional index. 
    The result will be linearly interpolated between the lower and the upper index.
    @return 0.0 if idx < 0.0 or idx > values.size() - 1 */
double GetAt(const std::vector<double>& values, double idx);

