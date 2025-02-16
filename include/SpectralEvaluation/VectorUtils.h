#pragma once

// ---------------------------------------------------------------------------------------------------------------
// --------------------- This file contains misc. utility methods for working on std::vector ---------------------
// ---------------------------------------------------------------------------------------------------------------

#include <vector>
#include <stddef.h>

/** @return the maximum value in the provided vector.
    @param idx will be filled with the index where the maximum value is found.
    If values.size() == 0 then 0.0 is returned and idx is set to zero. */
double Max(const std::vector<double>& values, size_t& idx);

/** @return the maximum value in the provided vector.
    If values.size() == 0 then 0.0 is returned. */
double Max(const std::vector<double>& values);

/** @return the maximum value in the provided range of indices.
    If values.size() == 0 or fromIdx >= values.size() or length == 0 then 0.0 is returned. */
double Max(std::vector<double>::const_iterator start, std::vector<double>::const_iterator end);

/** @return the maximum absolute value in the provided vector.
    If values.size() == 0 then 0.0 is returned. */
double MaxAbs(const std::vector<double>& values);

/** @return the minimum value in the provided vector.
    @param idx will be filled with the index where the minimum value is found.
    If values.size() == 0 then 0.0 is returned and idx is set to zero. */
double Min(const std::vector<double>& values, size_t& idx);

/** @return the minimum value in the provided vector.
    If values.size() == 0 then 0.0 is returned. */
double Min(const std::vector<double>& values);

/** @return the minimum value in the provided range of indices.
    If values.size() == 0 or fromIdx >= values.size() or length == 0 then 0.0 is returned. */
double Min(std::vector<double>::const_iterator start, std::vector<double>::const_iterator end);

/** @return the minimum _and_ maximum value in the provided vector.
    If values.size() == 0 then (0.0, 0.0) is returned. */
std::pair<double, double> MinMax(const std::vector<double>& values);

/** @return the minimum _and_ maximum value in the provided vector.
    @param idx will be filled with the index where the minimum and maximum values are found.
    If values.size() == 0 then (0.0, 0.0) is returned. */
std::pair<double, double> MinMax(const std::vector<double>& values, std::pair<size_t, size_t>& idx);

/** @return the minimum _and_ maximum value in the provided range of indices.
    @param idx will be filled with the index where the minimum and maximum values are found, relative to 'start'.
    If values.size() == 0 then (0.0, 0.0) is returned. */
std::pair<double, double> MinMax(std::vector<double>::const_iterator start, std::vector<double>::const_iterator end, std::pair<size_t, size_t>& idx);

/** @return the sum of all the given values
    If values.size() == 0 then 0.0 is returned. */
double Sum(const std::vector<double>& values);

/** @return the sum of all the values between start and end.
    If values.size() == 0 then 0.0 is returned. */
double Sum(std::vector<double>::const_iterator start, std::vector<double>::const_iterator end);

/** @return the sum of the absolutes of all the given values.
    If values.size() == 0 then 0.0 is returned. */
double SumAbs(const std::vector<double>& values);

/** @return the sum the squared differences between the two provided vectors.
    If a.size() == 0 then 0.0 is returned.
    If a.size() != b.size() then -1.0 is returned. */
double SumOfSquaredDifferences(const std::vector<double>& a, const std::vector<double>& b);

/** Multiplies all values in the provided vector with the provided factor */
void Mult(std::vector<double>& values, double factor);

// Divides each value in 'first' with the corresponding value in 'second'.
// The result is stored in 'first'.
// All divisions by zero will be equal to zero.
// @throws std::invalid_argument if first.size() != second.size()
void Div(std::vector<double>& first, const std::vector<double>& second);

/** Adds the provided value to each value in the vector */
void Add(std::vector<double>& values, double factor);

/** Inverts all values in the provided vector, i.e. values[ii] = 1.0/values[ii] */
void Invert(std::vector<double>& values);

/** Reverts all values in the provided vector, such that the first value will be the last */
void Reverse(std::vector<double>& values);

/** Multiplies all values in the provided first vector with the corresponding value in the second vector.
    The results are stored in the second vector.
    @throws std::invalid_argument if the two vectors have different length. */
void Mult(const std::vector<double>& firstVector, std::vector<double>& secondVectorAndResult);

/** Calculates the exponent of the provided values */
void Exp(std::vector<double>& values);

/** Calculates the natural logarithm (base e) of the provided values */
void Log(std::vector<double>& values);

/** Performs a high-pass binomial filtering on the given values */
void HighPassBinomial(std::vector<double>& values, int nIterations);

/** Performs a low-pass binomial filtering on the given values */
void LowPassBinomial(std::vector<double>& values, int nIterations);

/** @return the Average of all the given values
    If values.size() == 0 then 0.0 is returned. */
double Average(const std::vector<double>& values);

/** @return the Average of all the values between start and end.
    If values.size() == 0 then 0.0 is returned. */
double Average(std::vector<double>::const_iterator start, std::vector<double>::const_iterator end);

/** This function calculates the variance of all the elements
    If values.size() <= 1 then 0.0 is returned. */
double Variance(const std::vector<double>& values);

/** This function calculates the standard deviation of all the elements
    If values.size() <= 1 then 0.0 is returned. */
double Stdev(const std::vector<double>& values);

/** @return the average of all the given values weighted with the inverse of each of the error.
    If values.size() == 0 or values.size() != errors.size() then 0.0 is returned. */
double WeightedAverage(const std::vector<double>& values, const std::vector<double>& errors);

/** Calculates the average of all the given values and subtracts it from the data. */
void RemoveMean(std::vector<double>& values);

/** Fits a line to all the given values and subtracts it from the data.
    Similar to 'RemoveMean' but uses a polynomial of order one instead of zero. */
void RemoveSlope(std::vector<double>& values);

/** @return the Median of all the given values.
    This will sort the provided vector.
    If values.size() == 0 then 0.0 is returned. */
double Median(std::vector<double>& values);

/** @return the area under the provided function, assuming it is sampled
    on a uniform grid with x-axis step size 'xStep'.
    If values.size() == 0 then 0.0 is returned. */
double Area(const std::vector<double>& values, double xStep);

/** Calculates the index value which corresponds to the center of mass for the
    given input dataset. */
double Centroid(const std::vector<double>& values);

/** Finds the 'N' lowest values in the input vector and fills them into the result vector
    If N == 0 or input.size() == 0 then result will be an empty vector upon return. */
void FindNLowest(const std::vector<double>& input, size_t N, std::vector<double>& result);

/** Normalizes a vector of values such that the highest value will be 1.0 and the lowest 0.0.
    If input.size() == 0 then output.size() will also be zero. */
void Normalize(const std::vector<double>& input, std::vector<double>& output);

/** Normalizes a vector of values such that the highest value will be 1.0 and the lowest 0.0. */
void Normalize(std::vector<double>& values);

/** Normalizes a vector of values such that the lowest value will be 0.0 and the
    sum of all the values will be 1.0.
    If input.size() == 0 then output.size() will also be zero. */
void NormalizeArea(const std::vector<double>& input, std::vector<double>& output);

/** Normalizes a function sampled on a uniform grid, with x-axis step size of 'xStep'
    such that the lowest value will be 0.0 and area under the graph will be 1.0.
    If input.size() == 0 then output.size() will also be zero. */
void NormalizeArea(const std::vector<double>& input, double xStep, std::vector<double>& output);

/** Finds the (fractional) index where the values in the provided vector crosses the y='valueToFind' line
    in the index-range [startIdx, stopIdx[.
    This assumes that the provided vector of values is sorted in increasing order.
     Returns -1.0 if the value cannot be found OR startIdx >= stopIdx OR startIdx >= values.size() */
double FindValue(const std::vector<double>& values, double valueToFind, size_t startIdx, size_t stopIdx);

/** Returns the value of the given vector at the fractional index.
    The result will be linearly interpolated between the lower and the upper index.
    @return 0.0 if idx < 0.0 or idx > values.size() - 1 */
double GetAt(const std::vector<double>& values, double idx);

/** @return the average resolution of the provided wavelength grid.
    This is the average delta lambda between two neigboring pixels. */
double Resolution(const std::vector<double>& wavelGrid);

/** Generates a linear mapping from minValue to maxValue with 'length' values. */
std::vector<double> GenerateVector(double minValue, double maxValue, size_t length);

/** Returns true if the provided vector contains the provided value. */
bool Contains(const std::vector<size_t>& data, size_t value);

