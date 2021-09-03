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

/** @return the minimum _and_ maximum value in the provided vector.
    If values.size() == 0 then (0.0, 0.0) is returned. */
std::pair<double, double> MinMax(const std::vector<double>& values);

/** @return the minimum _and_ maximum value in the provided vector.
    @param idx will be filled with the index where the minimum and maximum values are found.
    If values.size() == 0 then (0.0, 0.0) is returned. */
std::pair<double, double> MinMax(const std::vector<double>& values, std::pair<size_t, size_t>& idx);

/** @return the sum of all the given values
    If values.size() == 0 then 0.0 is returned. */
double Sum(const std::vector<double>& values);

/** @return the sum of the absolutes of all the given values.
    If values.size() == 0 then 0.0 is returned. */
double SumAbs(const std::vector<double>& values);

/** @return the sum the squared differences between the two provided vectors.
    If a.size() == 0 then 0.0 is returned.
    If a.size() != b.size() then -1.0 is returned. */
double SumOfSquaredDifferences(const std::vector<double>& a, const std::vector<double>& b);

/** Multiplies all values in the provided vector with the provided factor */
void Mult(std::vector<double>& values, double factor);

/** Inverts all values in the provided vector */
void Invert(std::vector<double>& values);

/** Multiplies all values in the provided first vector with the corresponding value in the second vector.
    The results are stored in the second vector.
    @throws std::invalid_argument if the two vectors have different length. */
void Mult(const std::vector<double>& firstVector, std::vector<double>& secondVectorAndResult);

/** Calculates the exponent of the provided values */
void Exp(std::vector<double>& values);

/** @return the Average of all the given values
    If values.size() == 0 then 0.0 is returned. */
double Average(const std::vector<double>& values);

/** @return the average of all the given values weighted with the inverse of each of the error.
    If values.size() == 0 or values.size() != errors.size() then 0.0 is returned. */
double WeightedAverage(const std::vector<double>& values, const std::vector<double>& errors);

/** Calculates the average of all the given values and subtracts it from the data. */
void RemoveMean(std::vector<double>& values);

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

