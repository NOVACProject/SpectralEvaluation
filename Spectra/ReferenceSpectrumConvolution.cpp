#include "ReferenceSpectrumConvolution.h"
#include <iostream>
#include <assert.h>

bool ConvolveReference(
    const std::vector<double>& /*pixelToWavelengthMapping*/,
    const std::vector<double>& /*slf*/,
    const std::vector<double>& /*highResReference*/,
    std::vector<double>& /*result*/)
{

    return true;
}

/* void CreateCVector(const SimpleSpectrum& input, MathFit::CVector& xValues, MathFit::CVector& yValues)
{
    const int inputSize = (int)input.wavelength.size();

    xValues.SetSize(inputSize);
    yValues.SetSize(inputSize);

    for (int idx = 0; idx < inputSize; ++idx)
    {
        xValues.SetAt(idx, input.wavelength[idx]);
        yValues.SetAt(idx, input.value[idx]);
    }
} */

double Max(const std::vector<double>& values, size_t& idx)
{
    double m = 0.0;
    idx = 0;

    for (size_t ii = 0; ii < values.size(); ++ii)
    {
        if (values[ii] > m)
        {
            m = values[ii];
            idx = ii;
        }
    }

    return m;
}

double Max(const std::vector<double>& values)
{
    size_t idx = 0U;
    return Max(values, idx);
}

void Normalize(const std::vector<double>& input, std::vector<double>& output)
{
    output.resize(input.size());

    const double maxValue = Max(input);

    for (size_t ii = 0; ii < input.size(); ++ii)
    {
        output[ii] = input[ii] / maxValue;
    }
}

/** Finds the first occurrence of the provided value in the given index-range of the vector. */
double FindValue(const std::vector<double>& values, double valueToFind, size_t startIdx, size_t stopIdx)
{
    if (stopIdx <= startIdx)
    {
        return (double)startIdx; // invalid range
    }

    for (size_t idx = startIdx + 1; idx < stopIdx; ++idx)
    {
        const double lastValue = values[idx - 1];
        const double thisValue = values[idx];

        if ((thisValue >= valueToFind && lastValue < valueToFind))
        {
            const double alpha = (thisValue - valueToFind) / (thisValue - lastValue);
            return (double)idx - alpha;
        }
        else if ((thisValue <= valueToFind && lastValue > valueToFind))
        {
            const double alpha = (lastValue - valueToFind) / (lastValue - thisValue);
            return (double)idx - 1 + alpha;
        }
    }

    // The vector doesn't contain the value to find.
    return (double)stopIdx;
}

double GetAt(const std::vector<double>& values, double idx)
{
    if (idx < 0.0)
    {
        return 0.0;
    }
    if (idx > (double)(values.size() - 1))
    {
        return (double)(values.size() - 1);
    }

    // linear interpolation between floor(idx) and ceil(idx)
    double x1 = values.at((int)std::floor(idx));
    double x2 = values.at((int)std::ceil(idx));

    double alpha = idx - std::floor(idx);

    return x1 * (1 - alpha) + x2 * alpha;
}

double CalculateFhwm(const SimpleSpectrum& slf)
{
    // 1. Find the maximum value
    size_t idxOfMax = 0U;
    const double maxValue = Max(slf.value, idxOfMax);

    // 2. Find the (fractional) index to the left and to the right where the amplitude has fallen to half.
    const double leftIdx  = FindValue(slf.value, maxValue * 0.5, 0U, idxOfMax - 1);
    const double rightIdx = FindValue(slf.value, maxValue * 0.5, idxOfMax + 1, slf.value.size());

    // 3. Get the x-axis values at the respective indices
    const double leftX  = GetAt(slf.wavelength, leftIdx);
    const double rightX = GetAt(slf.wavelength, rightIdx);

    // 4. We now have the FWHM
    return (rightX - leftX);
}

bool Convolve(
    const SimpleSpectrum& slf,
    const SimpleSpectrum& highResReference,
    std::vector<double>& result)
{
    if (slf.wavelength.size() != slf.value.size())
    {
        std::cout << " Error in call to 'convolve', the SLF must have as many values as wavelength values." << std::endl;
        return false;
    }
    if(highResReference.wavelength.size() != highResReference.value.size())
    {
        std::cout << " Error in call to 'convolve', the reference must have as many values as wavelength values." << std::endl;
        return false;
    }

    const size_t refSize = highResReference.value.size();
    const size_t coreSize = slf.value.size();

    // To preserve the energy, we need to normalize the slit-function to the range 0->1
    std::vector<double> normalizedSlf;
    Normalize(slf.value, normalizedSlf);
    assert(normalizedSlf.size() == slf.value.size());

    // create an intermediate result which has the correct number of values for the calculations
    std::vector<double> intermediate(refSize + coreSize - 1, 0.0);

    // The actual convolution
    for (size_t n = 0; n < refSize + coreSize - 1; ++n)
    {
        intermediate[n] = 0;

        size_t kmin = (n >= coreSize - 1) ? n - (coreSize - 1) : 0;
        size_t kmax = (n < refSize - 1) ? n : refSize - 1;

        for (size_t k = kmin; k <= kmax; k++)
        {
            intermediate[n] += highResReference.value[k] * normalizedSlf[n - k];
        }
    }

    // Cut down the result to the same size as the output is supposed to be
    result = std::vector<double>(begin(intermediate) + coreSize / 2 + 1, begin(intermediate) + refSize + coreSize / 2 + 1);

    return true;
}