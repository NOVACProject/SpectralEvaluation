#include "ReferenceSpectrumConvolution.h"
#include <iostream>
#include <algorithm>
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

double Max(const std::vector<double>& values)
{
    double m = 0.0;

    for (double v : values)
    {
        m = std::max(m, v);
    }

    return m;
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