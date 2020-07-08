#include <SpectralEvaluation/Spectra/SpectrumUtils.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Interpolation.h>
#include <cassert>

enum class Sign
{
    Undetermined = 0,
    Positive,
    Negative
};

/**
 * @brief Returns the sign of the provided value. This is Undeterminied if abs(value) < nearZero
 * @param value The numerical value to calculate the sign of.
 * @param nearZero A limit for values close to zero, this must be greater than 0
 * @return The sign of the provided value
*/
inline Sign SignOf(double value, double nearZero)
{
    // Assuming nearZero > epsilon
    return value > nearZero ? Sign::Positive : (value < -nearZero ? Sign::Negative : Sign::Undetermined);
}

void FindPeaks(const CSpectrum& spectrum, double minimumPeakIntensity, std::vector<SpectrumDataPoint>& result)
{
    result.clear();

    if (spectrum.m_length < 3)
    {
        return;
    }

    // Calculate the first and second order derivatives
    std::vector<double> ddx;
    if (!Derivative(spectrum.m_data, spectrum.m_length, ddx))
    {
        return;
    }

    // Locate all points where the derivative changes sign from positive to negative
    Sign lastSignOfDerivative = Sign::Undetermined;
    size_t idxOfLastSignificantDerivativeSign = 0; // the idx where the sign of the derivative was last known
    for (size_t ii = 1; ii < (size_t)spectrum.m_length - 1; ++ii)
    {
        Sign currentSignOfDerivative = SignOf(ddx[ii], 1e-5);
        if (currentSignOfDerivative == Sign::Negative && lastSignOfDerivative == Sign::Positive)
        {
            SpectrumDataPoint pt;

            // do a small linear interpolation to find the location of the zero crossing more precisely
            const double alpha = ddx[idxOfLastSignificantDerivativeSign] / (ddx[idxOfLastSignificantDerivativeSign] - ddx[ii]);
            const double idx = (double)idxOfLastSignificantDerivativeSign + alpha * ((double)ii - (double)idxOfLastSignificantDerivativeSign);

            assert(idx > (double)idxOfLastSignificantDerivativeSign && idx < (double)ii);

            pt.pixel = idx;
            if (spectrum.m_wavelength.size() == (size_t)spectrum.m_length)
            {
                LinearInterpolation(spectrum.m_wavelength, idx, pt.wavelength);
            }
            LinearInterpolation(spectrum.m_data, (size_t)spectrum.m_length, idx, pt.intensity);

            if (pt.intensity > minimumPeakIntensity)
            {
                result.push_back(pt);
            }
        }

        if (currentSignOfDerivative != Sign::Undetermined)
        {
            lastSignOfDerivative = currentSignOfDerivative;
            idxOfLastSignificantDerivativeSign = ii;
        }
    }

}


bool Derivative(const double* data, size_t dataLength, std::vector<double>& result)
{
    if (data == nullptr || dataLength < 3 || dataLength > 1024 * 1024) // check such that the data isn't unreasonably large
    {
        return false;
    }

    result.resize(dataLength);
    result[0] = 0.0;

    for (size_t ii = 1; ii < dataLength - 2; ++ii)
    {
        result[ii] = data[ii + 1] - data[ii - 1];
    }
    result[dataLength - 1] = 0.0;

    return true;
}
