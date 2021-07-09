#include <SpectralEvaluation/Spectra/SpectrumUtils.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Interpolation.h>
#include <SpectralEvaluation/Evaluation/BasicMath.h>
#include <cassert>
#include <algorithm>

namespace novac
{

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

/// <summary>
/// Calculates and returns the center of mass of the values in the region [firstIndex, lastIndex]
/// </summary>
double CenterOfMass(const double* values, size_t firstIndex, size_t lastIndex, double baseline = 0.0)
{
    double weights = 0.0;
    double weightedSum = 0.0;
    for (size_t ii = firstIndex; ii < lastIndex; ++ii)
    {
        const double weight = std::abs(values[ii] - baseline);
        weightedSum += ii * weight;
        weights += weight;
    }
    return weightedSum / weights;
}

void FindPeaks(const CSpectrum& spectrum, double minimumIntensity, std::vector<SpectrumDataPoint>& result)
{
    result.clear();

    if (spectrum.m_length < 3)
    {
        return;
    }

    const double intensityRange = spectrum.MaxValue() - spectrum.MinValue();
    const double minPeakHeight = intensityRange / 3000.0; // TODO: Find a reasonable value here!
    const size_t minPeakWidth = 5; // TODO: Find a reasonable value here!

    std::vector<double> lowPassFilteredSpectrum(spectrum.m_data, spectrum.m_data + spectrum.m_length);
    CBasicMath math;
    math.LowPassBinomial(lowPassFilteredSpectrum.data(), spectrum.m_length, 3);

    // Calculate the first and second order derivatives
    std::vector<double> ddx;
    if (!Derivative(lowPassFilteredSpectrum.data(), spectrum.m_length, ddx))
    {
        return;
    }
    std::vector<double> ddx2;
    if (!Derivative(lowPassFilteredSpectrum.data(), spectrum.m_length, 2, ddx2))
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

            // Find a 'basis' region for the peak by locating the area around it where ddx2 is negative.
            size_t startIdx = ii;
            size_t endIdx = ii;
            while (startIdx > 0 && ddx2[startIdx] < 0.0)
            {
                --startIdx;
            }
            while (endIdx < spectrum.m_length && ddx2[endIdx] < 0.0)
            {
                ++endIdx;
            }

            const double peakHeight = lowPassFilteredSpectrum[ii] - std::max(lowPassFilteredSpectrum[startIdx], lowPassFilteredSpectrum[endIdx]);
            const size_t peakWidth = endIdx - startIdx;

            if (peakHeight > minPeakHeight && peakWidth >= minPeakWidth)
            {
                // do a small linear interpolation to find the location of the zero crossing more precisely
                const double alpha = ddx[idxOfLastSignificantDerivativeSign] / (ddx[idxOfLastSignificantDerivativeSign] - ddx[ii]);
                const double idx = (double)idxOfLastSignificantDerivativeSign + alpha * ((double)ii - (double)idxOfLastSignificantDerivativeSign);

                assert(idx > (double)idxOfLastSignificantDerivativeSign && idx < (double)ii);

                pt.pixel = idx;

                // extract the intensity of the spectrum at this fractional pixel point by linear interpolation
                LinearInterpolation(spectrum.m_data, (size_t)spectrum.m_length, idx, pt.intensity);

                if (pt.intensity > minimumIntensity)
                {
                    if (spectrum.m_wavelength.size() == (size_t)spectrum.m_length)
                    {
                        LinearInterpolation(spectrum.m_wavelength, idx, pt.wavelength);
                    }

                    result.push_back(pt);
                }
            }
        }

        if (currentSignOfDerivative != Sign::Undetermined)
        {
            lastSignOfDerivative = currentSignOfDerivative;
            idxOfLastSignificantDerivativeSign = ii;
        }
    }
}

// TODO: Combine FindPeaks and FindValleys into one method, since they do the same processing...
void FindValleys(const CSpectrum& spectrum, double minimumIntensity, std::vector<SpectrumDataPoint>& result)
{
    result.clear();

    if (spectrum.m_length < 3)
    {
        return;
    }

    const double intensityRange = spectrum.MaxValue() - spectrum.MinValue();
    const double minValleyDepth = intensityRange / 3000.0; // TODO: Find a reasonable value here!
    const size_t minValleyidth = 5; // TODO: Find a reasonable value here!

    std::vector<double> lowPassFilteredSpectrum(spectrum.m_data, spectrum.m_data + spectrum.m_length);
    CBasicMath math;
    math.LowPassBinomial(lowPassFilteredSpectrum.data(), spectrum.m_length, 3);

    // Calculate the first and second order derivatives
    std::vector<double> ddx;
    if (!Derivative(lowPassFilteredSpectrum.data(), spectrum.m_length, ddx))
    {
        return;
    }
    std::vector<double> ddx2;
    if (!Derivative(lowPassFilteredSpectrum.data(), spectrum.m_length, 2, ddx2))
    {
        return;
    }

    // Locate all points where the derivative changes sign from negative to positive
    Sign lastSignOfDerivative = Sign::Undetermined;
    size_t idxOfLastSignificantDerivativeSign = 0; // the idx where the sign of the derivative was last known
    for (size_t ii = 1; ii < (size_t)spectrum.m_length - 1; ++ii)
    {
        Sign currentSignOfDerivative = SignOf(ddx[ii], 1e-5);
        if (currentSignOfDerivative == Sign::Positive && lastSignOfDerivative == Sign::Negative)
        {
            SpectrumDataPoint pt;

            // Find a 'basis' region for the peak by locating the area around it where ddx2 is positive.
            size_t startIdx = ii;
            size_t endIdx = ii;
            while (startIdx > 0 && ddx2[startIdx] > 0.0)
            {
                --startIdx;
            }
            while (endIdx < spectrum.m_length && ddx2[endIdx] > 0.0)
            {
                ++endIdx;
            }

            const double valleyDepth = std::min(lowPassFilteredSpectrum[startIdx], lowPassFilteredSpectrum[endIdx]) - lowPassFilteredSpectrum[ii];
            const size_t valleyWidth = endIdx - startIdx;

            if (valleyDepth > minValleyDepth && valleyWidth >= minValleyidth)
            {
                // do a small linear interpolation to find the location of the zero crossing more precisely
                const double alpha = ddx[idxOfLastSignificantDerivativeSign] / (ddx[idxOfLastSignificantDerivativeSign] - ddx[ii]);
                const double idx = (double)idxOfLastSignificantDerivativeSign + alpha * ((double)ii - (double)idxOfLastSignificantDerivativeSign);

                assert(idx > (double)idxOfLastSignificantDerivativeSign && idx < (double)ii);

                pt.pixel = idx;

                // extract the intensity of the spectrum at this fractional pixel point by linear interpolation
                LinearInterpolation(spectrum.m_data, (size_t)spectrum.m_length, idx, pt.intensity);

                if (pt.intensity > minimumIntensity)
                {
                    if (spectrum.m_wavelength.size() == (size_t)spectrum.m_length)
                    {
                        LinearInterpolation(spectrum.m_wavelength, idx, pt.wavelength);
                    }

                    result.push_back(pt);
                }
            }
        }

        if (currentSignOfDerivative != Sign::Undetermined)
        {
            lastSignOfDerivative = currentSignOfDerivative;
            idxOfLastSignificantDerivativeSign = ii;
        }
    }
}

void FindKeypointsInSpectrum(const CSpectrum& spectrum, double minimumIntensity, std::vector<SpectrumDataPoint>& result)
{
    // Locate all peaks and valleys in the measured spectrum
    //  These have the correct pixel position (for the spectrometer we're trying to calibrate)
    std::vector< SpectrumDataPoint> peaks;
    FindPeaks(spectrum, minimumIntensity, peaks);

    std::vector< SpectrumDataPoint> valleys;
    FindValleys(spectrum, minimumIntensity, valleys);

    if (peaks.size() == 0 || valleys.size() == 0)
    {
        std::cout << "Failed to find peaks/valleys in the measured spectrum." << std::endl;
        result.clear();
        return;
    }

    result.clear();
    result.reserve(valleys.size() + peaks.size());

    for (SpectrumDataPoint& pt : valleys)
    {
        pt.type = -1;
        result.push_back(pt);
    }

    for (SpectrumDataPoint& pt : peaks)
    {
        pt.type = +1;
        result.push_back(pt);
    }

    std::sort(begin(result), end(result), [](const SpectrumDataPoint& p1, const SpectrumDataPoint& p2) { return p1.pixel < p2.pixel; });
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

bool Derivative(const double* data, size_t dataLength, int order, std::vector<double>& result)
{
    if (order == 1)
    {
        return Derivative(data, dataLength, result);
    }
    else if (order == 2)
    {
        if (data == nullptr || dataLength < 3 || dataLength > 1024 * 1024) // check such that the data isn't unreasonably large
        {
            return false;
        }

        result.resize(dataLength);
        result[0] = 0.0;

        for (size_t ii = 1; ii < dataLength - 2; ++ii)
        {
            result[ii] = data[ii + 1] - 2 * data[ii] + data[ii - 1];
        }
        result[dataLength - 1] = 0.0;

        return true;
    }
    else
    {
        return false;
    }
}

struct Point
{
    Point(double xLocation, double yLocation)
        : x(xLocation), y(yLocation)
    {
    }
    double x;
    double y;
};

// This is not an optimal way of calculating the envelope but does its job
std::vector<Point> Reduce(const std::vector<Point>& selectedPoints, const double* const data)
{
    std::vector<Point> reducedPoints;
    reducedPoints.push_back(selectedPoints[0]);
    for (size_t ii = 1; ii < selectedPoints.size() - 1; ++ii)
    {
        // Draw a line from selectedPoints[ii - 1] to selectedPoints[ii + 1]. 
        //  If this line goes _above_ all the intermediate points in the data then selectedPoints[ii] can be removed.
        bool intersectsData = false;
        const size_t startIdx = static_cast<size_t>(selectedPoints[ii - 1].x);
        const size_t stopIdx = static_cast<size_t>(selectedPoints[ii + 1].x);
        for (size_t pixelIdx = startIdx; pixelIdx < stopIdx; ++pixelIdx)
        {
            const double alpha = (pixelIdx - selectedPoints[ii - 1].x) / (selectedPoints[ii + 1].x - selectedPoints[ii - 1].x);
            const double y = selectedPoints[ii - 1].y + alpha * (selectedPoints[ii + 1].y - selectedPoints[ii - 1].y) / (selectedPoints[ii + 1].x - selectedPoints[ii - 1].x);
            if (y < data[pixelIdx])
            {
                intersectsData = true;
                break;
            }
        }

        if (intersectsData)
        {
            reducedPoints.push_back(selectedPoints[ii]);
        }
    }
    reducedPoints.push_back(selectedPoints[selectedPoints.size() - 1]);

    return reducedPoints;
}


bool GetEnvelope(const CSpectrum& spectrum, std::vector<double>& pixel, std::vector<double>& intensity)
{
    // Version 1. Selecting peaks
    /* std::vector<Point> selectedPoints;

    // Select all local maxima
    for (size_t ii = 1; ii < static_cast<unsigned long long>(spectrum.m_length) - 1; ++ii)
    {
        if (spectrum.m_data[ii] >= spectrum.m_data[ii - 1] && spectrum.m_data[ii] > spectrum.m_data[ii + 1])
        {
            selectedPoints.push_back(Point{ static_cast<double>(ii), spectrum.m_data[ii] });
        }
    }

    // Now reduce the points by seeing which point can be skipped by drawing a line from the point before and after (but don't remove the first and last point)
    std::vector<Point> reducedPoints;
    while (true)
    {
        reducedPoints = Reduce(selectedPoints, spectrum.m_data);

        if (reducedPoints.size() < selectedPoints.size())
        {
            selectedPoints = reducedPoints;
        }
        else
        {
            break;
        }
    }

    pixel.clear();
    intensity.clear();
    pixel.reserve(selectedPoints.size());
    intensity.reserve(selectedPoints.size());
    // always include the first point
    pixel.push_back(0);
    intensity.push_back(spectrum.m_data[0]);
    for (size_t ii = 1; ii < selectedPoints.size(); ++ii)
    {
        pixel.push_back(selectedPoints[ii].x);
        intensity.push_back(selectedPoints[ii].y);
    }
    // always include the last point
    pixel.push_back(static_cast<double>(spectrum.m_length) - 1);
    intensity.push_back(spectrum.m_data[spectrum.m_length - 1]);

    return true; */

    // Version 2. Selecting _all_ peaks and low pass filtering the resulting data
    pixel.clear();
    intensity.clear();

    // Select all local maxima
    for (size_t ii = 1; ii < static_cast<unsigned long long>(spectrum.m_length) - 1; ++ii)
    {
        if (spectrum.m_data[ii] >= spectrum.m_data[ii - 1] && spectrum.m_data[ii] > spectrum.m_data[ii + 1])
        {
            pixel.push_back(static_cast<double>(ii));
            intensity.push_back(spectrum.m_data[ii]);
        }
    }

    // Low pass filter to reduce fast variations
    CBasicMath basicMath;
    basicMath.LowPassBinomial(intensity.data(), static_cast<int>(intensity.size()), 15);

    // done.
    return true;

}

}
