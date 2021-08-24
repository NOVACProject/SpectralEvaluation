#include <SpectralEvaluation/Spectra/SpectrumUtils.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Interpolation.h>
#include <SpectralEvaluation/Evaluation/BasicMath.h>
#include <SpectralEvaluation/VectorUtils.h>

#include <algorithm>
#include <cmath>

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

    CBasicMath math;
    std::vector<double> lowPassFilteredSpectrum(spectrum.m_data, spectrum.m_data + spectrum.m_length);
    math.LowPassBinomial(lowPassFilteredSpectrum.data(), spectrum.m_length, 10);

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
            while (startIdx > 1 && ddx2[startIdx] < 0.0)
            {
                --startIdx;
            }
            --startIdx;

            while (endIdx < static_cast<size_t>(spectrum.m_length) - 1 && ddx2[endIdx] < 0.0)
            {
                ++endIdx;
            }
            ++endIdx;

            const double peakHeight = lowPassFilteredSpectrum[ii] - std::max(lowPassFilteredSpectrum[startIdx], lowPassFilteredSpectrum[endIdx]);
            const size_t peakWidth = endIdx - startIdx;

            if (peakHeight > minPeakHeight && peakWidth >= minPeakWidth)
            {
                // Get the centroid position of the peak
                std::vector<double> peakValues{ lowPassFilteredSpectrum.data() + startIdx, lowPassFilteredSpectrum.data() + endIdx };
                std::vector<double> normalizedPeak;
                ::Normalize(peakValues, normalizedPeak);
                const double centroid = Centroid(normalizedPeak);
                pt.pixel = centroid + startIdx;

                // extract the intensity of the spectrum at this fractional pixel point by linear interpolation
                LinearInterpolation(spectrum.m_data, (size_t)spectrum.m_length, pt.pixel, pt.intensity);

                if (pt.intensity > minimumIntensity)
                {
                    if (spectrum.m_wavelength.size() == (size_t)spectrum.m_length)
                    {
                        LinearInterpolation(spectrum.m_wavelength, pt.pixel, pt.wavelength);
                    }
                    pt.leftPixel = static_cast<double>(startIdx);
                    pt.rightPixel = static_cast<double>(endIdx);

                    const int nearestPixel = (int)std::round(pt.pixel);
                    if (ddx2[nearestPixel + 1] + -2.0 * ddx2[nearestPixel] + ddx2[nearestPixel - 1] < 0.0) // third derivative is negative
                    {
                        pt.flatTop = true;
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
    math.LowPassBinomial(lowPassFilteredSpectrum.data(), spectrum.m_length, 10);

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

            // Find a 'basis' region for the valley by locating the area around it where ddx2 is positive.
            size_t startIdx = ii;
            size_t endIdx = ii;
            while (startIdx > 1 && ddx2[startIdx] > 0.0)
            {
                --startIdx;
            }
            --startIdx;

            while (endIdx < static_cast<size_t>(spectrum.m_length) - 1 && ddx2[endIdx] > 0.0)
            {
                ++endIdx;
            }
            ++endIdx;

            const double valleyDepth = std::min(lowPassFilteredSpectrum[startIdx], lowPassFilteredSpectrum[endIdx]) - lowPassFilteredSpectrum[ii];
            const size_t valleyWidth = endIdx - startIdx;

            if (valleyDepth > minValleyDepth && valleyWidth >= minValleyidth)
            {
                // Get the centroid position of the valley
                std::vector<double> peakValues{ lowPassFilteredSpectrum.data() + startIdx, lowPassFilteredSpectrum.data() + endIdx };
                std::vector<double> normalizedValley;
                ::Normalize(peakValues, normalizedValley);
                // Invert the valley to get a peak
                for (size_t jj = 0; jj < normalizedValley.size(); ++jj)
                {
                    normalizedValley[jj] = 1.0 - normalizedValley[jj];
                }
                const double centroid = Centroid(normalizedValley);
                pt.pixel = centroid + startIdx;

                // extract the intensity of the spectrum at this fractional pixel point by linear interpolation
                LinearInterpolation(spectrum.m_data, (size_t)spectrum.m_length, pt.pixel, pt.intensity);

                if (pt.intensity > minimumIntensity)
                {
                    if (spectrum.m_wavelength.size() == (size_t)spectrum.m_length)
                    {
                        LinearInterpolation(spectrum.m_wavelength, pt.pixel, pt.wavelength);
                    }
                    pt.leftPixel = static_cast<double>(startIdx);
                    pt.rightPixel = static_cast<double>(endIdx);

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
        pt.type = SpectrumDataPointType::Valley;
        result.push_back(pt);
    }

    for (SpectrumDataPoint& pt : peaks)
    {
        pt.type = SpectrumDataPointType::Peak;
        result.push_back(pt);
    }

    std::sort(begin(result), end(result), [](const SpectrumDataPoint& p1, const SpectrumDataPoint& p2) { return p1.pixel < p2.pixel; });
}

// TODO: Move
double Average(const double* data, size_t size)
{
    double sum = data[0];
    for (size_t ii = 1; ii < size; ++ii)
    {
        sum += data[ii];
    }
    return sum / (double)size;
}

inline double LeftWidth(const SpectrumDataPoint& point)
{
    return point.pixel - point.leftPixel;
}
inline double RightWidth(const SpectrumDataPoint& point)
{
    return point.rightPixel - point.pixel;
}

bool PeakIsFullyResolved(const CSpectrum& spectrum, const SpectrumDataPoint& point)
{
    // Extract a portion of the spectrum surrounding the found peak and calculate the second derivative there.
    const size_t leftPixel = (size_t)std::round(std::max(0.0, point.pixel - 3 * LeftWidth(point)));
    const size_t rightPixel = std::min(static_cast<size_t>(spectrum.m_length), static_cast<size_t>(std::round(point.pixel + 3 * RightWidth(point))));
    const size_t cutOutLength = rightPixel - leftPixel;

    CBasicMath math;
    std::vector<double> lowPassFilteredSpectrum(spectrum.m_data + leftPixel, spectrum.m_data + rightPixel);
    math.LowPassBinomial(lowPassFilteredSpectrum.data(), static_cast<int>(cutOutLength), 10);

    std::vector<double> ddx2;
    if (!Derivative(lowPassFilteredSpectrum.data(), cutOutLength, 2, ddx2))
    {
        return true;
    }

    // There should now be exactly two zero-crossings (of significant amplitude) in this cut out second-derivative spectrum.
    Sign lastSignOfDerivative = Sign::Positive;
    const double leastSignificantAmplitude = 0.01 * std::max(Max(ddx2), -Min(ddx2));
    int numberOfZeroCrossings = 0;
    for (size_t ii = 1; ii < ddx2.size() - 1; ++ii)
    {
        const Sign currentSignOfDerivative = SignOf(ddx2[ii], leastSignificantAmplitude);
        if (currentSignOfDerivative == Sign::Positive && lastSignOfDerivative == Sign::Negative ||
            currentSignOfDerivative == Sign::Negative && lastSignOfDerivative == Sign::Positive)
        {
            lastSignOfDerivative = currentSignOfDerivative;
            ++numberOfZeroCrossings;
        }
    }

    return (numberOfZeroCrossings < 3);
}

std::vector<SpectrumDataPoint> FilterByType(const std::vector<SpectrumDataPoint>& input, SpectrumDataPointType typeToFind)
{
    std::vector<SpectrumDataPoint> result;
    result.reserve(input.size());

    for each (auto point in input)
    {
        if (point.type == typeToFind)
        {
            result.push_back(point);
        }
    }

    return result;
}

void FindEmissionLines(const CSpectrum& spectrum, std::vector<SpectrumDataPoint>& result, bool includeNotClearlyResolvedLines)
{
    // Find the 10% lowest values and use these as a baseline
    // if we assume tha there are a small number of narrow peaks, then this will be a relatively ok level for the baseline
    std::vector<double> spectrumValues{ spectrum.m_data, spectrum.m_data + spectrum.m_length };
    std::sort(begin(spectrumValues), end(spectrumValues));

    const size_t numberOfPointsForBaseline = spectrum.m_length / 10;
    const double baseline = Average(spectrumValues.data(), numberOfPointsForBaseline);
    const double maximumIntensity = spectrumValues.back(); // the vector is now sorted in increasing order, the maximum value is in the back.

    const double threshold = baseline + 0.02 * (maximumIntensity - baseline); // threshold at 2% above the baseline.

    FindPeaks(spectrum, threshold, result);

    if (result.size() == 0)
    {
        return;
    }

    // Filter out not fully resolved emission lines.
    {
        const size_t distance = 3; // Minimum distance between two peaks (in terms of peak-half-width) for a peak to be judged to be isolated.
        for (size_t peakIdx = 0; peakIdx < result.size(); ++peakIdx)
        {
            if (peakIdx < result.size() - 1)
            {
                const double thisRightWidth = RightWidth(result[peakIdx]);
                const double nextLeftWidth = LeftWidth(result[peakIdx + 1]);
                if (result[peakIdx].pixel + distance * thisRightWidth >=
                    result[peakIdx + 1].pixel - distance * nextLeftWidth)
                {
                    // This peak overlaps the next one. Both are unresolved
                    result[peakIdx].type = SpectrumDataPointType::UnresolvedPeak;
                    result[peakIdx + 1].type = SpectrumDataPointType::UnresolvedPeak;
                    ++peakIdx;
                    continue;
                }
            }

            if (!PeakIsFullyResolved(spectrum, result[peakIdx]))
            {
                result[peakIdx].type = SpectrumDataPointType::UnresolvedPeak;
                continue;
            }

            // all tests passed, the peak seems good.
            result[peakIdx].type = SpectrumDataPointType::Peak;
        }

        if (!includeNotClearlyResolvedLines)
        {
            // Extract only the fully resolved peaks.
            result = FilterByType(result, SpectrumDataPointType::Peak);
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

bool GetEnvelope(const CSpectrum& spectrum, std::vector<double>& pixel, std::vector<double>& intensity)
{
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
