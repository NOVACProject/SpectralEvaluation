#include <SpectralEvaluation/Calibration/InstrumentLineShapeEstimation.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Calibration/FraunhoferSpectrumGeneration.h>
#include <SpectralEvaluation/Evaluation/BasicMath.h>
#include <SpectralEvaluation/Spectra/SpectrumUtils.h>
#include <SpectralEvaluation/VectorUtils.h>

#undef min
#undef max

#include <algorithm>
#include <memory>
#include <cmath>

namespace novac
{
/// --------------------------- Estimating the instrument line shape --------------------------------

static double GaussianSigmaToFwhm(double sigma)
{
    return sigma * 2.0 * std::sqrt(2.0 * std::log(2.0));
}

InstrumentLineShapeEstimation::LineShapeEstimationState InstrumentLineShapeEstimation::EstimateInstrumentLineShapeFromKeypointDistance(IFraunhoferSpectrumGenerator& fraunhoferSpectrumGen, const CSpectrum& spectrum, novac::CCrossSectionData& estimatedLineShape, double& fwhm)
{
    // Prepare a high-pass filtered version of the measured spectrum, removing all baseline
    CBasicMath math;
    std::unique_ptr<CSpectrum> filteredSpectrum = std::make_unique<CSpectrum>(spectrum);
    math.LowPassBinomial(filteredSpectrum->m_data, (int)filteredSpectrum->m_length, 5);
    math.HighPassBinomial(filteredSpectrum->m_data, (int)filteredSpectrum->m_length, 500);
    Normalize(*filteredSpectrum);

    LineShapeEstimationState state;

    state.medianPixelDistanceInMeas = GetMedianKeypointDistanceFromSpectrum(*filteredSpectrum, "Meas");
    const double pixelDistanceFromInitialCalibration = std::abs(this->pixelToWavelengthMapping.back() - this->pixelToWavelengthMapping.front()) / (double)this->pixelToWavelengthMapping.size();
    const double medianMeasKeypointDistanceInWavelength = state.medianPixelDistanceInMeas * pixelDistanceFromInitialCalibration;
    std::cout << "Measured spectrum has keypoint distance: " << state.medianPixelDistanceInMeas << std::endl;

    // Guess for a gaussian line shape and create a convolved solar spectrum with this setup
    double estimatedGaussianSigma;
    if (this->HasInitialLineShape())
    {
        estimatedGaussianSigma = GetFwhm(*this->initialLineShapeEstimation) / GaussianSigmaToFwhm(1.0);
    }
    else
    {
        estimatedGaussianSigma = medianMeasKeypointDistanceInWavelength / GaussianSigmaToFwhm(1.0);
    }
    double lowerSigmaLimit = estimatedGaussianSigma * 0.25;
    double medianPixelDistanceAtLowerSigmaLimit = std::numeric_limits<double>::max();
    double upperSigmaLimit = 4.0 * estimatedGaussianSigma;
    double medianPixelDistanceAtUpperSigmaLimit = 0.0;

    // Find an upper and a lower limit for the Gaussian sigma by pure search

    // Lower limit for sigma.
    while (medianPixelDistanceAtLowerSigmaLimit > state.medianPixelDistanceInMeas)
    {
        novac::CCrossSectionData ils;
        CreateGaussian(lowerSigmaLimit, pixelDistanceFromInitialCalibration, ils);

        auto solarSpectrum = fraunhoferSpectrumGen.GetFraunhoferSpectrum(this->pixelToWavelengthMapping, ils);
        math.HighPassBinomial(solarSpectrum->m_data, (int)solarSpectrum->m_length, 500);
        Normalize(*solarSpectrum);

        medianPixelDistanceAtLowerSigmaLimit = GetMedianKeypointDistanceFromSpectrum(*solarSpectrum, "LowLimit");

        std::cout << "Gaussian width: " << lowerSigmaLimit << " gives keypoint distance: " << medianPixelDistanceAtLowerSigmaLimit << std::endl;
        state.attempts.push_back(std::pair<double, double>{lowerSigmaLimit, medianPixelDistanceAtLowerSigmaLimit});

        if (std::isnan(medianPixelDistanceAtLowerSigmaLimit) || medianPixelDistanceAtLowerSigmaLimit > state.medianPixelDistanceInMeas)
        {
            lowerSigmaLimit /= 10.0;
        }
    }

    // Upper limit for sigma
    while (medianPixelDistanceAtUpperSigmaLimit < state.medianPixelDistanceInMeas)
    {
        novac::CCrossSectionData ils;
        CreateGaussian(upperSigmaLimit, pixelDistanceFromInitialCalibration, ils);

        auto solarSpectrum = fraunhoferSpectrumGen.GetFraunhoferSpectrum(this->pixelToWavelengthMapping, ils);
        math.HighPassBinomial(solarSpectrum->m_data, (int)solarSpectrum->m_length, 500);
        Normalize(*solarSpectrum);

        medianPixelDistanceAtUpperSigmaLimit = GetMedianKeypointDistanceFromSpectrum(*solarSpectrum, "HighLimit");

        std::cout << "Gaussian width: " << upperSigmaLimit << " gives keypoint distance: " << medianPixelDistanceAtUpperSigmaLimit << std::endl;
        state.attempts.push_back(std::pair<double, double>{upperSigmaLimit, medianPixelDistanceAtUpperSigmaLimit});

        if (std::isnan(medianPixelDistanceAtUpperSigmaLimit) || std::abs(medianPixelDistanceAtUpperSigmaLimit) < 0.1 )
        {
            upperSigmaLimit /= 2.0;
            medianPixelDistanceAtUpperSigmaLimit = 0.0;
        }
        else if (medianPixelDistanceAtUpperSigmaLimit < state.medianPixelDistanceInMeas)
        {
            upperSigmaLimit *= 10.0;
        }
    }

    double medianPixelDistanceInSolarSpectrum = std::numeric_limits<double>::max();

    // Do a binary search to find the sigma which produces the best estimation of the instrument line shape
    //  This can probably be made faster by using more of the data and interpolate the "best" new guess. However, i have so far not gotten better results with that approach.
    int interationNum = 0;
    while (true)
    {
        ++interationNum;
        CreateGaussian(estimatedGaussianSigma, pixelDistanceFromInitialCalibration, estimatedLineShape);

        auto solarSpectrum = fraunhoferSpectrumGen.GetFraunhoferSpectrum(this->pixelToWavelengthMapping, estimatedLineShape);
        math.HighPassBinomial(solarSpectrum->m_data, (int)solarSpectrum->m_length, 500);
        Normalize(*solarSpectrum);

        medianPixelDistanceInSolarSpectrum = GetMedianKeypointDistanceFromSpectrum(*solarSpectrum, "Theory");

        std::cout << "Gaussian width: " << estimatedGaussianSigma << " gives keypoint distance: " << medianPixelDistanceInSolarSpectrum << std::endl;
        state.attempts.push_back(std::pair<double, double>{estimatedGaussianSigma, medianPixelDistanceInSolarSpectrum});

        if (medianPixelDistanceInSolarSpectrum > state.medianPixelDistanceInMeas)
        {
            upperSigmaLimit = estimatedGaussianSigma;
            medianPixelDistanceAtUpperSigmaLimit = medianPixelDistanceInSolarSpectrum;
        }
        else if (medianPixelDistanceInSolarSpectrum < state.medianPixelDistanceInMeas)
        {
            lowerSigmaLimit = estimatedGaussianSigma;
            medianPixelDistanceAtLowerSigmaLimit = medianPixelDistanceInSolarSpectrum;
        }

        if (std::abs(medianPixelDistanceInSolarSpectrum - state.medianPixelDistanceInMeas) < 0.01 * state.medianPixelDistanceInMeas ||
            std::abs(upperSigmaLimit - lowerSigmaLimit) < 0.01 * upperSigmaLimit)
        {
            break;
        }

        estimatedGaussianSigma = 0.5 * (lowerSigmaLimit + upperSigmaLimit);
        // estimatedGaussianSigma = lowerSigmaLimit + (state.medianPixelDistanceInMeas - medianPixelDistanceAtLowerSigmaLimit) * (upperSigmaLimit - lowerSigmaLimit) / (medianPixelDistanceAtUpperSigmaLimit - medianPixelDistanceAtLowerSigmaLimit);
    }

    // Calculate the FWHM of the gaussian as well
    fwhm = GaussianSigmaToFwhm(estimatedGaussianSigma);

    std::cout << "Final instrument line shape estimation gave Gaussian sigma of: " << estimatedGaussianSigma << " and a fwhm of: " << fwhm << std::endl;

    return state;
}

bool InstrumentLineShapeEstimation::HasInitialLineShape() const
{
    return this->initialLineShapeEstimation != nullptr && this->initialLineShapeEstimation->GetSize() > 0;
}

bool LineIntersects(const std::vector<double>& data, size_t index, double threshold)
{
    return data[index] > threshold && data[index - 1] < threshold ||
        data[index] < threshold && data[index - 1] > threshold;
}

double InstrumentLineShapeEstimation::GetMedianKeypointDistanceFromSpectrum(const CSpectrum& spectrum, const std::string& /*spectrumName*/) const
{
    int version = 2;

    if (version == 2)
    {
        // Version 2, getting the median distance between zero crossings of the spectrum in the region [measuredPixelStart, measuredPixelStop]
        // start by normalizing the data by removing the median value
        const size_t start = (spectrum.m_length < this->measuredPixelStart) ? (spectrum.m_length / 10) : this->measuredPixelStart;
        const size_t length = (spectrum.m_length < this->measuredPixelStart) ? (spectrum.m_length - 2 * start) : std::min(static_cast<size_t>(spectrum.m_length - this->measuredPixelStart), this->measuredPixelStop - this->measuredPixelStart);
        std::vector<double> normalizedData{ spectrum.m_data + start, spectrum.m_data + start + length };
        std::vector<double> copyOfData{ begin(normalizedData), end(normalizedData) };
        auto median = Median(copyOfData);
        for (size_t ii = 0; ii < normalizedData.size(); ++ii)
        {
            normalizedData[ii] -= median;
        }
        auto minMax = MinMax(normalizedData);
        const double scaleF = 2.0 / (minMax.second - minMax.first);
        for (size_t ii = 0; ii < normalizedData.size(); ++ii)
        {
            normalizedData[ii] *= scaleF;
        }

        // Now find the locations where this normalized spectrum crosses; 1) the y=0.1 , 2) y=0.0 and 3) y=-0.1
        struct intersectionPoint
        {
            intersectionPoint(double px, int t)
                : pixel(px), type(t)
            { }

            double pixel = 0.0;
            int type = 0; // +1, 0 or -1
        };
        std::vector<intersectionPoint> intersections;
        for (size_t ii = 1; ii < normalizedData.size(); ++ii)
        {
            if (LineIntersects(normalizedData, ii, 0.0))
            {
                intersectionPoint pt(ii - 1 - normalizedData[ii - 1] / (normalizedData[ii] - normalizedData[ii - 1]), 0);
                intersections.push_back(pt);
            }
            if (LineIntersects(normalizedData, ii, +0.1))
            {
                intersectionPoint pt(ii - 1 + (0.1 - normalizedData[ii - 1]) / (normalizedData[ii] - normalizedData[ii - 1]), +1);
                intersections.push_back(pt);
            }
            if (LineIntersects(normalizedData, ii, -0.1))
            {
                intersectionPoint pt(ii - 1 + (-0.1 - normalizedData[ii - 1]) / (normalizedData[ii] - normalizedData[ii - 1]), -1);
                intersections.push_back(pt);
            }
        }

        std::vector<double> zeroCrossingDistances;
        int lastIntersectionIdx = -1;
        for (size_t ii = 1; ii < intersections.size() - 1; ++ii)
        {
            if (intersections[ii].type == 0 && intersections[ii - 1].type != intersections[ii].type && intersections[ii + 1].type != intersections[ii].type)
            {
                if (lastIntersectionIdx >= 0)
                {
                    zeroCrossingDistances.push_back(intersections[ii].pixel - intersections[lastIntersectionIdx].pixel);
                }
                lastIntersectionIdx = static_cast<int>(ii);
            }
        }

        // {
        //     std::ofstream outFile{ "D:/NOVAC/SpectrometerCalibration/ZeroCrossings.csv", std::ios::app };
        //     outFile << "------- " << spectrumName << "--------------" << std::endl;
        //     for (const auto& p : intersections)
        //     {
        //         outFile << p.pixel << std::endl;
        //     }
        // }

        double medianPixelDistance = Average(zeroCrossingDistances);
        return medianPixelDistance;
    }
    else
    {
        // Version 1, getting the median distance between keypoints.
        // Find all keypoints in the prepared spectrum
        std::vector<SpectrumDataPoint> allKeypoints;
        FindKeypointsInSpectrum(spectrum, 0.01, allKeypoints);

        // Filter out the keypoints to only be the points in the allowed range [measuredPixelStart, measuredPixelStop]
        std::vector<SpectrumDataPoint> keypoints;
        keypoints.reserve(allKeypoints.size());
        for (size_t ii = 0; ii < allKeypoints.size(); ++ii)
        {
            if (static_cast<size_t>(allKeypoints[ii].pixel) >= this->measuredPixelStart && static_cast<size_t>(allKeypoints[ii].pixel) <= this->measuredPixelStop)
            {
                keypoints.push_back(allKeypoints[ii]);
            }
        }

        if (keypoints.size() <= 2)
        {
            // Cannot determine the distance if there's only one point
            return std::numeric_limits<double>::quiet_NaN();
        }

        // {
        //     std::ofstream outFile{ "D:/NOVAC/SpectrometerCalibration/Keypoints.csv", std::ios::app };
        //     outFile << "------- " << spectrumName << "--------------" << std::endl;
        //     for (const auto& pt : keypoints)
        //     {
        //         outFile << pt.leftPixel << "\t" << pt.pixel << "\t" << pt.rightPixel << "\t" << pt.type << std::endl;
        //     }
        // }

        // Estimate the width by getting the median keypoint distance
        std::vector<double> keypointDistances(keypoints.size() - 1); // the distance between adjacent keypoints, in pixels
        std::vector<double> keypointWidth(keypoints.size()); // the width of each keypoint
        for (size_t ii = 0; ii < keypoints.size() - 1; ++ii)
        {
            keypointDistances[ii] = keypoints[ii + 1].pixel - keypoints[ii].pixel;
            keypointWidth[ii] = keypoints[ii].rightPixel - keypoints[ii].leftPixel;
        }
        keypointWidth[keypoints.size() - 1] = keypoints[keypoints.size() - 1].rightPixel - keypoints[keypoints.size() - 1].leftPixel;
        const double medianPixelDistance = Median(keypointDistances);
        const double medianPixelWidth = Median(keypointWidth);

        return medianPixelDistance;
    }
}

double GetFwhm(const novac::CCrossSectionData& lineshape)
{
    if (lineshape.m_crossSection.size() <= 1)
    {
        // empty input
        return 0.0;
    }

    std::pair<size_t, size_t> minMaxIdx;
    const auto minMax = MinMax(lineshape.m_crossSection, minMaxIdx);

    // Set a threshold which is half the maximum value.
    const double threshold = minMax.first + 0.5 * (minMax.second - minMax.first);

    // from the center (position of maximum value, get the point to the left and to the right where the lineshape crosses the threshold.
    const double leftIdx = FindValue(lineshape.m_crossSection, threshold, 0U, minMaxIdx.second);
    const double rightIdx = FindValue(lineshape.m_crossSection, threshold, minMaxIdx.second, lineshape.m_crossSection.size());

    if (leftIdx < 0.0 || rightIdx < leftIdx)
    {
        return 0.0;
    }

    if (lineshape.m_waveLength.size() > 0)
    {
        // convert index to wavelength
        return std::abs(GetAt(lineshape.m_waveLength, rightIdx) - GetAt(lineshape.m_waveLength, leftIdx));
    }
    else
    {
        return (rightIdx - leftIdx);
    }
}

}