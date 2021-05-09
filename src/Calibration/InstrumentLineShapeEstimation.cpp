#include <SpectralEvaluation/Calibration/InstrumentLineShapeEstimation.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Calibration/FraunhoferSpectrumGeneration.h>
#include <SpectralEvaluation/Evaluation/BasicMath.h>
#include <SpectralEvaluation/Spectra/SpectrumUtils.h>
#include <SpectralEvaluation/VectorUtils.h>

#include <memory>

namespace novac
{
/// --------------------------- Estimating the instrument line shape --------------------------------

void InstrumentLineShapeEstimation::EstimateInstrumentLineShape(const CSpectrum& spectrum, novac::CCrossSectionData& estimatedLineShape, double& fwhm)
{
    // Prepare a high-pass filtered version of the measured spectrum, removing all baseline
    CBasicMath math;
    std::unique_ptr<CSpectrum> hpFilteredSpectrum = std::make_unique<CSpectrum>(spectrum);
    math.HighPassBinomial(hpFilteredSpectrum->m_data, (int)hpFilteredSpectrum->m_length, 500);
    Normalize(*hpFilteredSpectrum);

    const double medianPixelDistanceInMeas = GetMedianKeypointDistanceFromSpectrum(*hpFilteredSpectrum);
    const double pixelDistanceFromInitialCalibration = std::abs(this->pixelToWavelengthMapping.back() - this->pixelToWavelengthMapping.front()) / (double)this->pixelToWavelengthMapping.size();
    const double medianMeasKeypointDistanceInWavelength = medianPixelDistanceInMeas * pixelDistanceFromInitialCalibration;

    // Guess for a gaussian line shape and create a convolved solar spectrum with this setup
    std::vector<std::pair<std::string, double>> noCrossSections;
    FraunhoferSpectrumGeneration fraunhoferSpectrumGen{ this->highResolutionSolarAtlas, noCrossSections };

    double estimatedGaussianSigma = medianMeasKeypointDistanceInWavelength;
    double lowerSigmaLimit = medianMeasKeypointDistanceInWavelength * 0.1;
    double medianPixelDistanceAtLowerSigmaLimit = std::numeric_limits<double>::max();
    double upperSigmaLimit = 10.0 * medianMeasKeypointDistanceInWavelength;
    double medianPixelDistanceAtUpperSigmaLimit = 0.0;

    // Find an upper and a lower limit for the Gaussian sigma by pure search

    // Lower limit for sigma.
    while (medianPixelDistanceAtLowerSigmaLimit > medianPixelDistanceInMeas)
    {
        novac::CCrossSectionData ils;
        CreateGaussian(lowerSigmaLimit, pixelDistanceFromInitialCalibration, ils);

        auto solarSpectrum = fraunhoferSpectrumGen.GetFraunhoferSpectrum(this->pixelToWavelengthMapping, ils);
        math.HighPassBinomial(solarSpectrum->m_data, (int)solarSpectrum->m_length, 500);
        Normalize(*solarSpectrum);

        medianPixelDistanceAtLowerSigmaLimit = GetMedianKeypointDistanceFromSpectrum(*solarSpectrum);

        if (std::isnan(medianPixelDistanceAtLowerSigmaLimit) || medianPixelDistanceAtLowerSigmaLimit > medianPixelDistanceInMeas)
        {
            lowerSigmaLimit /= 10.0;
        }
    }

    // Upper limit for sigma
    while (medianPixelDistanceAtUpperSigmaLimit < medianPixelDistanceInMeas)
    {
        novac::CCrossSectionData ils;
        CreateGaussian(upperSigmaLimit, pixelDistanceFromInitialCalibration, ils);

        auto solarSpectrum = fraunhoferSpectrumGen.GetFraunhoferSpectrum(this->pixelToWavelengthMapping, ils);
        math.HighPassBinomial(solarSpectrum->m_data, (int)solarSpectrum->m_length, 500);
        Normalize(*solarSpectrum);

        medianPixelDistanceAtUpperSigmaLimit = GetMedianKeypointDistanceFromSpectrum(*solarSpectrum);

        if (std::isnan(medianPixelDistanceAtUpperSigmaLimit))
        {
            upperSigmaLimit /= 2.0;
            medianPixelDistanceAtUpperSigmaLimit = 0.0;
        }
        else if (medianPixelDistanceAtUpperSigmaLimit < medianPixelDistanceInMeas)
        {
            upperSigmaLimit *= 10.0;
        }
    }

    double medianPixelDistanceInSolarSpectrum = std::numeric_limits<double>::max();

    // Do a binary search to find the sigma which produces the best estimation of the instrument line shape
    //  This can probably be made faster by using more of the data and interpolate the "best" new guess. However, i have so far not gotten better results with that approach.
    int interationNum = 0;
    while (std::abs(medianPixelDistanceInSolarSpectrum - medianPixelDistanceInMeas) > 0.1 * medianMeasKeypointDistanceInWavelength &&
        std::abs(upperSigmaLimit - lowerSigmaLimit) > 0.01 * upperSigmaLimit)
    {
        ++interationNum;
        CreateGaussian(estimatedGaussianSigma, pixelDistanceFromInitialCalibration, estimatedLineShape);

        auto solarSpectrum = fraunhoferSpectrumGen.GetFraunhoferSpectrum(this->pixelToWavelengthMapping, estimatedLineShape);
        math.HighPassBinomial(solarSpectrum->m_data, (int)solarSpectrum->m_length, 500);
        Normalize(*solarSpectrum);

        medianPixelDistanceInSolarSpectrum = GetMedianKeypointDistanceFromSpectrum(*solarSpectrum);

        std::cout << "Gaussian width: " << estimatedGaussianSigma << " gives keypoint distance: " << medianPixelDistanceInSolarSpectrum << std::endl;

        if (medianPixelDistanceInSolarSpectrum > medianPixelDistanceInMeas)
        {
            upperSigmaLimit = estimatedGaussianSigma;
            medianPixelDistanceAtUpperSigmaLimit = medianPixelDistanceInSolarSpectrum;

            estimatedGaussianSigma = 0.5 * (lowerSigmaLimit + upperSigmaLimit);
        }
        else if (medianPixelDistanceInSolarSpectrum < medianPixelDistanceInMeas)
        {
            lowerSigmaLimit = estimatedGaussianSigma;
            medianPixelDistanceAtLowerSigmaLimit = medianPixelDistanceInSolarSpectrum;

            estimatedGaussianSigma = 0.5 * (lowerSigmaLimit + upperSigmaLimit);
        }
    }

    std::cout << "Final instrument line shape estimation gave Gaussian sigma of: " << estimatedGaussianSigma << std::endl;

    // Calculate the FWHM of the gaussian as well
    fwhm = estimatedGaussianSigma * 2.0 * std::sqrt(2.0 * std::log(2.0));
}

double InstrumentLineShapeEstimation::GetMedianKeypointDistanceFromSpectrum(const CSpectrum& spectrum) const
{
    // Find all keypoints in the prepared spectrum
    std::vector<SpectrumDataPoint> keypoints;
    FindKeypointsInSpectrum(spectrum, 0.01, keypoints);

    if (keypoints.size() <= 2)
    {
        // Cannot determine the distance if there's only one point
        return std::numeric_limits<double>::quiet_NaN();
    }

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