#include <SpectralEvaluation/Calibration/InstrumentLineShapeEstimation.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Calibration/FraunhoferSpectrumGeneration.h>
#include <SpectralEvaluation/Evaluation/BasicMath.h>
#include <SpectralEvaluation/Spectra/SpectrumUtils.h>
#include <SpectralEvaluation/Spectra/Scattering.h>
#include <SpectralEvaluation/VectorUtils.h>
#include <SpectralEvaluation/Evaluation/FitWindow.h>
#include <SpectralEvaluation/Evaluation/DoasFit.h>
#include <SpectralEvaluation/Calibration/InstrumentLineShape.h>

// TODO: Remove, this is included for debugging only.
#include <SpectralEvaluation/File/File.h>

#undef min
#undef max

#include <algorithm>
#include <memory>
#include <cmath>

namespace novac
{

static double GaussianSigmaToFwhm(double sigma)
{
    return sigma * 2.0 * std::sqrt(2.0 * std::log(2.0));
}

/// ------------------------------------------------------------------------------------------
/// --------------------------- InstrumentLineShapeEstimation --------------------------------
/// ------------------------------------------------------------------------------------------

bool InstrumentLineShapeEstimation::HasInitialLineShape() const
{
    return this->initialLineShapeEstimation != nullptr && this->initialLineShapeEstimation->GetSize() > 0;
}



/// --------------------------------------------------------------------------------------------------------------
/// --------------------------- InstrumentLineShapeEstimationFromKeypointDistance --------------------------------
/// --------------------------------------------------------------------------------------------------------------

InstrumentLineShapeEstimationFromKeypointDistance::LineShapeEstimationState InstrumentLineShapeEstimationFromKeypointDistance::EstimateInstrumentLineShape(IFraunhoferSpectrumGenerator& fraunhoferSpectrumGen, const CSpectrum& spectrum, novac::CCrossSectionData& estimatedLineShape, double& fwhm)
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

        if (std::isnan(medianPixelDistanceAtUpperSigmaLimit) || std::abs(medianPixelDistanceAtUpperSigmaLimit) < 0.1)
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

    std::cout << "Final instrument line shape estimation gave Gaussian w of: " << estimatedGaussianSigma << " and a fwhm of: " << fwhm << std::endl;

    return state;
}

bool LineIntersects(const std::vector<double>& data, size_t index, double threshold)
{
    return data[index] > threshold && data[index - 1] < threshold ||
        data[index] < threshold && data[index - 1] > threshold;
}

double InstrumentLineShapeEstimationFromKeypointDistance::GetMedianKeypointDistanceFromSpectrum(const CSpectrum& spectrum, const std::string& /*spectrumName*/) const
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

inline bool IsZero(double value)
{
    return std::abs(value) < std::numeric_limits<double>::epsilon();
}


/// --------------------------------------------------------------------------------------------------
/// --------------------------- InstrumentLineshapeEstimationFromDoas --------------------------------
/// --------------------------------------------------------------------------------------------------

InstrumentLineshapeEstimationFromDoas::LineShapeEstimationResult InstrumentLineshapeEstimationFromDoas::EstimateInstrumentLineShape(IFraunhoferSpectrumGenerator& fraunhoferSpectrumGen, const CSpectrum& measuredSpectrum)
{
    if (this->initialLineShapeEstimation == nullptr)
    {
        throw std::invalid_argument("Cannot estimate the instrument line shape without an initial estimate.");
    }
    if (this->pixelToWavelengthMapping.size() == 0 || this->pixelToWavelengthMapping.size() != measuredSpectrum.m_length)
    {
        throw std::invalid_argument("Invalid setup of InstrumentLineshapeEstimationFromDoas, the initial pixel-to-wavelength mapping must have the same length as the measured spectrum.");
    }

    double stepSize = 1.0;

    LineShapeEstimationResult result;

    // Get an estimation of parameters of the Super-Gaussian from the initial line shape estimation
    SuperGaussianLineShape parameterizedLineShape;
    if (FUNCTION_FIT_RETURN_CODE::SUCCESS != FitInstrumentLineShape(*this->initialLineShapeEstimation, parameterizedLineShape))
    {
        // TODO: change type of exception
        throw std::invalid_argument("Failed to fit a super gaussian line shape to measured data.");
    }

    std::cout << " Initial super gaussian is (w: " << parameterizedLineShape.w << ", k: " << parameterizedLineShape.k << ")" << std::endl;

    InstrumentLineshapeEstimationFromDoas::LineShapeUpdate lastUpdate;
    lastUpdate.residualSize = std::numeric_limits<double>::max();
    lastUpdate.error = std::numeric_limits<double>::max();
    SuperGaussianLineShape lastLineShape = parameterizedLineShape;

    LineShapeEstimationAttempt optimumResult;
    optimumResult.error = std::numeric_limits<double>::max();
    optimumResult.shift = 0.0;
    optimumResult.lineShape = parameterizedLineShape;

    int iterationCount = 0;
    while (iterationCount < 100) // TODO: Determine a limit here
    {
        ++iterationCount;

        auto update = CalculateGradient(fraunhoferSpectrumGen, measuredSpectrum, parameterizedLineShape);

        LineShapeEstimationAttempt currentAttempt;
        currentAttempt.error = update.error;
        currentAttempt.shift = update.shift;
        currentAttempt.lineShape = parameterizedLineShape;
        result.attempts.push_back(currentAttempt);

        std::cout << " Super gaussian (w: " << parameterizedLineShape.w << ", k: " << parameterizedLineShape.k << ") gave error: " << update.error << " (with pa: " << update.residualSize << ")" << std::endl;
        std::cout << "    Delta is (w: " << update.parameterDelta[0] << ", k: " << update.parameterDelta[1] << ") " << std::endl;

        if (update.error < optimumResult.error)
        {
            optimumResult.error = update.error;
            optimumResult.shift = update.shift;
            optimumResult.lineShape = parameterizedLineShape;
        }

        if (Max(update.parameterDelta) < 1e-6 || std::abs(lastUpdate.error - update.error) < 1e-6)
        {
            // we're done
            break;
        }

        lastLineShape = parameterizedLineShape;
        lastUpdate = update;

        // use the gradient to modify the parameterizedLineShape
        for (int parameterIdx = 0; parameterIdx < static_cast<int>(update.parameterDelta.size()); ++parameterIdx)
        {
            const double currentValue = GetParameterValue(lastLineShape, parameterIdx);
            SetParameterValue(parameterizedLineShape, parameterIdx, currentValue + stepSize * update.parameterDelta[parameterIdx]);
        }
    }

    result.result = optimumResult;

    return result;
}

InstrumentLineshapeEstimationFromDoas::LineShapeUpdate InstrumentLineshapeEstimationFromDoas::CalculateGradient(
    IFraunhoferSpectrumGenerator& fraunhoferSpectrumGen,
    const CSpectrum& measuredSpectrum,
    const SuperGaussianLineShape& currentLineShape)
{
    CBasicMath math;

    std::vector<double> instrumentLineShapeWavelength;
    {
        // Super-sample the SLF, in order to not miss any finer details in the derivatives.
        const size_t length = 1 + 2 * this->initialLineShapeEstimation->m_waveLength.size();
        const double range = this->initialLineShapeEstimation->m_waveLength.back() - this->initialLineShapeEstimation->m_waveLength.front();
        const double minValue = -range * 0.5;
        const double delta = range / static_cast<double>(length - 1);
        for (size_t ii = 0; ii < length; ++ii)
        {
            instrumentLineShapeWavelength.push_back(minValue + ii * delta);
        }
    }

    CFitWindow doasFitSetup;
    doasFitSetup.fitLow = 1200;     // TODO: Input parameter
    doasFitSetup.fitHigh = 1500;    // TODO: Input parameter
    doasFitSetup.polyOrder = 3;

    // Sample this Super-Gaussian to get the line shape to convolve with.
    auto superGaussianLineShape = std::make_unique<novac::CCrossSectionData>();
    superGaussianLineShape->m_waveLength = instrumentLineShapeWavelength;
    superGaussianLineShape->m_crossSection = SampleInstrumentLineShape(currentLineShape, instrumentLineShapeWavelength, 0.0, 1.0);
    NormalizeArea(superGaussianLineShape->m_crossSection, superGaussianLineShape->m_crossSection);
    const double fwhm = GetFwhm(*superGaussianLineShape);

    auto currentFraunhoferSpectrum = fraunhoferSpectrumGen.GetFraunhoferSpectrum(this->pixelToWavelengthMapping, *superGaussianLineShape, false);

    // Log the Fraunhofer reference (to get Optical Depth)
    std::vector<double> filteredFraunhoferSpectrum{ currentFraunhoferSpectrum->m_data, currentFraunhoferSpectrum->m_data + currentFraunhoferSpectrum->m_length };
    math.Log(filteredFraunhoferSpectrum.data(), currentFraunhoferSpectrum->m_length);

    // 1. Include the Fraunhofer reference as the first reference
    doasFitSetup.nRef = 1;
    doasFitSetup.ref[0].m_data = std::make_unique<novac::CCrossSectionData>(filteredFraunhoferSpectrum);
    doasFitSetup.ref[0].m_columnOption = novac::SHIFT_FIX;
    doasFitSetup.ref[0].m_columnValue = 1.0;
    doasFitSetup.ref[0].m_squeezeOption = novac::SHIFT_FIX;
    doasFitSetup.ref[0].m_squeezeValue = 1.0;
    doasFitSetup.ref[0].m_shiftOption = novac::SHIFT_FREE;
    doasFitSetup.ref[0].m_shiftValue = 0.0;

    // 2. Include the Ring spectrum as the second reference
    auto ringSpectrum = Doasis::Scattering::CalcRingSpectrum(*currentFraunhoferSpectrum);
    std::vector<double> filteredRingSpectrum{ ringSpectrum.m_data, ringSpectrum.m_data + ringSpectrum.m_length };
    math.Log(filteredRingSpectrum.data(), ringSpectrum.m_length);
    doasFitSetup.ref[doasFitSetup.nRef].m_data = std::make_unique<novac::CCrossSectionData>(filteredRingSpectrum);
    doasFitSetup.ref[doasFitSetup.nRef].m_columnOption = novac::SHIFT_FREE;
    doasFitSetup.ref[doasFitSetup.nRef].m_squeezeOption = novac::SHIFT_FIX;
    doasFitSetup.ref[doasFitSetup.nRef].m_squeezeValue = 1.0;
    doasFitSetup.ref[doasFitSetup.nRef].m_shiftOption = novac::SHIFT_LINK;
    doasFitSetup.ref[doasFitSetup.nRef].m_shiftValue = 0.0;
    doasFitSetup.nRef += 1;

    // 3. Include the pseudo-absorbers derived from the derivative of the instrument-line-shape function.
    for (int parameterIdx = 0; parameterIdx < 2; ++parameterIdx)
    {
        auto diffSampledLineShape = std::make_unique<novac::CCrossSectionData>();
        diffSampledLineShape->m_waveLength = instrumentLineShapeWavelength;
        diffSampledLineShape->m_crossSection = PartialDerivative(currentLineShape, instrumentLineShapeWavelength, parameterIdx);

        auto diffFraunhofer = fraunhoferSpectrumGen.GetDifferentialFraunhoferSpectrum(this->pixelToWavelengthMapping, *diffSampledLineShape, fwhm);

        auto pseudoAbsorber = std::make_unique<novac::CCrossSectionData>();
        pseudoAbsorber->m_crossSection.resize(currentFraunhoferSpectrum->m_length);
        pseudoAbsorber->m_waveLength = std::vector<double>(begin(this->pixelToWavelengthMapping), end(this->pixelToWavelengthMapping));
        for (size_t ii = 0; ii < pseudoAbsorber->m_crossSection.size(); ++ii)
        {
            pseudoAbsorber->m_crossSection[ii] = IsZero(currentFraunhoferSpectrum->m_data[ii]) ? 0.0 : diffFraunhofer->m_data[ii] / currentFraunhoferSpectrum->m_data[ii];
        }

        if (parameterIdx == 0)
        {
            novac::SaveDataToFile("D:/NOVAC/SpectrometerCalibration/PseudoAbsorberSpectrum0.txt", pseudoAbsorber->m_crossSection);
        }
        else
        {
            novac::SaveDataToFile("D:/NOVAC/SpectrometerCalibration/PseudoAbsorberSpectrum1.txt", pseudoAbsorber->m_crossSection);
        }

        doasFitSetup.ref[doasFitSetup.nRef].m_data = std::move(pseudoAbsorber);
        doasFitSetup.ref[doasFitSetup.nRef].m_columnOption = novac::SHIFT_FREE;
        doasFitSetup.ref[doasFitSetup.nRef].m_squeezeOption = novac::SHIFT_FIX;
        doasFitSetup.ref[doasFitSetup.nRef].m_squeezeValue = 1.0;
        doasFitSetup.ref[doasFitSetup.nRef].m_shiftOption = novac::SHIFT_LINK;
        doasFitSetup.ref[doasFitSetup.nRef].m_shiftValue = 0.0;
        doasFitSetup.nRef += 1;
    }

    DoasFit doas;
    doas.Setup(doasFitSetup);

    // Log the Measured spectrum (to get Optical Depth)
    std::vector<double> filteredMeasuredSpectrum{ measuredSpectrum.m_data, measuredSpectrum.m_data + measuredSpectrum.m_length };
    math.Log(filteredMeasuredSpectrum.data(), measuredSpectrum.m_length);

    DoasResult doasResult;
    // TODO: this throws an exception if the fit fails. Catch it!
    doas.Run(filteredMeasuredSpectrum.data(), static_cast<size_t>(measuredSpectrum.m_length), doasResult);

    // Save the setup such that we can debug and inspect
    {
        novac::SaveDataToFile("D:/NOVAC/SpectrometerCalibration/FraunhoferSpectrum.txt", filteredFraunhoferSpectrum);
        novac::SaveDataToFile("D:/NOVAC/SpectrometerCalibration/MeasuredSpectrum.txt", filteredMeasuredSpectrum);
        novac::SaveDataToFile("D:/NOVAC/SpectrometerCalibration/RingSpectrum.txt", filteredRingSpectrum);
    }

    LineShapeUpdate update;
    update.parameterDelta = std::vector<double>{ doasResult.referenceResult[2].column, doasResult.referenceResult[3].column };
    update.residualSize = doasResult.chiSquare;
    update.shift = doasResult.referenceResult[0].shift;

    // Now also estimate the error by doing the DOAS fit without the Pseudo-absorbers
    {
        doasFitSetup.ref[2].m_columnOption = novac::SHIFT_FIX;
        doasFitSetup.ref[2].m_columnValue = 0.0;

        doasFitSetup.ref[3].m_columnOption = novac::SHIFT_FIX;
        doasFitSetup.ref[3].m_columnValue = 0.0;

        doas.Setup(doasFitSetup);

        DoasResult newDoasResult;
        doas.Run(filteredMeasuredSpectrum.data(), static_cast<size_t>(measuredSpectrum.m_length), newDoasResult);
        update.error = newDoasResult.chiSquare;
    }


    return update;
}

}