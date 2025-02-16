#include <SpectralEvaluation/Calibration/InstrumentLineShapeEstimation.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Calibration/FraunhoferSpectrumGeneration.h>
#include <SpectralEvaluation/Calibration/CrossSectionSpectrumGenerator.h>
#include <SpectralEvaluation/Evaluation/BasicMath.h>
#include <SpectralEvaluation/Spectra/SpectrumUtils.h>
#include <SpectralEvaluation/Spectra/Scattering.h>
#include <SpectralEvaluation/Spectra/WavelengthRange.h>
#include <SpectralEvaluation/VectorUtils.h>
#include <SpectralEvaluation/Evaluation/FitWindow.h>
#include <SpectralEvaluation/Evaluation/DoasFit.h>
#include <SpectralEvaluation/Calibration/InstrumentLineShape.h>
#include <SpectralEvaluation/Math/IndexRange.h>

#undef min
#undef max

#include <assert.h>
#include <algorithm>
#include <memory>
#include <cmath>
#include <sstream>
#include <limits>

namespace novac
{

    static double GaussianSigmaToFwhm(double sigma)
    {
        return sigma * 2.0 * std::sqrt(2.0 * std::log(2.0));
    }


    /// ------------------------------------------------------------------------------------------
    /// --------------------------- InstrumentLineShapeEstimation --------------------------------
    /// ------------------------------------------------------------------------------------------

    InstrumentLineShapeEstimation::InstrumentLineShapeEstimation(const std::vector<double>& initialPixelToWavelengthMapping)
        : pixelToWavelengthMapping(initialPixelToWavelengthMapping)
    {
    }

    InstrumentLineShapeEstimation::InstrumentLineShapeEstimation(const std::vector<double>& initialPixelToWavelengthMapping, const novac::CCrossSectionData& initialLineShape)
        : pixelToWavelengthMapping(initialPixelToWavelengthMapping)
    {
        initialLineShapeEstimation = std::make_unique<novac::CCrossSectionData>(initialLineShape);
    }

    void InstrumentLineShapeEstimation::UpdateInitialLineShape(const novac::CCrossSectionData& newInitialLineShape)
    {
        this->initialLineShapeEstimation = std::make_unique<novac::CCrossSectionData>(newInitialLineShape);
    }

    bool InstrumentLineShapeEstimation::HasInitialLineShape() const
    {
        return this->initialLineShapeEstimation != nullptr && this->initialLineShapeEstimation->GetSize() > 0;
    }



    /// --------------------------------------------------------------------------------------------------------------
    /// --------------------------- InstrumentLineShapeEstimationFromKeypointDistance --------------------------------
    /// --------------------------------------------------------------------------------------------------------------

    InstrumentLineShapeEstimationFromKeypointDistance::LineShapeEstimationState InstrumentLineShapeEstimationFromKeypointDistance::EstimateInstrumentLineShape(IFraunhoferSpectrumGenerator& fraunhoferSpectrumGen, const CSpectrum& spectrum, novac::CCrossSectionData& estimatedLineShape, double& fwhm)
    {
        // Check the range over which we can compare the spectra
        const WavelengthRange comparisonWavelengthRange = fraunhoferSpectrumGen.GetFraunhoferRange(this->pixelToWavelengthMapping);
        if (comparisonWavelengthRange.Empty())
        {
            throw std::invalid_argument("Cannot create an initial estimate for the instrument line shape, no overlap between the instruments wavelength range and the solar spectrum's wavelength range.");
        }

        std::vector<double> spectrumData{ spectrum.m_data, spectrum.m_data + spectrum.m_length };
        const double maximumSpectrumValue = Max(spectrumData);
        const IndexRange goodIntensityRange = Threshold(spectrumData, 0.2 * maximumSpectrumValue);

        const IndexRange comparisonIndexRange
        {
            std::max(goodIntensityRange.from, WavelengthToPixel(this->pixelToWavelengthMapping, comparisonWavelengthRange.low)),
            std::min(goodIntensityRange.to, WavelengthToPixel(this->pixelToWavelengthMapping, comparisonWavelengthRange.high))
        };
        if (comparisonIndexRange.Empty())
        {
            throw std::invalid_argument("Cannot create an initial estimate for the instrument line shape, no overlay between the goood intensity range of the spectrum and the solar spectrum's wavelength range.");
        }

        // Prepare a high-pass filtered version of the measured spectrum, removing all baseline
        CBasicMath math;
        std::unique_ptr<CSpectrum> filteredSpectrum = std::make_unique<CSpectrum>(spectrum);
        math.LowPassBinomial(filteredSpectrum->m_data, (int)filteredSpectrum->m_length, 5);
        math.HighPassBinomial(filteredSpectrum->m_data, (int)filteredSpectrum->m_length, 500);
        Normalize(*filteredSpectrum);

        LineShapeEstimationState state;

        state.medianPixelDistanceInMeas = GetMedianKeypointDistanceFromSpectrum(*filteredSpectrum, comparisonIndexRange, "Meas");
        const double pixelDistanceFromInitialCalibration = std::abs(this->pixelToWavelengthMapping.back() - this->pixelToWavelengthMapping.front()) / (double)this->pixelToWavelengthMapping.size();
        const double medianMeasKeypointDistanceInWavelength = state.medianPixelDistanceInMeas * pixelDistanceFromInitialCalibration;

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

            medianPixelDistanceAtLowerSigmaLimit = GetMedianKeypointDistanceFromSpectrum(*solarSpectrum, comparisonIndexRange, "LowLimit");

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

            medianPixelDistanceAtUpperSigmaLimit = GetMedianKeypointDistanceFromSpectrum(*solarSpectrum, comparisonIndexRange, "HighLimit");

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

            medianPixelDistanceInSolarSpectrum = GetMedianKeypointDistanceFromSpectrum(*solarSpectrum, comparisonIndexRange, "Theory");

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
                std::abs(upperSigmaLimit - lowerSigmaLimit) < 0.1 * upperSigmaLimit)
            {
                break;
            }

            estimatedGaussianSigma = 0.5 * (lowerSigmaLimit + upperSigmaLimit);
        }

        // Save the results
        state.lineShape.center = 0.0;
        state.lineShape.sigma = estimatedGaussianSigma;
        fwhm = GaussianSigmaToFwhm(estimatedGaussianSigma);

        return state;
    }

    IndexRange InstrumentLineShapeEstimationFromKeypointDistance::Threshold(const std::vector<double>& spectrum, double threshold)
    {
        IndexRange range{ 0, 0 };
        bool foundFirstValueAboveThreshold = false;
        const size_t halfAverageLength = 7;
        for (size_t idx = halfAverageLength; idx < spectrum.size() - halfAverageLength - 1; ++idx)
        {
            const double curAverage = Average(begin(spectrum) + idx - halfAverageLength, begin(spectrum) + idx + halfAverageLength);
            if (curAverage > threshold)
            {
                if (!foundFirstValueAboveThreshold)
                {
                    range.from = idx;
                    foundFirstValueAboveThreshold = true;
                }
                range.to = idx;
            }
        }

        return range;
    }

    bool LineIntersects(const std::vector<double>& data, size_t index, double threshold)
    {
        return (data[index] > threshold && data[index - 1] < threshold) ||
            (data[index] < threshold && data[index - 1] > threshold);
    }

    double InstrumentLineShapeEstimationFromKeypointDistance::GetMedianKeypointDistanceFromSpectrum(const CSpectrum& spectrum, const IndexRange& pixelRange, const std::string& /*spectrumName*/) const
    {
        // Version 2, getting the median distance between zero crossings of the spectrum in the region [measuredPixelStart, measuredPixelStop]
        // start by normalizing the data by removing the median value
        std::vector<double> normalizedData{ spectrum.m_data + pixelRange.from, spectrum.m_data + pixelRange.to };
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

        if (intersections.size() == 0)
        {
            return 0.0;
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

        if (zeroCrossingDistances.size() > 0)
        {
            double medianPixelDistance = Average(zeroCrossingDistances);
            return medianPixelDistance;
        }
        else
        {
            return 0.0;
        }
    }

    /// <summary>
    /// Investigates the provided lineshape and returns the left and right index 
    /// defining the FWHM.
    /// If no such points could be found then (0; 0) is returned.
    /// </summary>
    std::pair<double, double> GetFwhm(const std::vector<double>& lineShape)
    {
        if (lineShape.size() <= 1)
        {
            // empty input
            return std::make_pair(0.0, 0.0);
        }

        std::pair<size_t, size_t> minMaxIdx;
        const auto minMax = MinMax(lineShape, minMaxIdx);

        // Set a threshold which is half the maximum value.
        const double threshold = minMax.first + 0.5 * (minMax.second - minMax.first);

        // from the center (position of maximum value, get the point to the left and to the right where the lineshape crosses the threshold.
        const double leftIdx = FindValue(lineShape, threshold, 0U, minMaxIdx.second);
        const double rightIdx = FindValue(lineShape, threshold, minMaxIdx.second, lineShape.size());

        if (leftIdx < 0.0 || rightIdx < leftIdx)
        {
            return std::make_pair(0.0, 0.0);
        }

        return std::make_pair(leftIdx, rightIdx);
    }

    double GetFwhm(const std::vector<double>& lineshapeWavelength, const std::vector<double>& lineShapeIntensity)
    {
        if (lineShapeIntensity.size() <= 1)
        {
            // empty input
            return 0.0;
        }

        const auto leftAndRightIdx = GetFwhm(lineShapeIntensity);

        if (leftAndRightIdx.first < 0.0 || leftAndRightIdx.second < leftAndRightIdx.first)
        {
            return 0.0;
        }

        if (lineshapeWavelength.size() > 0)
        {
            // convert index to wavelength
            return std::abs(GetAt(lineshapeWavelength, leftAndRightIdx.second) - GetAt(lineshapeWavelength, leftAndRightIdx.first));
        }
        else
        {
            return (leftAndRightIdx.second - leftAndRightIdx.first);
        }
    }

    double GetFwhm(const novac::CCrossSectionData& lineshape)
    {
        if (lineshape.m_crossSection.size() <= 1)
        {
            // empty input
            return 0.0;
        }
        return GetFwhm(lineshape.m_waveLength, lineshape.m_crossSection);
    }

    inline bool IsZero(double value)
    {
        return std::abs(value) < std::numeric_limits<double>::epsilon();
    }


    /// --------------------------------------------------------------------------------------------------
    /// --------------------------- InstrumentLineshapeEstimationFromDoas --------------------------------
    /// --------------------------------------------------------------------------------------------------

    InstrumentLineshapeEstimationFromDoas::InstrumentLineshapeEstimationFromDoas(const std::vector<double>& initialPixelToWavelengthMapping, const novac::CCrossSectionData& initialLineShape, bool addDebugOutput)
        : InstrumentLineShapeEstimation(initialPixelToWavelengthMapping, initialLineShape), debugOutput(addDebugOutput)
    {
        // Get an estimation of parameters of the Super-Gaussian from the initial line shape estimation
        if (FUNCTION_FIT_RETURN_CODE::SUCCESS != FitInstrumentLineShape(*this->initialLineShapeEstimation, this->initialLineShapeFunction))
        {
            throw InstrumentLineShapeEstimationException("Failed to fit a super gaussian line shape to measured data.");
        }
    }

    InstrumentLineshapeEstimationFromDoas::InstrumentLineshapeEstimationFromDoas(const std::vector<double>& initialPixelToWavelengthMapping, const novac::SuperGaussianLineShape& initialLineShape, bool addDebugOutput)
        : InstrumentLineShapeEstimation(initialPixelToWavelengthMapping), debugOutput(addDebugOutput)
    {
        this->initialLineShapeFunction = initialLineShape;
        this->initialLineShapeEstimation = std::make_unique<novac::CCrossSectionData>(SampleInstrumentLineShape(initialLineShape));
    }

    InstrumentLineshapeEstimationFromDoas::LineShapeEstimationResult InstrumentLineshapeEstimationFromDoas::EstimateInstrumentLineShape(
        const CSpectrum& measuredSpectrum,
        const LineShapeEstimationSettings& settings,
        IFraunhoferSpectrumGenerator& fraunhoferSpectrumGen,
        ICrossSectionSpectrumGenerator* ozoneSpectrumGen)
    {
        if (this->initialLineShapeEstimation == nullptr)
        {
            throw std::invalid_argument("Cannot estimate the instrument line shape without an initial estimate.");
        }
        if (this->pixelToWavelengthMapping.size() == 0 || this->pixelToWavelengthMapping.size() != static_cast<size_t>(measuredSpectrum.m_length))
        {
            throw std::invalid_argument("Invalid setup of Instrument Lineshape Estimation, the initial pixel-to-wavelength mapping must have the same length as the measured spectrum.");
        }
        if (settings.endPixel <= settings.startPixel)
        {
            throw std::invalid_argument("Invalid setup of Instrument Lineshape Estimation, a valid wavelength range must be defined.");
        }

        const int maximumNumberOfIterations = 100; // slightly arbitrary, but we need to stop at some point.

        double stepSize = 1.0;

        LineShapeEstimationResult result;

        SuperGaussianLineShape parameterizedLineShape = this->initialLineShapeFunction;
        if (debugOutput)
        {
            std::cout << " Initial super gaussian is (w: " << parameterizedLineShape.w << ", k: " << parameterizedLineShape.k << " giving fwhm: " << parameterizedLineShape.Fwhm() << ")" << std::endl;
        }

        InstrumentLineshapeEstimationFromDoas::LineShapeUpdate lastUpdate;
        lastUpdate.residualSize = std::numeric_limits<double>::max();
        lastUpdate.currentError = std::numeric_limits<double>::max();
        SuperGaussianLineShape lastLineShape = parameterizedLineShape;

        LineShapeEstimationAttempt optimumResult; // i.e. the best result we have ever gotten
        optimumResult.error = std::numeric_limits<double>::max();
        optimumResult.lineShape = parameterizedLineShape;

        int iterationCount = 0;
        bool successfullyConverged = false;
        bool allowSpectrumShift = true;
        while (iterationCount < maximumNumberOfIterations)
        {
            ++iterationCount;

            auto update = GetGradient(fraunhoferSpectrumGen, ozoneSpectrumGen, measuredSpectrum, parameterizedLineShape, settings, allowSpectrumShift);

            LineShapeEstimationAttempt currentAttempt;
            currentAttempt.error = update.currentError;
            currentAttempt.shift = update.shift;
            currentAttempt.lineShape = parameterizedLineShape;
            result.attempts.push_back(currentAttempt);

            if (debugOutput)
            {
                std::cout << " Super gaussian (w: " << parameterizedLineShape.w << ", k: " << parameterizedLineShape.k << ", fwhm: " << parameterizedLineShape.Fwhm() << ") gave currentError : " << update.currentError << " (with pa : " << update.residualSize << ") and shift : " << update.shift << std::endl;
                std::cout << "    Delta is (w: " << update.parameterDelta[0] << ", k: " << update.parameterDelta[1] << ") " << std::endl;
            }

            if (update.currentError < optimumResult.error)
            {
                optimumResult.error = update.currentError;
                optimumResult.shift = update.shift;
                optimumResult.lineShape = parameterizedLineShape;
            }

            if (MaxAbs(update.parameterDelta) < 1e-6 || std::abs(lastUpdate.currentError - update.currentError) < 1e-8)
            {
                // we're done
                successfullyConverged = true;
                if (debugOutput)
                {
                    std::cout << "Instrument line shape estimation successfully converged after " << iterationCount << " iterations." << std::endl;
                }
                break;
            }

            if (update.currentError > lastUpdate.currentError ||
                !UpdatedParametersAreValid(lastLineShape, update.parameterDelta, stepSize))
            {
                // The update resulted in a worse solution than the one we had earlier, or an invalid solution.
                //  Step back and make a smaller step size.
                stepSize *= 0.5;

                if (debugOutput)
                {
                    std::cout << " Updated resulted in a worse result than last iteration, reducing step size to " << stepSize << " and attempting again." << std::endl;
                }

                for (int parameterIdx = 0; parameterIdx < static_cast<int>(lastUpdate.parameterDelta.size()); ++parameterIdx)
                {
                    const double currentValue = GetParameterValue(lastLineShape, parameterIdx);
                    const double newValue = currentValue - stepSize * lastUpdate.parameterDelta[parameterIdx];
                    SetParameterValue(parameterizedLineShape, parameterIdx, newValue);
                }
            }
            else
            {
                lastLineShape = parameterizedLineShape;
                lastUpdate = update;

                // Limit how much we are allowed to update each of the parameters if the magnitude of the gradient is large.
                const double gradientSize = SumAbs(update.parameterDelta);
                const double currentStepSize = std::min(stepSize, 0.1 / gradientSize);

                // use the gradient to modify the parameterizedLineShape
                for (int parameterIdx = 0; parameterIdx < static_cast<int>(update.parameterDelta.size()); ++parameterIdx)
                {
                    const double currentValue = GetParameterValue(lastLineShape, parameterIdx);
                    SetParameterValue(parameterizedLineShape, parameterIdx, currentValue - currentStepSize * update.parameterDelta[parameterIdx]);
                }
            }
        }

        if (debugOutput && !successfullyConverged)
        {
            std::cout << "Instrument line shape estimation did not reach an optimum solution during the " << iterationCount << " iterations." << std::endl;
        }

        result.result = optimumResult;

        return result;
    }

    bool InstrumentLineshapeEstimationFromDoas::UpdatedParametersAreValid(const SuperGaussianLineShape& currentLineShape, const std::vector<double>& parameterDelta, double stepSize)
    {
        for (int parameterIdx = 0; parameterIdx < static_cast<int>(parameterDelta.size()); ++parameterIdx)
        {
            const double currentValue = GetParameterValue(currentLineShape, parameterIdx);
            const double newValue = currentValue - stepSize * parameterDelta[parameterIdx];

            if (newValue <= 0.0)
            {
                return false;
            }
        }

        return true;
    }

    InstrumentLineshapeEstimationFromDoas::LineShapeUpdate InstrumentLineshapeEstimationFromDoas::GetGradient(
        IFraunhoferSpectrumGenerator& fraunhoferSpectrumGen,
        ICrossSectionSpectrumGenerator* ozoneSpectrumGen,
        const CSpectrum& measuredSpectrum,
        const SuperGaussianLineShape& currentLineShape,
        const LineShapeEstimationSettings& settings,
        bool& allowSpectrumShift)
    {
        InstrumentLineshapeEstimationFromDoas::LineShapeUpdate update;
        bool exceptionHappened = false;
        try
        {
            update = CalculateGradientAndCurrentError(fraunhoferSpectrumGen, ozoneSpectrumGen, measuredSpectrum, currentLineShape, settings, allowSpectrumShift);
        }
        catch (DoasFitException& e)
        {
            std::cout << "Exception caught while calculating instrument line shape gradient: " << e.what() << std::endl;
            exceptionHappened = true;
        }

        // If the fit is really bad of failed entirely. Go back and attempt again.
        if (exceptionHappened || (update.residualSize > 2.0 && std::abs(update.shift) > 2.0))
        {
            if (allowSpectrumShift)
            {
                allowSpectrumShift = false;
            }
            else
            {
                throw InstrumentLineShapeEstimationException("Instrument line shape estimation failed, DOAS fit failed.");
            }

            // Chi2 > 2 is indicative of a _really_ bad DOAS fit. Attempt to do the DOAS fit again, but this time without allowing the shift to happen.
            // Also, this time we don't catch the exception. If an exception happens here then abort the instrument line shape estimation.
            update = CalculateGradientAndCurrentError(fraunhoferSpectrumGen, ozoneSpectrumGen, measuredSpectrum, currentLineShape, settings, false);
        }

        return update;
    }

    InstrumentLineshapeEstimationFromDoas::LineShapeUpdate InstrumentLineshapeEstimationFromDoas::CalculateGradientAndCurrentError(
        IFraunhoferSpectrumGenerator& fraunhoferSpectrumGen,
        ICrossSectionSpectrumGenerator* ozoneSpectrumGen,
        const CSpectrum& measuredSpectrum,
        const SuperGaussianLineShape& currentLineShape,
        const LineShapeEstimationSettings& settings,
        bool allowShift)
    {
        CBasicMath math;

        const novac::SHIFT_TYPE shiftOption = allowShift ? novac::SHIFT_TYPE::SHIFT_FREE : novac::SHIFT_TYPE::SHIFT_FIX;

        CFitWindow doasFitSetup;
        doasFitSetup.fitLow = static_cast<int>(settings.startPixel);
        doasFitSetup.fitHigh = static_cast<int>(settings.endPixel);
        doasFitSetup.polyOrder = 3;
        doasFitSetup.fitType = FIT_TYPE::FIT_POLY;

        // Sample this Super-Gaussian to get the line shape to convolve with.
        auto sampledLineShape = SampleInstrumentLineShape(currentLineShape);
        const double fwhm = currentLineShape.Fwhm();

        // Trick here: this algorithm is made faster by not calculating the fraunhofer spectrum on the _entire_ pixel-to-wavelength mapping
        // but only a wavelength range slightly larger than the one used for the DOAS fit.
        const size_t margin = 100;
        const IndexRange indexRangeToUse{
            settings.startPixel < margin ? (size_t)0 : settings.startPixel - margin,
            std::min(pixelToWavelengthMapping.size() - 1, (size_t)settings.endPixel + margin) };

        if (indexRangeToUse.from > 0)
        {
            doasFitSetup.fitLow -= static_cast<int>(indexRangeToUse.from);
            doasFitSetup.fitHigh -= static_cast<int>(indexRangeToUse.from);
        }

        const std::vector<double> selectedPixelToWavelengthMapping(begin(pixelToWavelengthMapping) + indexRangeToUse.from, begin(pixelToWavelengthMapping) + indexRangeToUse.to);

        auto currentFraunhoferSpectrum = fraunhoferSpectrumGen.GetFraunhoferSpectrum(selectedPixelToWavelengthMapping, sampledLineShape, fwhm, false);
        if (currentFraunhoferSpectrum == nullptr || currentFraunhoferSpectrum->m_length == 0)
        {
            std::stringstream message;
            message << "Failed to generate synthetic fraunhofer spectrum for instrument line shape: ";
            message << "(w: " << currentLineShape.w << ", k: " << currentLineShape.k << ")";

            throw InstrumentLineShapeEstimationException(message.str());
        }

        // Log the Fraunhofer reference (to get Optical Depth)
        std::vector<double> filteredFraunhoferSpectrum{ currentFraunhoferSpectrum->m_data, currentFraunhoferSpectrum->m_data + currentFraunhoferSpectrum->m_length };
        math.Log(filteredFraunhoferSpectrum.data(), currentFraunhoferSpectrum->m_length);

        // 1. Include the Fraunhofer reference as the first reference
        doasFitSetup.reference.clear();
        CReferenceFile fraunhofer;
        fraunhofer.m_data = std::make_unique<novac::CCrossSectionData>(filteredFraunhoferSpectrum);
        fraunhofer.m_columnOption = novac::SHIFT_TYPE::SHIFT_FIX;
        fraunhofer.m_columnValue = -1.0;
        fraunhofer.m_squeezeOption = novac::SHIFT_TYPE::SHIFT_FIX;
        fraunhofer.m_squeezeValue = 1.0;
        fraunhofer.m_shiftOption = shiftOption;
        fraunhofer.m_shiftValue = 0.0;
        doasFitSetup.reference.push_back(fraunhofer);

        // 2. Include the Ring spectrum as the second reference.
        auto ringSpectrum = std::make_unique<novac::CSpectrum>(Doasis::Scattering::CalcRingSpectrum(*currentFraunhoferSpectrum));
        std::vector<double> filteredRingSpectrum{ ringSpectrum->m_data, ringSpectrum->m_data + ringSpectrum->m_length };
        math.Log(filteredRingSpectrum.data(), ringSpectrum->m_length);
        CReferenceFile ringReference;
        ringReference.m_data = std::make_unique<novac::CCrossSectionData>(filteredRingSpectrum);
        ringReference.m_columnOption = novac::SHIFT_TYPE::SHIFT_FREE;
        ringReference.m_squeezeOption = novac::SHIFT_TYPE::SHIFT_FIX;
        ringReference.m_squeezeValue = 1.0;
        ringReference.m_shiftOption = novac::SHIFT_TYPE::SHIFT_LINK;
        ringReference.m_shiftValue = 0.0;
        doasFitSetup.reference.push_back(ringReference);

        // 3. Include the Ozone spectrum as the third reference (if required)
        if (ozoneSpectrumGen != nullptr)
        {
            CReferenceFile ozoneReference;
            auto currentOzoneSpectrum = ozoneSpectrumGen->GetCrossSection(selectedPixelToWavelengthMapping, sampledLineShape, fwhm, false);
            ozoneReference.m_data = std::make_unique<novac::CCrossSectionData>(*currentOzoneSpectrum);
            ozoneReference.m_columnOption = novac::SHIFT_TYPE::SHIFT_FREE;
            ozoneReference.m_squeezeOption = novac::SHIFT_TYPE::SHIFT_FIX;
            ozoneReference.m_squeezeValue = 1.0;
            ozoneReference.m_shiftOption = novac::SHIFT_TYPE::SHIFT_LINK;
            ozoneReference.m_shiftValue = 0.0;
            doasFitSetup.reference.push_back(ozoneReference);
        }

        // 4. Include the pseudo-absorbers derived from the derivative of the instrument-line-shape function.
        std::vector<size_t> parameterIndices;
        for (int parameterIdx = 0; parameterIdx < 2; ++parameterIdx)
        {
            auto diffSampledLineShape = std::make_unique<novac::CCrossSectionData>();
            diffSampledLineShape->m_waveLength = sampledLineShape.m_waveLength;
            diffSampledLineShape->m_crossSection = PartialDerivative(currentLineShape, sampledLineShape.m_waveLength, parameterIdx);

            auto diffFraunhofer = fraunhoferSpectrumGen.GetDifferentialFraunhoferSpectrum(selectedPixelToWavelengthMapping, *diffSampledLineShape, fwhm);

            auto pseudoAbsorber = std::make_unique<novac::CCrossSectionData>();
            pseudoAbsorber->m_crossSection.resize(currentFraunhoferSpectrum->m_length);
            pseudoAbsorber->m_waveLength = std::vector<double>(begin(selectedPixelToWavelengthMapping), end(selectedPixelToWavelengthMapping));
            for (size_t ii = 0; ii < pseudoAbsorber->m_crossSection.size(); ++ii)
            {
                pseudoAbsorber->m_crossSection[ii] = IsZero(currentFraunhoferSpectrum->m_data[ii]) ? 0.0 : diffFraunhofer->m_data[ii] / currentFraunhoferSpectrum->m_data[ii];
            }

            parameterIndices.push_back(doasFitSetup.reference.size()); // remember the index to pick out the results from

            CReferenceFile newReference;
            newReference.m_data = std::move(pseudoAbsorber);
            newReference.m_columnOption = novac::SHIFT_TYPE::SHIFT_FREE;
            newReference.m_squeezeOption = novac::SHIFT_TYPE::SHIFT_FIX;
            newReference.m_squeezeValue = 1.0;
            newReference.m_shiftOption = novac::SHIFT_TYPE::SHIFT_LINK;
            newReference.m_shiftValue = 0.0;
            doasFitSetup.reference.push_back(newReference);
        }

        DoasFit doas;
        doas.Setup(doasFitSetup);

        // Log the Measured spectrum (to get Optical Depth)
        std::vector<double> filteredMeasuredSpectrum{ measuredSpectrum.m_data + indexRangeToUse.from, measuredSpectrum.m_data + indexRangeToUse.to };
        math.Log(filteredMeasuredSpectrum.data(), static_cast<int>(filteredMeasuredSpectrum.size()));

        DoasResult doasResult;
        doas.Run(filteredMeasuredSpectrum.data(), indexRangeToUse.Length(), doasResult);

        LineShapeUpdate update;
        update.parameterDelta = std::vector<double>{ doasResult.referenceResult[parameterIndices.front()].column, doasResult.referenceResult[parameterIndices.back()].column };
        update.residualSize = doasResult.chiSquare;
        update.shift = doasResult.referenceResult[0].shift;

        // Now also estimate the currentError by doing the DOAS fit without the Pseudo-absorbers
        {
            doasFitSetup.reference[parameterIndices.front()].m_columnOption = novac::SHIFT_TYPE::SHIFT_FIX;
            doasFitSetup.reference[parameterIndices.front()].m_columnValue = 0.0;

            doasFitSetup.reference[parameterIndices.back()].m_columnOption = novac::SHIFT_TYPE::SHIFT_FIX;
            doasFitSetup.reference[parameterIndices.back()].m_columnValue = 0.0;

            doas.Setup(doasFitSetup);

            DoasResult newDoasResult;
            doas.Run(filteredMeasuredSpectrum.data(), indexRangeToUse.Length(), newDoasResult);
            update.currentError = newDoasResult.chiSquare;
        }

        return update;
    }

}