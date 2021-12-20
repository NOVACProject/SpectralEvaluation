#pragma once

#include <memory>
#include <vector>
#include <string>
#include <SpectralEvaluation/Calibration/InstrumentLineShape.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>

namespace novac
{
    class CSpectrum;
    class IFraunhoferSpectrumGenerator;
    class ICrossSectionSpectrumGenerator;
    class DoasFit;
    struct IndexRange;

    class InstrumentLineShapeEstimationException : public std::exception
    {
    private:
        const char* const m_msg = "";

    public:
        InstrumentLineShapeEstimationException(const char* msg) :
            m_msg(msg)
        {}

        InstrumentLineShapeEstimationException(const std::string& msg) :
            m_msg(msg.c_str())
        {}

        const char* what() const noexcept override final { return m_msg; }
    };

    /** Abstract base class for the different instrument line shape estimators. */
    class InstrumentLineShapeEstimation
    {
    public:

        /** Sets the pixel to wavelength mapping for the spectrometer,
             this is the inital state before doing the calibration. */
        void UpdateInitialCalibration(const std::vector<double>& newPixelToWavelengthMapping)
        {
            this->pixelToWavelengthMapping = newPixelToWavelengthMapping;
        }

        /** Sets initial estimation of the instrument line shape. Used as a starting point for the estimation routine. */
        void UpdateInitialLineShape(const novac::CCrossSectionData& newInitialLineShape);

        bool HasInitialLineShape() const;

    protected:

        InstrumentLineShapeEstimation(const std::vector<double>& initialPixelToWavelengthMapping);

        InstrumentLineShapeEstimation(const std::vector<double>& initialPixelToWavelengthMapping, const novac::CCrossSectionData& initialLineShape);

        /** The assumed pixel-to-wavelength mapping for the device. */
        std::vector<double> pixelToWavelengthMapping;

        /** The initial, starting guess for the instrument line shape. */
        std::unique_ptr<novac::CCrossSectionData> initialLineShapeEstimation;
    };

    /** This is a helper class for estimating the instrument line shape of an instrument
        using a measured spectrum with a (reasonably well known) pixel-to-wavelength calibration
        by measuring the distance between keypoints in the meaured and synthetic spectra.
        This is a rather blunt tool and may be used to get an initial estimate if the instrument line shape is not at all known. */
    class InstrumentLineShapeEstimationFromKeypointDistance : public InstrumentLineShapeEstimation
    {
    public:
        /** A helper structure to describe the internal state of this estimator */
        struct LineShapeEstimationState
        {
            double medianPixelDistanceInMeas = 0.0;

            GaussianLineShape lineShape;

            std::vector<std::pair<double, double>> attempts;
        };

        InstrumentLineShapeEstimationFromKeypointDistance(const std::vector<double>& initialPixelToWavelengthMapping)
            : InstrumentLineShapeEstimation(initialPixelToWavelengthMapping)
        {
        }

        InstrumentLineShapeEstimationFromKeypointDistance(const std::vector<double>& initialPixelToWavelengthMapping, const novac::CCrossSectionData& initialLineShape)
            : InstrumentLineShapeEstimation(initialPixelToWavelengthMapping, initialLineShape)
        {
        }

        /** Creates a rough estimation of the instrument line shape as the Gaussian line shape which best fits to the measured spectrum.
             If this->HasInitialLineShape() is true, then the initial estimation will be used.
             This method has the advantage that the wavelength calibration does not have to be very accurate, but does take quite a few iterations to succeed.
             @param measuredSpectrum A measured sky spectrum containing Fraunhofer lines.
             @param estimatedLineShape Will on successful return be filled with the estimated line shape.
             @param gaussianWidth Will on successful return be filled with the estimated FWHM of the line shape.
             @return A structure showing how the result was achieved. */
        LineShapeEstimationState EstimateInstrumentLineShape(IFraunhoferSpectrumGenerator& fraunhoferSpectrumGen, const CSpectrum& measuredSpectrum, novac::CCrossSectionData& estimatedLineShape, double& fwhm);

    private:

        double GetMedianKeypointDistanceFromSpectrum(const CSpectrum& spectrum, const IndexRange& pixelRange, const std::string& spectrumName) const;

        /** Returns the first and the last index value where the spectrum is consistently above the provided threshold */
        static IndexRange Threshold(const std::vector<double>& spectrum, double threshold);

    };


    /** This is a helper class for estimating the instrument line shape of an instrument
        using a measured spectrum with a well known pixel-to-wavelength calibration by
        pseudo-absorbers representing the error in measured instrument line shape into a DOAS evaluation. */
    class InstrumentLineshapeEstimationFromDoas : public InstrumentLineShapeEstimation
    {
    public:

        InstrumentLineshapeEstimationFromDoas(const std::vector<double>& initialPixelToWavelengthMapping, const novac::CCrossSectionData& initialLineShape, bool addDebugOutput = false);

        InstrumentLineshapeEstimationFromDoas(const std::vector<double>& initialPixelToWavelengthMapping, const novac::SuperGaussianLineShape& initialLineShape, bool addDebugOutput = false);

        struct LineShapeEstimationAttempt
        {
            SuperGaussianLineShape lineShape;

            double error = 0.0;

            double shift = 0.0;
        };

        struct LineShapeEstimationSettings
        {
            /** The first pixel (inclusive) in the range which will be used to estimate the line shape.
                Must be smaller than endPixel. */
            size_t startPixel = 0;

            /** The last pixel after the range which will be used to estimate the line shape.
                Must be smaller than the length of the measured spectrum. */
            size_t endPixel = 0;
        };

        /** A helper structure to provide the result from the estimation as well as the
             measurements of the goodness-of-fit. */
        struct LineShapeEstimationResult
        {
            LineShapeEstimationAttempt result;

            std::vector<LineShapeEstimationAttempt> attempts;
        };

        /** Estimates the instrument line shape by fitting a Super Gaussian to the measured spectrum
         *   @return The fitted line shape
         *   @throw std::invalid_argument if the initial line shape hasn't been provided */
        LineShapeEstimationResult EstimateInstrumentLineShape(
            const CSpectrum& measuredSpectrum,
            const LineShapeEstimationSettings& settings,
            IFraunhoferSpectrumGenerator& fraunhoferSpectrumGen,
            ICrossSectionSpectrumGenerator* ozoneSpectrumGen = nullptr);

    private:

        /** Set to true to enable debugging output to stdout */
        bool debugOutput = false;

        struct LineShapeUpdate
        {
            double currentError = 0.0;    // This is the chi2 of the DOAS fit without the parameter adjustment.
            double residualSize = 0.0;    // This is the chi2 of the DOAS fit with the parameter adjustment.
            double shift = 0.0;           // This is the shift of the DOAS fit with the parameter adjustment.
            std::vector<double> parameterDelta; // The retrieved parameter adjustment.
        };

        /** The initial, starting guess for the instrument line shape.
             This contains the parameterized line shape, as apart from 'initialLineShapeEstimation'
             which contains the sampled data. */
        SuperGaussianLineShape initialLineShapeFunction;

        /** Attempts to calculate the gradient of the parameters of the instrument line shape fit
             by calling CalculateGradientAndCurrentError. */
        LineShapeUpdate GetGradient(
            IFraunhoferSpectrumGenerator& fraunhoferSpectrumGen,
            ICrossSectionSpectrumGenerator* ozoneSpectrumGen,
            const CSpectrum& measuredSpectrum,
            const SuperGaussianLineShape& currentLineShape,
            const LineShapeEstimationSettings& settings,
            bool& allowSpectrumShift);

        /** Calculates the gradient of the parameters of the instrument line shape using a DOAS fit.
             At the same time an currentError measure at the current location is calculated (by reusing the same objects).
             @throws DoasFitException if the fit fails. */
        LineShapeUpdate CalculateGradientAndCurrentError(
            IFraunhoferSpectrumGenerator& fraunhoferSpectrumGen,
            ICrossSectionSpectrumGenerator* ozoneSpectrumGen,
            const CSpectrum& measuredSpectrum,
            const SuperGaussianLineShape& currentLineShape,
            const LineShapeEstimationSettings& settings,
            bool allowShift = true);

        /** Returns false if the provided parameterDelta and stepSize will result in invalid parameter setttings (typically negative values). */
        static bool UpdatedParametersAreValid(const SuperGaussianLineShape& currentLineShape, const std::vector<double>& parameterDelta, double stepSize);

    };

    /** Estimates the Full Width at Half Maximum of a given lineshape function. */
    double GetFwhm(const novac::CCrossSectionData& lineshape);

    /** Estimates the Full Width at Half Maximum of a given lineshape function. */
    double GetFwhm(const std::vector<double>& lineshapeWavelength, const std::vector<double>& lineShapeIntensity);

}