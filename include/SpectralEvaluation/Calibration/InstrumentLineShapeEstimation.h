#pragma once

#include <memory>
#include <vector>
#include <string>
#include <SpectralEvaluation/Calibration/InstrumentLineShape.h>

namespace novac
{
class CSpectrum;
class CCrossSectionData;
class IFraunhoferSpectrumGenerator;
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

    const char* what() const noexcept override final { return m_msg; }
};

/// <summary>
/// Abstract base class for the different instrument line shape estimators.
/// </summary>
class InstrumentLineShapeEstimation
{
public:

    /// <summary>
    /// Sets the pixel to wavelength mapping for the spectrometer,
    /// this is the inital state before doing the calibration.
    /// </summary>
    void UpdateInitialCalibration(const std::vector<double>& newPixelToWavelengthMapping)
    {
        this->pixelToWavelengthMapping = newPixelToWavelengthMapping;
    }

    /// <summary>
    /// Sets initial estimation of the instrument line shape. Used as a starting point for the estimation routine.
    /// </summary>
    void UpdateInitialLineShape(const novac::CCrossSectionData& newInitialLineShape)
    {
        this->initialLineShapeEstimation = std::make_unique<novac::CCrossSectionData>(newInitialLineShape);
    }

    bool HasInitialLineShape() const;

protected:

    InstrumentLineShapeEstimation(const std::vector<double>& initialPixelToWavelengthMapping)
        : pixelToWavelengthMapping(initialPixelToWavelengthMapping)
    {
    }

    InstrumentLineShapeEstimation(const std::vector<double>& initialPixelToWavelengthMapping, const novac::CCrossSectionData& initialLineShape)
        : pixelToWavelengthMapping(initialPixelToWavelengthMapping)
    {
        initialLineShapeEstimation = std::make_unique<novac::CCrossSectionData>(initialLineShape);
    }

    /// <summary>
    /// The assumed pixel-to-wavelength mapping for the device.
    /// </summary>
    std::vector<double> pixelToWavelengthMapping;

    /// <summary>
    /// The initial, starting guess for the instrument line shape.
    /// </summary>
    std::unique_ptr<novac::CCrossSectionData> initialLineShapeEstimation;
};

/// <summary>
/// This is a helper class for estimating the instrument line shape of an instrument 
/// using a measured spectrum with a (reasonably well known) pixel-to-wavelength calibration
/// by measuring the distance between keypoints in the meaured and synthetic spectra.
/// </summary>
class InstrumentLineShapeEstimationFromKeypointDistance : public InstrumentLineShapeEstimation
{
public:
    /// <summary>
    /// A helper structure to describe the internal state of this estimator
    /// </summary>
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

    /// <summary>
    /// Creates a rough estimation of the instrument line shape as the Gaussian line shape which best fits to the measured spectrum.
    /// If this->HasInitialLineShape() is true, then the initial estimation will be used.
    /// This method has the advantage that the wavelength calibration does not have to be very accurate, but does take quite a few iterations to succeed.
    /// </summary>
    /// <param name="measuredSpectrum">A measured sky spectrum containing Fraunhofer lines</param>
    /// <param name="estimatedLineShape">Will on successful return be filled with the estimated line shape</param>
    /// <param name="gaussianWidth">Will on successful return be filled with the estimated FWHM of the line shape</param>
    /// <returns>A structure showing how the result was achieved</returns>
    LineShapeEstimationState EstimateInstrumentLineShape(IFraunhoferSpectrumGenerator& fraunhoferSpectrumGen, const CSpectrum& measuredSpectrum, novac::CCrossSectionData& estimatedLineShape, double& fwhm);

private:

    double GetMedianKeypointDistanceFromSpectrum(const CSpectrum& spectrum, const IndexRange& pixelRange, const std::string& spectrumName) const;
};


/// <summary>
/// This is a helper class for estimating the instrument line shape of an instrument
/// using a measured spectrum with a well known pixel-to-wavelength calibration by
/// pseudo-absorbers representing the error in measured instrument line shape into a DOAS evaluation.
/// </summary>
class InstrumentLineshapeEstimationFromDoas : public InstrumentLineShapeEstimation
{
public:

    InstrumentLineshapeEstimationFromDoas(const std::vector<double>& initialPixelToWavelengthMapping, const novac::CCrossSectionData& initialLineShape);

    InstrumentLineshapeEstimationFromDoas(const std::vector<double>& initialPixelToWavelengthMapping, const novac::SuperGaussianLineShape& initialLineShape);

    struct LineShapeEstimationAttempt
    {
        SuperGaussianLineShape lineShape;

        double error = 0.0;

        double shift = 0.0;
    };

    struct LineShapeEstimationSettings
    {
        /// <summary>
        /// The first pixel (inclusive) in the range which will be used to estimate the line shape.
        /// Must be smaller than endPixel.
        /// </summary>
        size_t startPixel = 0;

        /// <summary>
        /// The last pixel after the range which will be used to estimate the line shape.
        /// Must be smaller than the length of the measured spectrum.
        /// </summary>
        size_t endPixel = 0;
    };

    /// <summary>
    /// A helper structure to provide the result from the estimation as well as the 
    /// measurements of the goodness-of-fit.
    /// </summary>
    struct LineShapeEstimationResult
    {
        LineShapeEstimationAttempt result;

        std::vector<LineShapeEstimationAttempt> attempts;
    };

    /// <summary>
    /// Estimates the instrument line shape by fitting a Super Gaussian to the measured spectrum
    /// </summary>
    /// <returns>The fitted line shape</returns>
    /// <throws>std::invalid_argument if the initial line shape hasn't been provided.</throws>
    LineShapeEstimationResult EstimateInstrumentLineShape(
        IFraunhoferSpectrumGenerator& fraunhoferSpectrumGen,
        const CSpectrum& measuredSpectrum,
        const LineShapeEstimationSettings& settings);

private:

    struct LineShapeUpdate
    {
        double currentError;    // This is the chi2 of the DOAS fit without the parameter adjustment.
        double residualSize;    // This is the chi2 of the DOAS fit with the parameter adjustment.
        double shift;           // This is the shift of the DOAS fit with the parameter adjustment.
        std::vector<double> parameterDelta; // The retrieved parameter adjustment.
    };

    /// <summary>
    /// The initial, starting guess for the instrument line shape.
    /// This contains the parameterized line shape, as apart from 'initialLineShapeEstimation'
    /// which contains the sampled data.
    /// </summary>
    SuperGaussianLineShape initialLineShapeFunction;

    /// <summary>
    /// Attempts to calculate the gradient of the parameters of the instrument line shape fit
    /// by calling CalculateGradientAndCurrentError.
    /// </summary>
    LineShapeUpdate GetGradient(
        IFraunhoferSpectrumGenerator& fraunhoferSpectrumGen,
        const CSpectrum& measuredSpectrum,
        const SuperGaussianLineShape& currentLineShape,
        const LineShapeEstimationSettings& settings,
        bool& allowSpectrumShift);

    /// <summary>
    /// Calculates the gradient of the parameters of the instrument line shape using a DOAS fit.
    /// At the same time an currentError measure at the current location is calculated (by reusing the same objects).
    /// @throws DoasFitException if the fit fails.
    /// </summary>
    LineShapeUpdate CalculateGradientAndCurrentError(
        IFraunhoferSpectrumGenerator& fraunhoferSpectrumGen,
        const CSpectrum& measuredSpectrum,
        const SuperGaussianLineShape& currentLineShape,
        const LineShapeEstimationSettings& settings,
        bool allowShift = true);

};

/// <summary>
/// Estimates the Full Width at Half Maximum of a given lineshape function.
/// </summary>
double GetFwhm(const novac::CCrossSectionData& lineshape);

/// <summary>
/// Estimates the Full Width at Half Maximum of a given lineshape function.
/// </summary>
double GetFwhm(const std::vector<double>& lineshapeWavelength, const std::vector<double>& lineShapeIntensity);

}