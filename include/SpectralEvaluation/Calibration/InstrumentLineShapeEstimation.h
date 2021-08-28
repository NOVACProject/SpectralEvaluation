#pragma once

#include <memory>
#include <vector>
#include <string>

namespace novac
{
class CSpectrum;
class CCrossSectionData;
class IFraunhoferSpectrumGenerator;

/// <summary>
/// This is a helper class for estimating the instrument line shape of an instrument 
/// using a measured spectrum with a (reasonably well known) pixel-to-wavelength calibration.
/// 
/// </summary>
class InstrumentLineShapeEstimation
{
public:

    /// <summary>
    /// A helper structure to describe the internal state of this estimator
    /// </summary>
    struct LineShapeEstimationState
    {
        double medianPixelDistanceInMeas = 0.0;

        std::vector<std::pair<double, double>> attempts;
    };

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
    /// Creates a rough estimation of the instrument line shape as the Gaussian line shape which best fits to the measured spectrum.
    /// If this->HasInitialLineShape() is true, then the initial estimation will be used.
    /// This method has the advantage that the wavelength calibration does not have to be very accurate, but does take quite a few iterations to succeed.
    /// </summary>
    /// <param name="measuredSpectrum">A measured sky spectrum containing Fraunhofer lines</param>
    /// <param name="estimatedLineShape">Will on successful return be filled with the estimated line shape</param>
    /// <param name="gaussianWidth">Will on successful return be filled with the estimated FWHM of the line shape</param>
    /// <returns>A structure showing how the result was achieved</returns>
    LineShapeEstimationState EstimateInstrumentLineShape(IFraunhoferSpectrumGenerator& fraunhoferSpectrumGen, const CSpectrum& measuredSpectrum, novac::CCrossSectionData& estimatedLineShape, double& fwhm);

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

private:
    /// <summary>
    /// The assumed pixel-to-wavelength mapping for the device.
    /// </summary>
    std::vector<double> pixelToWavelengthMapping;

    /// <summary>
    /// The initial, starting guess for the instrument line shape.
    /// </summary>
    std::unique_ptr<novac::CCrossSectionData> initialLineShapeEstimation;

    /// <summary>
    /// The first pixel to include when checking the properties of the spectrum. 
    /// Often do the signal in the spectra decline at short wavelengths and this is a means to disregard points with low intensity.
    /// TODO: Settable or adaptable to input
    /// </summary>
    const size_t measuredPixelStart = 1000;

    /// <summary>
    /// The last pixel to include when checking the properties of the spectrum. 
    /// Often do the signal in the spectra decline at long wavelengths and this is a means to disregard points with low intensity.
    /// This must be larger than measuredPixelStart.
    /// TODO: Settable or adaptable to input
    /// </summary>
    const size_t measuredPixelStop = 4095;

    double GetMedianKeypointDistanceFromSpectrum(const CSpectrum& spectrum, const std::string& spectrumName) const;
};

/// <summary>
/// Estimates the Full Width at Half Maximum of a given lineshape function.
/// </summary>
double GetFwhm(const novac::CCrossSectionData& lineshape);

}