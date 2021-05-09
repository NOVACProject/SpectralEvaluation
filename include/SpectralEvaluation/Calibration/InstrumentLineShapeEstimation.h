#pragma once

#include <vector>
#include <string>

namespace novac
{
class CSpectrum;
class CCrossSectionData;

/// <summary>
/// This is a helper class for estimating the instrument line shape of an instrument 
/// using a measured spectrum with a (reasonably well known) pixel-to-wavelength calibration.
/// 
/// </summary>
class InstrumentLineShapeEstimation
{
public:
    InstrumentLineShapeEstimation(const std::string& highResolutionSolarAtlas, const std::vector<double>& pixelToWavelengthMapping)
        : highResolutionSolarAtlas(highResolutionSolarAtlas), pixelToWavelengthMapping(pixelToWavelengthMapping)
    {
    }

    /// <summary>
    /// Creates a rough estimation of the instrument line shape as the Gaussian line shape which best fits to the measured spectrum.
    /// </summary>
    /// <param name="measuredSpectrum">A measured sky spectrum containing Fraunhofer lines</param>
    /// <param name="estimatedLineShape">Will on successful return be filled with the estimated line shape</param>
    /// <param name="gaussianWidth">Will on successful return be filled with the estimated FWHM of the line shape</param>
    void EstimateInstrumentLineShape(const CSpectrum& measuredSpectrum, novac::CCrossSectionData& estimatedLineShape, double& fwhm);

    /// <summary>
    /// Sets the pixel to wavelength mapping for the spectrometer,
    /// this is the inital state before doing the calibration.
    /// </summary>
    void UpdateInitialCalibration(const std::vector<double>& newPixelToWavelengthMapping)
    {
        this->pixelToWavelengthMapping = newPixelToWavelengthMapping;
    }

private:
    const std::string highResolutionSolarAtlas;

    std::vector<double> pixelToWavelengthMapping;

    double GetMedianKeypointDistanceFromSpectrum(const CSpectrum& spectrum) const;
};
}