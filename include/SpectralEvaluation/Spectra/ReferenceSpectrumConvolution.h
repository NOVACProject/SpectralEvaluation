#pragma once

#include <vector>
#include <string>

// ---------------------------------------------------------------------------------------------------------------
// -------------- This header contains methods used to prepare reference spectra for the evaluation --------------
// ---------------------------------------------------------------------------------------------------------------

namespace Evaluation
{
    class CCrossSectionData;

    /** Performs a convolution of the high resolution reference function with the given slf (slit function, the convolution core)
        and resamples the result to the given pixelToWavelengthMapping.
        This expects the slf to be shifted to have the center in the middle of the vector. */
    bool ConvolveReference(const std::string& pixelToWavelengthMapping, const std::string& slf, const std::string& highResReference, CCrossSectionData& result);

    /** Performs a convolution of the high resolution reference function with the given gaussian slf (slit function, the convolution core)
        and resamples the result to the given pixelToWavelengthMapping.
        This expects the slf to be shifted to have the center in the middle of the vector. */
    bool ConvolveReferenceGaussian(const std::string& pixelToWavelengthMapping, double gaussianSigma, const std::string& highResReference, CCrossSectionData& result);

    /** Performs a convolution of the high resolution reference function with the given slf (slit function, the convolution core)
        and resamples the result to the given pixelToWavelengthMapping.
        This expects the slf to be shifted to have the center in the middle of the vector. */
    bool ConvolveReference(const std::vector<double>& pixelToWavelengthMapping, const CCrossSectionData& slf, const CCrossSectionData& highResReference, std::vector<double>& result);

    /** Performs a convolution of the high resolution reference function with the given slf (slit function, the convolution core).
        The result will be sampled on the same wavelength grid as the highResReference.
        This expects the slf to be shifted to have the center in the middle of the vector. */
    bool Convolve(const CCrossSectionData& slf, const CCrossSectionData& highResReference, std::vector<double>& result);

    /** Takes a simple slit-function as input and returns the Full Width at Half Maximum (FWHM), in the x-axis unit of the SLF. */
    double CalculateFhwm(const CCrossSectionData& slf);

}
