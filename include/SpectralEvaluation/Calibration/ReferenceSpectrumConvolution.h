#pragma once

#include <vector>
#include <string>

// ---------------------------------------------------------------------------------------------------------------
// -------------- This header contains methods used to prepare reference spectra for the evaluation --------------
// ---------------------------------------------------------------------------------------------------------------

namespace novac
{
enum class WavelengthConversion
{
    None,
    VacuumToAir
};

enum class ConvolutionMethod
{
    Direct,
    Fft
};

class CCrossSectionData;

/** Performs a convolution of the high resolution reference function with the given slf (slit function, the convolution core)
    and resamples the result to the given pixelToWavelengthMapping.
    This expects the slf to be shifted to have the center in the middle of the vector. */
bool ConvolveReference(
    const std::string& pixelToWavelengthMapping,
    const std::string& slf,
    const std::string& highResReference,
    CCrossSectionData& result,
    WavelengthConversion conversion = WavelengthConversion::None);

/** Performs a convolution of the high resolution reference function with the given gaussian slf (slit function, the convolution core)
    and resamples the result to the given pixelToWavelengthMapping.
    This expects the slf to be shifted to have the center in the middle of the vector. */
bool ConvolveReferenceGaussian(const std::string& pixelToWavelengthMapping, double gaussianSigma, const std::string& highResReference, CCrossSectionData& result, WavelengthConversion conversion = WavelengthConversion::None);

/** Performs a convolution of the high resolution reference function with the given slf (slit function, the convolution core)
    and resamples the result to the given pixelToWavelengthMapping.
    This expects the slf to be shifted to have the center in the middle of the vector.
    If fwhmOfInstrumentLineShape is provided (larger than zero) than this will be used to determine the convolution grid,
        otherwise this will be determined from the shape of the slf.
    @throws std::invalid_argument if the slf or highResReference does not have a valid pixel-to-wavelength mapping, or the calculation of the fwhm failed. */
void ConvolveReference(
    const std::vector<double>& pixelToWavelengthMapping,
    const CCrossSectionData& slf,
    const CCrossSectionData& highResReference,
    std::vector<double>& result,
    WavelengthConversion conversion = WavelengthConversion::None,
    ConvolutionMethod method = ConvolutionMethod::Direct,
    double fwhmOfInstrumentLineShape = 0.0,
    bool normalizeSlf = true);

/** Performs a convolution of the high resolution reference function with the given slf (slit function, the convolution core).
    The result will be sampled on the same wavelength grid as the highResReference.
    This expects the slf to be shifted to have the center in the middle of the vector. */
bool Convolve(const CCrossSectionData& slf, const CCrossSectionData& highResReference, std::vector<double>& result, WavelengthConversion conversion = WavelengthConversion::None);

/** Takes a simple slit-function as input and returns the Full Width at Half Maximum (FWHM), in the x-axis unit of the SLF. */
double CalculateFhwm(const CCrossSectionData& slf);

}
