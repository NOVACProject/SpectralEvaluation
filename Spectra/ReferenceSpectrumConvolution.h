#pragma once

#include <vector>

// ---------------------------------------------------------------------------------------------------------------
// -------------- This header contains methods used to prepare reference spectra for the evaluation --------------
// ---------------------------------------------------------------------------------------------------------------

/** A simple struct used to represent a spectrum sampled on a specific wavelength grid.
    The vectors wavelength and value are expected to be of equal length. */
struct SimpleSpectrum
{
    SimpleSpectrum() { }

    explicit SimpleSpectrum(size_t length)
        : wavelength(length, 0.0), value(length, 0.0)
    {
    }

    std::vector<double> wavelength;
    std::vector<double> value;
};

/** Performs a convolution of the high resolution reference function with the given slf (slit function, the convolution core)
    and resamples the result to the given pixelToWavelengthMapping. 
    This expects the slf to be normalized (amplitude range 0->1) and shifted to have the center in the middle of the vector. */
bool ConvolveReference(const std::vector<double>& pixelToWavelengthMapping, const std::vector<double>& slf, const std::vector<double>& highResReference, std::vector<double>& result);

/** Performs a convolution of the high resolution reference function with the given slf (slit function, the convolution core).
    The result will be sampled on the same wavelength grid as the highResReference. 
    This expects the slf to be normalized (amplitude range 0->1) and shifted to have the center in the middle of the vector. */
bool Convolve(const SimpleSpectrum& slf, const SimpleSpectrum& highResReference, std::vector<double>& result);
