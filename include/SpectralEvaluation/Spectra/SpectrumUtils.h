#pragma once

#include <vector>
#include <stddef.h>

namespace novac
{

class CSpectrum;

enum class SpectrumDataPointType
{
    Peak = 0,   // Represents a point of a local maximum in the spectrum.
    Valley = 1, // Represents a point of a local minimum in the spectrum.
    UnresolvedPeak = 2 // Represents a not fully resolved emission line in the spectrum.
};

// Simple structure used to represent a point in a spectrum
struct SpectrumDataPoint
{
    // The spectrometer pixel where this point is found
    double pixel = 0.0;

    // The wavelength (in nm) where this point is found, set to zero if not known
    double wavelength = 0.0;

    // The intensity value of the spectrum at this point
    double intensity = 0.0;

    // A type enumerator, may be used to classify points
    SpectrumDataPointType type = SpectrumDataPointType::Peak;

    // An associated left point
    double leftPixel = 0.0;

    // An associated right point 
    double rightPixel = 0.0;

    // True if this is a flat-top region (indication of saturation)
    bool flatTop = false;
};

/**
 * @brief Locates all significant peaks in the provided spectrum and returns the result in the provided vector
 * @param spectrum The spectrum in which peaks should be found,
 *  if this has a wavelength calibration then the resulting points will have a wavelength filled in.
 * @param minimumIntensity Only peaks with a maximum intensity above this value will be returned.
 * @param result Will on return be filled with the found peaks.
 * The resulting points will have their leftPixel and rightPixel set to the points where the peak is judged to start. */
void FindPeaks(const CSpectrum& spectrum, double minimumIntensity, std::vector<SpectrumDataPoint>& result);

/**
 * @brief Locates all significant valleys in the provided spectrum and returns the result in the provided vector
 * @param spectrum The spectrum in which valleys should be found,
 *  if this has a wavelength calibration then the resulting points will have a wavelength filled in.
 * @param minimumIntensity Only valleys with a maximum intensity above this value will be returned.
 * @param result Will on return be filled with the found valleys.
 * The resulting points will have their leftPixel and rightPixel set to the points where the valley is judged to start. */
void FindValleys(const CSpectrum& spectrum, double minimumIntensity, std::vector<SpectrumDataPoint>& result);

/**
 * @brief Locates all significant peaks _and_ valleys in the provided spectrum and returns the result in the provided vector
 * @param spectrum The spectrum in which keypoints should be found,
 *  if this has a wavelength calibration then the resulting points will have a wavelength filled in.
 * @param minimumIntensity Only keypoints with an intensity above this value will be returned.
 * @param result Will on return be filled with the found keypoints.
 * The valleys found will have type set to -1 and peaks found will have type set to +1
 * The result is sorted with increasing pixel values. */
void FindKeypointsInSpectrum(const CSpectrum& spectrum, double minimumIntensity, std::vector<SpectrumDataPoint>& result);

/**
 * @brief Locates all significant peaks in the provided spectrum, assuming that the provided spectrum is an emission spectrum
 *  containing a number of well separated emission lines.
 * @param spectrum The spectrum in which emission lines should be found,
 *  if this has a wavelength calibration then the resulting points will have a wavelength filled in.
 * @param result Will on return be filled with the found peaks.
 * The resulting points will have their leftPixel and rightPixel set to the points where the peak is judged to start. */
void FindEmissionLines(const CSpectrum& spectrum, std::vector<SpectrumDataPoint>& result, bool includeNotClearlyResolvedLines = true);

/**
    @brief Filters the provided input vector by the SpectrumDataPoint::type and returns
        a copy of the input containing only the input points which have the provided type. */
std::vector<SpectrumDataPoint> FilterByType(const std::vector<SpectrumDataPoint>& input, SpectrumDataPointType typeToFind);

/**
 * @brief Calculates a first order derivative on the provided data array wrt index value using finite difference.
 * @param data The data array (spectrum) to calculate the derivative of.
 * @param dataLength The length of data.
 * @param result Will on successful result be filled with the calculated first order derivative.
 *       This will be reallocated to dataLength long and the first and last elements will be set to zero
 * @returns true if all is ok */
bool Derivative(const double* data, size_t dataLength, std::vector<double>& result);

/**
 * @brief Calculates a first or second order derivative on the provided data array wrt index value using finite difference.
 * @param data The data array (spectrum) to calculate the derivative of.
 * @param dataLength The length of data.
 * @param result Will on successful result be filled with the calculated first order derivative.
 *       This will be reallocated to dataLength long and the first and last elements will be set to zero
 * @returns true if all is ok */
bool Derivative(const double* data, size_t dataLength, int order, std::vector<double>& result);

/**
 * @brief Locates significant points in the provided spectrum with the property that the points (pixel[ii], intensity[ii])
 *  lies on the spectrum and all points in the spectrum lies either at or below the lines between the found points.
 *  This is the envelope of the spectrum, or the complex-hull if you will.
 * @param spectrum The spectrum in which points should be found.
 * @param pixel Will on return be filled with the pixel-values of the found points.
 * @param intensity Will on return be filled with the intensity of the found points (in the same units as the spectrum). s*/
bool GetEnvelope(const CSpectrum& spectrum, std::vector<double>& pixel, std::vector<double>& intensity);

/** Searches for the pixel in the provided pixel-to-wavelength mapping where the provided
*   wavelength can be found. Will round to the nearest pixel.
*   If wavelength is smaller than the smalles value in pixelToWavelengthMapping then zero is returned.
*   If wavelength is larger than the largest value in pixelToWavelengthMapping, then pixelToWavelengthMapping.size() - 1 is returned. */
size_t WavelengthToPixel(const std::vector<double>& pixelToWavelengthMapping, double wavelength);

/** Searches for the pixel in the provided pixel-to-wavelength mapping where the provided
*   wavelength can be found and returns the fractional pixel-index where this is found.
*   If wavelength is smaller than the smalles value in pixelToWavelengthMapping then zero is returned.
*   If wavelength is larger than the largest value in pixelToWavelengthMapping, then pixelToWavelengthMapping.size() - 1 is returned. */
double WavelengthToFractionalPixel(const std::vector<double>& pixelToWavelengthMapping, double wavelength);

}
