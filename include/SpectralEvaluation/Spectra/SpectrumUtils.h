#pragma once

#include <vector>

namespace novac
{

class CSpectrum;

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
    int type = 0;
};

/**
 * @brief Locates all significant peaks in the provided spectrum and returns the result in the provided vector
 * @param spectrum The spectrum in which peaks should be found,
 *  if this has a wavelength calibration then the resulting points will have a wavelength filled in.
 * @param minimumIntensity Only peaks with a maximum intensity above this value will be returned.
 * @param result Will on return be filled with the found peaks.
*/
void FindPeaks(const CSpectrum& spectrum, double minimumIntensity, std::vector<SpectrumDataPoint>& result);

/**
 * @brief Locates all significant valleys in the provided spectrum and returns the result in the provided vector
 * @param spectrum The spectrum in which valleys should be found,
 *  if this has a wavelength calibration then the resulting points will have a wavelength filled in.
 * @param minimumIntensity Only valleys with a maximum intensity above this value will be returned.
 * @param result Will on return be filled with the found valleys.
*/
void FindValleys(const CSpectrum& spectrum, double minimumIntensity, std::vector<SpectrumDataPoint>& result);

/**
 * @brief Locates all significant peaks _and_ valleys in the provided spectrum and returns the result in the provided vector
 * @param spectrum The spectrum in which keypoints should be found,
 *  if this has a wavelength calibration then the resulting points will have a wavelength filled in.
 * @param minimumIntensity Only keypoints with an intensity above this value will be returned.
 * @param result Will on return be filled with the found keypoints.
*/
void FindKeypointsInSpectrum(const CSpectrum& spectrum, double minimumIntensity, std::vector<SpectrumDataPoint>& result);

/**
 * @brief Calculates a first order derivative on the provided data array wrt index value using finite difference.
 * @param data The data array (spectrum) to calculate the derivative of.
 * @param dataLength The length of data.
 * @param result Will on successful result be filled with the calculated first order derivative.
 *       This will be reallocated to dataLength long and the first and last elements will be set to zero
 * @returns true if all is ok
*/
bool Derivative(const double* data, size_t dataLength, std::vector<double>& result);

/**
 * @brief Calculates a first or second order derivative on the provided data array wrt index value using finite difference.
 * @param data The data array (spectrum) to calculate the derivative of.
 * @param dataLength The length of data.
 * @param result Will on successful result be filled with the calculated first order derivative.
 *       This will be reallocated to dataLength long and the first and last elements will be set to zero
 * @returns true if all is ok
*/
bool Derivative(const double* data, size_t dataLength, int order, std::vector<double>& result);

/**
 * @brief Locates significant points in the provided spectrum with the property that the points (pixel[ii], intensity[ii])
 *  lies on the spectrum and all points in the spectrum lies either at or below the lines between the found points.
 *  This is the envelope of the spectrum, or the complex-hull if you will.
 * @param spectrum The spectrum in which points should be found.
 * @param pixel Will on return be filled with the pixel-values of the found points.
 * @param intensity Will on return be filled with the intensity of the found points (in the same units as the spectrum).
*/
bool GetEnvelope(const CSpectrum& spectrum, std::vector<double>& pixel, std::vector<double>& intensity);

}
