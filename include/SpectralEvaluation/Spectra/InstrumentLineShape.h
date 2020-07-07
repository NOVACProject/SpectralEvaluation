#pragma once

// ---------------------------------------------------------------------------------------------------------------
// -- This header contains methods used to extract and characterize the instrument line shape of a spectrometer --
// ---------------------------------------------------------------------------------------------------------------

class CSpectrum;

namespace Evaluation
{

// --------- Possible representations of instrument line shapes ---------
// Symmetric gaussian line shape: exp(-x^2/(2 * sigma^2))
struct GaussianLineShape
{
    // The width parameter
    double sigma = 0.0;

    double Fwhm() const { return sigma * 2.35482; }
};

// Asymmetric gaussian line shape, consisting of a left and a right half with different widths.
//  The 'left' is used for x < center and 'right' used for x >= center
struct AsymmetricGaussianLineShape
{
    // The width parameter of the left gaussian
    double sigmaLeft = 0.0;

    // The width parameter of the right gaussian
    double sigmaRight = 0.0;
};

// Symmetric super-gaussian line shape: exp(-0.5 * [x/sigma]^P)
struct SuperGaussianLineShape
{
    // The width parameter
    double sigma = 0.0;

    // The exponent (P = 2.0 equals a 'regular' Gaussian)
    double P = 2.0;
};

enum class ILF_RETURN_CODE
{
    SUCCESS = 0,
    FIT_FAILURE,
    EMPTY_INPUT,
    MISSING_WAVELENGTH_CALIBRATION
};

// Fits a symmetrical Gaussian line to an extract of a mercury spectrum containing only one (full) mercury line.
// The measured spectrum needs to have a wavelength calibration
ILF_RETURN_CODE FitInstrumentLineShape(const CSpectrum& mercuryLine, GaussianLineShape& result);
ILF_RETURN_CODE FitInstrumentLineShape(const CSpectrum& mercuryLine, AsymmetricGaussianLineShape& result);
ILF_RETURN_CODE FitInstrumentLineShape(const CSpectrum& mercuryLine, SuperGaussianLineShape& result);

}
