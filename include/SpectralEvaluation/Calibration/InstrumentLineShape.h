#pragma once
#include <vector>

// ---------------------------------------------------------------------------------------------------------------
// -- This header contains methods used to extract and characterize the instrument line shape of a spectrometer --
// ---------------------------------------------------------------------------------------------------------------

namespace novac
{

class CSpectrum;

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
// The measured spectrum needs to be dark corrected and corrected such that any offset has been removed
ILF_RETURN_CODE FitInstrumentLineShape(const CSpectrum& mercuryLine, GaussianLineShape& result);
ILF_RETURN_CODE FitInstrumentLineShape(const CSpectrum& mercuryLine, AsymmetricGaussianLineShape& result);
ILF_RETURN_CODE FitInstrumentLineShape(const CSpectrum& mercuryLine, SuperGaussianLineShape& result);

/**  Calculates the value of the provided line shape on the provided x-axis grid
    The additional parameters required to calculate the value of the line shape are provided as additional parameters
    @param center The x-axis value around which the line shape should be centered
    @param amplitude The maximum value of the line shape (above the baseline)
    @param baseline This value is added to each output value  */
std::vector<double> SampleInstrumentLineShape(const GaussianLineShape& lineShape, const std::vector<double>& x, double center, double amplitude, double baseline = 0.0);
std::vector<double> SampleInstrumentLineShape(const AsymmetricGaussianLineShape& lineShape, const std::vector<double>& x, double center, double amplitude, double baseline = 0.0);
std::vector<double> SampleInstrumentLineShape(const SuperGaussianLineShape& lineShape, const std::vector<double>& x, double center, double amplitude, double baseline = 0.0);


enum class MercuryPeak
{
    Hg302nm = 0
};

// This takes a full measured (but dark corrected) mercury spectrum with a reasonably correct 
// wavelength calibration, extracts the given mercury peak and fits a Gaussian line shape to this.
ILF_RETURN_CODE GetInstrumentLineShape(const CSpectrum& fullMercurySpectrum, MercuryPeak peakToUse, GaussianLineShape& result);

}
