#pragma once
#include <vector>
#include <SpectralEvaluation/Math/FunctionFit.h>

// ---------------------------------------------------------------------------------------------------------------
// -- This header contains methods used to extract and characterize the instrument line shape of a spectrometer --
// ---------------------------------------------------------------------------------------------------------------

namespace novac
{

class CSpectrum;
class CCrossSectionData;

// --------- Possible representations of instrument line shapes ---------
// Symmetric gaussian line shape: exp(-x^2/(2 * sigma^2))
struct GaussianLineShape
{
    // The width parameter
    double sigma = 0.0;

    // The center location
    double center = 0.0;

    // The full width at half maximum.
    double Fwhm() const { return sigma * 2.35482004; }
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

// Symmetric super-gaussian line shape: exp(- abs(x/w)^k)
// A regular gaussian is a special case of this, with k=2.0 and w=sigma*sqrt(2)
struct SuperGaussianLineShape
{
    // The width parameter.
    double w = 0.0;

    // The exponent (k = 2.0 equals a 'regular' Gaussian)
    double k = 2.0;

    // The center location
    double center = 0.0;

    // The full width at half maximum.
    double Fwhm() const;
};

// Fits a symmetrical Gaussian line to an extract of a mercury spectrum containing only one (full) mercury line.
// The measured spectrum needs to have a wavelength calibration
// The measured spectrum needs to be dark corrected and corrected such that any offset has been removed
FUNCTION_FIT_RETURN_CODE FitInstrumentLineShape(const CSpectrum& mercuryLine, GaussianLineShape& result);
FUNCTION_FIT_RETURN_CODE FitInstrumentLineShape(const CSpectrum& mercuryLine, AsymmetricGaussianLineShape& result);
FUNCTION_FIT_RETURN_CODE FitInstrumentLineShape(const CSpectrum& mercuryLine, SuperGaussianLineShape& result);

FUNCTION_FIT_RETURN_CODE FitInstrumentLineShape(const CCrossSectionData& mercuryLine, SuperGaussianLineShape& result);

/**  Calculates the value of the provided line shape on the provided x-axis grid
    The additional parameters required to calculate the value of the line shape are provided as additional parameters
    @param center The x-axis value around which the line shape should be centered
    @param amplitude The maximum value of the line shape (above the baseline)
    @param baseline This value is added to each output value  */
std::vector<double> SampleInstrumentLineShape(const GaussianLineShape& lineShape, const std::vector<double>& x, double center, double amplitude, double baseline = 0.0);
std::vector<double> SampleInstrumentLineShape(const AsymmetricGaussianLineShape& lineShape, const std::vector<double>& x, double center, double amplitude, double baseline = 0.0);
std::vector<double> SampleInstrumentLineShape(const SuperGaussianLineShape& lineShape, const std::vector<double>& x, double center, double amplitude, double baseline = 0.0);

/** Calculates the value of the provided line shape on an auto determined x-axis grid. 
    The returned line shape will be centered on zero and have a normalized amplitude. */
CCrossSectionData SampleInstrumentLineShape(const SuperGaussianLineShape& lineShape);

/** Calculates the partial derivative with respect to the 'sigma' parameter of the provided Gaussian line shape using finite difference. */
std::vector<double> PartialDerivative(const GaussianLineShape& lineShape, const std::vector<double>& x);

/** Calculates the partial derivative with respect to the given parameter of the provided super Gaussian line shape using finite difference.
    @param parameter 0 corresponds to partial derivative with respect to 'w'. 1 corresponds to partial derivative wrt 'k'. */
std::vector<double> PartialDerivative(const SuperGaussianLineShape& lineShape, const std::vector<double>& x, int parameter);

/** Retrieves the value of the given parameter (w == 0 and k == 1) */
double GetParameterValue(const SuperGaussianLineShape& lineShape, int parameterIdx);

/** Sets the value of the given parameter (w == 0 and k == 1) */
void SetParameterValue(SuperGaussianLineShape& lineShape, int parameterIdx, double value);

}
