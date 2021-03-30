#pragma once

#include <vector>

// ---------------------------------------------------------------------
// -------------------- Elementary function fitting --------------------
// ---------------------------------------------------------------------

namespace MathFit
{
class CGaussFunction;
}

namespace novac
{
// TODO: Rename to something showing that this is the function-fitting-return-code
enum class ILF_RETURN_CODE
{
    SUCCESS = 0,
    FIT_FAILURE,
    EMPTY_INPUT,
    MISSING_WAVELENGTH_CALIBRATION
};

/// <summary>
/// Creates a basic initial estimate of a gaussian function which could fit to the provided data.
/// This does not do the actual fit but uses heuristics to get an initial estimation, making
/// the Gaussian fit itself more stable.
/// </summary>
void CreateInitialEstimate(const std::vector<double>& x, const std::vector<double>& y, MathFit::CGaussFunction& result);

/// <summary>
/// Fits a Gaussian function to the provided data.
/// </summary>
ILF_RETURN_CODE FitGaussian(std::vector<double>& x, std::vector<double>& y, MathFit::CGaussFunction& gaussian);

/// <summary>
/// Fits a Gaussian function to the provided data.
/// The x-axis data will be index.
/// </summary>
ILF_RETURN_CODE FitGaussian(std::vector<double>& y, MathFit::CGaussFunction& gaussian);

}