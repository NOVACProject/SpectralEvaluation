#include "catch.hpp"
#include <SpectralEvaluation/Spectra/InstrumentLineShape.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>

std::vector<double> CreatePixelToWavelengthMapping(double start, double stop, int size = 2048); // located elsewhere
std::vector<double> CreateGaussian(double sigma, const std::vector<double>& x); // located elsewhere


// -------- Fitting symmetrical Gaussians --------
TEST_CASE("FitInstrumentLineShape (Symmetric Gaussian): Simple Gaussian input returns correct fitted function", "[InstrumentLineShape]")
{
    const double sigma = 0.2; // [nm]
    CSpectrum idealGaussian;
    idealGaussian.m_wavelength = CreatePixelToWavelengthMapping(-2.0, +2.0, 43); // 4nm range with 43 sample pointss
    std::vector<double> spec = CreateGaussian(sigma, idealGaussian.m_wavelength);
    memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));
    idealGaussian.m_length = (long)spec.size();

    Evaluation::GaussianLineShape result;
    auto ret = Evaluation::FitInstrumentLineShape(idealGaussian, result);

    REQUIRE(ret == Evaluation::ILF_RETURN_CODE::SUCCESS);
    REQUIRE(fabs(result.sigma - sigma) < 0.005);
}
