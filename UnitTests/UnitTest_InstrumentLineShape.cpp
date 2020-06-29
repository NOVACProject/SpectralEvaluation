#include "catch.hpp"
#include <SpectralEvaluation/Spectra/InstrumentLineShape.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>

std::vector<double> CreatePixelToWavelengthMapping(double start, double stop, int size = 2048); // located elsewhere
std::vector<double> CreateGaussian(double sigma, const std::vector<double>& x); // located elsewhere

// creates an asymetric gaussian with width 'sigmaNegative' for x < 0 and 'sigmaPositive' for x > 0
std::vector<double> CreateAsymmetricGaussian(double sigmaNegative, double sigmaPositive, const std::vector<double>& x)
{
    std::vector<double> slf(x.size());
    const double s2n = sigmaNegative * sigmaNegative;
    const double s2p = sigmaPositive * sigmaPositive;

    size_t idx = 0;
    while (idx < x.size() && x[idx] < 0.0)
    {
        slf[idx] = std::exp(-(x[idx] * x[idx]) / (2.0 * s2n));
        ++idx;
    }
    while (idx < x.size())
    {
        slf[idx] = std::exp(-(x[idx] * x[idx]) / (2.0 * s2p));
        ++idx;
    }

    return slf;
}


// -------- Fitting symmetrical Gaussians --------
TEST_CASE("FitInstrumentLineShape (Symmetric Gaussian): Simple Gaussian input returns correct fitted function", "[InstrumentLineShape]")
{
    CSpectrum idealGaussian;
    idealGaussian.m_wavelength = CreatePixelToWavelengthMapping(-2.0, +2.0, 43); // 4nm range with 43 sample pointss
    idealGaussian.m_length = (long)idealGaussian.m_wavelength.size();

    SECTION("Narrow line width")
    {
        const double sigma = 0.2; // [nm]
        std::vector<double> spec = CreateGaussian(sigma, idealGaussian.m_wavelength);
        memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));

        Evaluation::GaussianLineShape result;
        auto ret = Evaluation::FitInstrumentLineShape(idealGaussian, result);

        REQUIRE(ret == Evaluation::ILF_RETURN_CODE::SUCCESS);
        REQUIRE(fabs(result.sigma - sigma) < 0.005);
    }

    SECTION("Wide line width")
    {
        const double sigma = 0.8; // [nm]
        std::vector<double> spec = CreateGaussian(sigma, idealGaussian.m_wavelength);
        memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));

        Evaluation::GaussianLineShape result;
        auto ret = Evaluation::FitInstrumentLineShape(idealGaussian, result);

        REQUIRE(ret == Evaluation::ILF_RETURN_CODE::SUCCESS);
        REQUIRE(fabs(result.sigma - sigma) < 0.005);
    }

    SECTION("Really wide line width")
    {
        const double sigma = 1.5; // [nm] this is very wide compared to the width of the wavelength window
        std::vector<double> spec = CreateGaussian(sigma, idealGaussian.m_wavelength);
        memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));

        Evaluation::GaussianLineShape result;
        auto ret = Evaluation::FitInstrumentLineShape(idealGaussian, result);

        REQUIRE(ret == Evaluation::ILF_RETURN_CODE::SUCCESS);
        REQUIRE(fabs(result.sigma - sigma) < 0.005);
    }
}


// -------- Fitting asymmetrical Gaussians --------
TEST_CASE("FitInstrumentLineShape (Asymmetric Gaussian): Simple Asymetric Gaussian input returns correct fitted function", "[InstrumentLineShape]")
{
    CSpectrum idealGaussian;
    idealGaussian.m_wavelength = CreatePixelToWavelengthMapping(-2.0, +2.0, 43); // 4nm range with 43 sample pointss
    idealGaussian.m_length = (long)idealGaussian.m_wavelength.size();

    SECTION("Symmetric input, narrow line width")
    {
        const double sigmaLeft = 0.2; // [nm]
        const double sigmaRight = 0.2; // [nm]
        std::vector<double> spec = CreateAsymmetricGaussian(sigmaLeft, sigmaRight, idealGaussian.m_wavelength);
        memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));

        Evaluation::AsymmetricGaussianLineShape result;
        auto ret = Evaluation::FitInstrumentLineShape(idealGaussian, result);

        REQUIRE(ret == Evaluation::ILF_RETURN_CODE::SUCCESS);
        REQUIRE(fabs(result.sigmaLeft - sigmaLeft) < 0.005);
        REQUIRE(fabs(result.sigmaRight - sigmaRight) < 0.005);
    }

    SECTION("Asymmetric input, narrow line width")
    {
        const double sigmaLeft = 0.2; // [nm]
        const double sigmaRight = 0.5; // [nm]
        std::vector<double> spec = CreateAsymmetricGaussian(sigmaLeft, sigmaRight, idealGaussian.m_wavelength);
        memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));

        Evaluation::AsymmetricGaussianLineShape result;
        auto ret = Evaluation::FitInstrumentLineShape(idealGaussian, result);

        REQUIRE(ret == Evaluation::ILF_RETURN_CODE::SUCCESS);
        REQUIRE(fabs(result.sigmaLeft - sigmaLeft) < 0.005);
        REQUIRE(fabs(result.sigmaRight - sigmaRight) < 0.005);
    }

    SECTION("Asymmetric input, narrow line width #2")
    {
        const double sigmaLeft = 0.5; // [nm]
        const double sigmaRight = 0.2; // [nm]
        std::vector<double> spec = CreateAsymmetricGaussian(sigmaLeft, sigmaRight, idealGaussian.m_wavelength);
        memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));

        Evaluation::AsymmetricGaussianLineShape result;
        auto ret = Evaluation::FitInstrumentLineShape(idealGaussian, result);

        REQUIRE(ret == Evaluation::ILF_RETURN_CODE::SUCCESS);
        REQUIRE(fabs(result.sigmaLeft - sigmaLeft) < 0.005);
        REQUIRE(fabs(result.sigmaRight - sigmaRight) < 0.005);
    }

    SECTION("Very asymmetric input, really wide line width")
    {
        const double sigmaLeft = 0.5; // [nm]
        const double sigmaRight = 1.5; // [nm] this is very wide compared to the width of the wavelength window
        std::vector<double> spec = CreateAsymmetricGaussian(sigmaLeft, sigmaRight, idealGaussian.m_wavelength);
        memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));

        Evaluation::AsymmetricGaussianLineShape result;
        auto ret = Evaluation::FitInstrumentLineShape(idealGaussian, result);

        REQUIRE(ret == Evaluation::ILF_RETURN_CODE::SUCCESS);
        REQUIRE(fabs(result.sigmaLeft - sigmaLeft) < 0.005);
        REQUIRE(fabs(result.sigmaRight - sigmaRight) < 0.005);
    }

    SECTION("Very asymmetric input, really wide line width #2")
    {
        const double sigmaLeft = 1.5; // [nm] this is very wide compared to the width of the wavelength window
        const double sigmaRight = 0.5; // [nm]
        std::vector<double> spec = CreateAsymmetricGaussian(sigmaLeft, sigmaRight, idealGaussian.m_wavelength);
        memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));

        Evaluation::AsymmetricGaussianLineShape result;
        auto ret = Evaluation::FitInstrumentLineShape(idealGaussian, result);

        REQUIRE(ret == Evaluation::ILF_RETURN_CODE::SUCCESS);
        REQUIRE(fabs(result.sigmaLeft - sigmaLeft) < 0.005);
        REQUIRE(fabs(result.sigmaRight - sigmaRight) < 0.005);
    }
}