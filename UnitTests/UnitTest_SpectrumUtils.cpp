#include "catch.hpp"
#include <SpectralEvaluation/Spectra/SpectrumUtils.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>

std::vector<double> CreateGaussian(double center, double sigma, const std::vector<double>& x); // located elsewhere


// -------- Calculating the derivative --------
TEST_CASE("Derivative", "[SpectrumUtils]")
{
    CSpectrum inputSpectrum;
    std::vector<double> result;

    SECTION("Flat Line - returns all zeros")
    {
        inputSpectrum.m_length = 31;
        for (int ii = 0; ii < inputSpectrum.m_length; ++ii)
        {
            inputSpectrum.m_data[ii] = 1.0;
        }

        bool ret = Derivative(inputSpectrum.m_data, inputSpectrum.m_length, result);

        REQUIRE(ret == true);
        REQUIRE(inputSpectrum.m_length == result.size());
        for (size_t ii = 0; ii < result.size(); ++ii)
        {
            REQUIRE(std::abs(result[ii] < 1e-9));
        }
    }
}


// -------- Locating peaks in the spectrum --------
TEST_CASE("FindPeaks", "[SpectrumUtils][FindPeaks]")
{
    CSpectrum inputSpectrum;
    inputSpectrum.m_length = 31;
    inputSpectrum.m_wavelength.resize(inputSpectrum.m_length);
    for (int ii = 0; ii < inputSpectrum.m_length; ++ii)
    {
        inputSpectrum.m_wavelength[ii] = 310.0 + ii * 0.3;
    }

    std::vector<SpectrumDataPoint> result;

    SECTION("Flat Line - returns no peaks")
    {
        for (int ii = 0; ii < inputSpectrum.m_length; ++ii)
        {
            inputSpectrum.m_data[ii] = 1.0;
        }

        FindPeaks(inputSpectrum, 0.0, result);

        REQUIRE(0 == result.size());
    }

    SECTION("Narrow Gaussian 1 - locates peak")
    {
        const double sigma = 0.2; // [nm]
        const double center = inputSpectrum.m_wavelength[inputSpectrum.m_length / 2]; 
        std::vector<double> spec = CreateGaussian(center, sigma, inputSpectrum.m_wavelength);
        memcpy(&inputSpectrum.m_data, spec.data(), spec.size() * sizeof(double));

        FindPeaks(inputSpectrum, 0.0, result);

        REQUIRE(1 == result.size());
        REQUIRE(result[0].pixel == inputSpectrum.m_length / 2);
        REQUIRE(result[0].wavelength == center);
        REQUIRE(result[0].intensity == 1.0);
    }

    SECTION("Narrow Gaussian 2 - locates peak")
    {
        const double sigma = 0.2; // [nm]
        const double center = inputSpectrum.m_wavelength[0] + 0.42 * (inputSpectrum.m_wavelength[inputSpectrum.m_length - 1] - inputSpectrum.m_wavelength[0]);
        std::vector<double> spec = CreateGaussian(center, sigma, inputSpectrum.m_wavelength);
        memcpy(&inputSpectrum.m_data, spec.data(), spec.size() * sizeof(double));

        FindPeaks(inputSpectrum, 0.0, result);

        REQUIRE(1 == result.size());
        REQUIRE(std::abs(result[0].wavelength - center) < 0.01);
    }

    SECTION("D2J2202 Measured Mercury Spectrum - Locates peaks")
    {

    }

}