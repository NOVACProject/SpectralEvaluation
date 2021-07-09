#include "catch.hpp"
#include <string.h>
#include <SpectralEvaluation/Calibration/InstrumentLineShape.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>

using namespace novac;

std::vector<double> CreatePixelToWavelengthMapping(double start, double stop, int size = 2048); // located elsewhere
std::vector<double> CreateGaussian(double sigma, const std::vector<double>& x); // located elsewhere
std::vector<double> CreateGaussian(double center, double sigma, const std::vector<double>& x); // located elsewhere

// creates an asymetric gaussian with width 'sigmaNegative' for x < center and 'sigmaPositive' for x > center
std::vector<double> CreateAsymmetricGaussian(double center, double sigmaNegative, double sigmaPositive, const std::vector<double>& x)
{
    std::vector<double> slf(x.size());
    const double s2n = sigmaNegative * sigmaNegative;
    const double s2p = sigmaPositive * sigmaPositive;

    size_t idx = 0;
    while (idx < x.size() && x[idx] < 0.0)
    {
        const double diff = x[idx] - center;
        slf[idx] = std::exp(-(diff * diff) / (2.0 * s2n));
        ++idx;
    }
    while (idx < x.size())
    {
        const double diff = x[idx] - center;
        slf[idx] = std::exp(-(diff * diff) / (2.0 * s2p));
        ++idx;
    }

    return slf;
}
std::vector<double> CreateAsymmetricGaussian(double sigmaNegative, double sigmaPositive, const std::vector<double>& x)
{
    return CreateAsymmetricGaussian(0.0, sigmaNegative, sigmaPositive, x);
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

        GaussianLineShape result;
        auto ret = FitInstrumentLineShape(idealGaussian, result);

        REQUIRE(ret == ILF_RETURN_CODE::SUCCESS);
        REQUIRE(fabs(result.sigma - sigma) < 0.005);
    }

    SECTION("Wide line width")
    {
        const double sigma = 0.8; // [nm]
        std::vector<double> spec = CreateGaussian(sigma, idealGaussian.m_wavelength);
        memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));

        GaussianLineShape result;
        auto ret = FitInstrumentLineShape(idealGaussian, result);

        REQUIRE(ret == ILF_RETURN_CODE::SUCCESS);
        REQUIRE(fabs(result.sigma - sigma) < 0.005);
    }

    SECTION("Really wide line width")
    {
        const double sigma = 1.5; // [nm] this is very wide compared to the width of the wavelength window
        std::vector<double> spec = CreateGaussian(sigma, idealGaussian.m_wavelength);
        memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));

        GaussianLineShape result;
        auto ret = FitInstrumentLineShape(idealGaussian, result);

        REQUIRE(ret == ILF_RETURN_CODE::SUCCESS);
        REQUIRE(fabs(result.sigma - sigma) < 0.005);
    }
}

TEST_CASE("FitInstrumentLineShape (Symmetric Gaussian): Simple not centered Gaussian input returns correct fitted function", "[InstrumentLineShape]")
{
    const double gaussianCenter = 300.0;
    CSpectrum idealGaussian;
    idealGaussian.m_wavelength = CreatePixelToWavelengthMapping(gaussianCenter - 2.0, gaussianCenter + 2.0, 43); // 4nm range with 43 sample pointss
    idealGaussian.m_length = (long)idealGaussian.m_wavelength.size();

    SECTION("Narrow line width")
    {
        const double sigma = 0.2; // [nm]
        std::vector<double> spec = CreateGaussian(gaussianCenter, sigma, idealGaussian.m_wavelength);
        memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));

        GaussianLineShape result;
        auto ret = FitInstrumentLineShape(idealGaussian, result);

        REQUIRE(ret == ILF_RETURN_CODE::SUCCESS);
        REQUIRE(fabs(result.sigma - sigma) < 0.005);
    }

    SECTION("Wide line width")
    {
        const double sigma = 0.8; // [nm]
        std::vector<double> spec = CreateGaussian(gaussianCenter, sigma, idealGaussian.m_wavelength);
        memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));

        GaussianLineShape result;
        auto ret = FitInstrumentLineShape(idealGaussian, result);

        REQUIRE(ret == ILF_RETURN_CODE::SUCCESS);
        REQUIRE(fabs(result.sigma - sigma) < 0.005);
    }

    SECTION("Really wide line width")
    {
        const double sigma = 1.5; // [nm] this is very wide compared to the width of the wavelength window
        std::vector<double> spec = CreateGaussian(gaussianCenter, sigma, idealGaussian.m_wavelength);
        memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));

        GaussianLineShape result;
        auto ret = FitInstrumentLineShape(idealGaussian, result);

        REQUIRE(ret == ILF_RETURN_CODE::SUCCESS);
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

        AsymmetricGaussianLineShape result;
        auto ret = FitInstrumentLineShape(idealGaussian, result);

        REQUIRE(ret == ILF_RETURN_CODE::SUCCESS);
        REQUIRE(fabs(result.sigmaLeft - sigmaLeft) < 0.005);
        REQUIRE(fabs(result.sigmaRight - sigmaRight) < 0.005);
    }

    SECTION("Asymmetric input, narrow line width")
    {
        const double sigmaLeft = 0.2; // [nm]
        const double sigmaRight = 0.5; // [nm]
        std::vector<double> spec = CreateAsymmetricGaussian(sigmaLeft, sigmaRight, idealGaussian.m_wavelength);
        memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));

        AsymmetricGaussianLineShape result;
        auto ret = FitInstrumentLineShape(idealGaussian, result);

        REQUIRE(ret == ILF_RETURN_CODE::SUCCESS);
        REQUIRE(fabs(result.sigmaLeft - sigmaLeft) < 0.005);
        REQUIRE(fabs(result.sigmaRight - sigmaRight) < 0.005);
    }

    SECTION("Asymmetric input, narrow line width #2")
    {
        const double sigmaLeft = 0.5; // [nm]
        const double sigmaRight = 0.2; // [nm]
        std::vector<double> spec = CreateAsymmetricGaussian(sigmaLeft, sigmaRight, idealGaussian.m_wavelength);
        memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));

        AsymmetricGaussianLineShape result;
        auto ret = FitInstrumentLineShape(idealGaussian, result);

        REQUIRE(ret == ILF_RETURN_CODE::SUCCESS);
        REQUIRE(fabs(result.sigmaLeft - sigmaLeft) < 0.005);
        REQUIRE(fabs(result.sigmaRight - sigmaRight) < 0.005);
    }

    SECTION("Very asymmetric input, really wide line width")
    {
        const double sigmaLeft = 0.5; // [nm]
        const double sigmaRight = 1.5; // [nm] this is very wide compared to the width of the wavelength window
        std::vector<double> spec = CreateAsymmetricGaussian(sigmaLeft, sigmaRight, idealGaussian.m_wavelength);
        memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));

        AsymmetricGaussianLineShape result;
        auto ret = FitInstrumentLineShape(idealGaussian, result);

        REQUIRE(ret == ILF_RETURN_CODE::SUCCESS);
        REQUIRE(fabs(result.sigmaLeft - sigmaLeft) < 0.005);
        REQUIRE(fabs(result.sigmaRight - sigmaRight) < 0.005);
    }

    SECTION("Very asymmetric input, really wide line width #2")
    {
        const double sigmaLeft = 1.5; // [nm] this is very wide compared to the width of the wavelength window
        const double sigmaRight = 0.5; // [nm]
        std::vector<double> spec = CreateAsymmetricGaussian(sigmaLeft, sigmaRight, idealGaussian.m_wavelength);
        memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));

        AsymmetricGaussianLineShape result;
        auto ret = FitInstrumentLineShape(idealGaussian, result);

        REQUIRE(ret == ILF_RETURN_CODE::SUCCESS);
        REQUIRE(fabs(result.sigmaLeft - sigmaLeft) < 0.005);
        REQUIRE(fabs(result.sigmaRight - sigmaRight) < 0.005);
    }
}

// -------- Fitting asymmetrical Gaussians --------
// TEST_CASE("FitInstrumentLineShape (Asymmetric Gaussian): Simple not centered Asymetric Gaussian input returns correct fitted function", "[InstrumentLineShape]")
// {
//     const double gaussianCenter = 300.0;
//     CSpectrum idealGaussian;
//     idealGaussian.m_wavelength = CreatePixelToWavelengthMapping(gaussianCenter - 2.0, gaussianCenter + 2.0, 43); // 4nm range with 43 sample pointss
//     idealGaussian.m_length = (long)idealGaussian.m_wavelength.size();
// 
//     SECTION("Symmetric input, narrow line width")
//     {
//         const double sigmaLeft = 0.2; // [nm]
//         const double sigmaRight = 0.2; // [nm]
//         std::vector<double> spec = CreateAsymmetricGaussian(gaussianCenter, sigmaLeft, sigmaRight, idealGaussian.m_wavelength);
//         memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));
// 
//         AsymmetricGaussianLineShape result;
//         auto ret = FitInstrumentLineShape(idealGaussian, result);
// 
//         REQUIRE(ret == ILF_RETURN_CODE::SUCCESS);
//         REQUIRE(fabs(result.sigmaLeft - sigmaLeft) < 0.005);
//         REQUIRE(fabs(result.sigmaRight - sigmaRight) < 0.005);
//     }
// 
//     SECTION("Asymmetric input, narrow line width")
//     {
//         const double sigmaLeft = 0.2; // [nm]
//         const double sigmaRight = 0.5; // [nm]
//         std::vector<double> spec = CreateAsymmetricGaussian(gaussianCenter, sigmaLeft, sigmaRight, idealGaussian.m_wavelength);
//         memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));
// 
//         AsymmetricGaussianLineShape result;
//         auto ret = FitInstrumentLineShape(idealGaussian, result);
// 
//         REQUIRE(ret == ILF_RETURN_CODE::SUCCESS);
//         REQUIRE(fabs(result.sigmaLeft - sigmaLeft) < 0.005);
//         REQUIRE(fabs(result.sigmaRight - sigmaRight) < 0.005);
//     }
// 
//     SECTION("Asymmetric input, narrow line width #2")
//     {
//         const double sigmaLeft = 0.5; // [nm]
//         const double sigmaRight = 0.2; // [nm]
//         std::vector<double> spec = CreateAsymmetricGaussian(gaussianCenter, sigmaLeft, sigmaRight, idealGaussian.m_wavelength);
//         memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));
// 
//         AsymmetricGaussianLineShape result;
//         auto ret = FitInstrumentLineShape(idealGaussian, result);
// 
//         REQUIRE(ret == ILF_RETURN_CODE::SUCCESS);
//         REQUIRE(fabs(result.sigmaLeft - sigmaLeft) < 0.005);
//         REQUIRE(fabs(result.sigmaRight - sigmaRight) < 0.005);
//     }
// 
//     SECTION("Very asymmetric input, really wide line width")
//     {
//         const double sigmaLeft = 0.5; // [nm]
//         const double sigmaRight = 1.5; // [nm] this is very wide compared to the width of the wavelength window
//         std::vector<double> spec = CreateAsymmetricGaussian(gaussianCenter, sigmaLeft, sigmaRight, idealGaussian.m_wavelength);
//         memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));
// 
//         AsymmetricGaussianLineShape result;
//         auto ret = FitInstrumentLineShape(idealGaussian, result);
// 
//         REQUIRE(ret == ILF_RETURN_CODE::SUCCESS);
//         REQUIRE(fabs(result.sigmaLeft - sigmaLeft) < 0.005);
//         REQUIRE(fabs(result.sigmaRight - sigmaRight) < 0.005);
//     }
// 
//     SECTION("Very asymmetric input, really wide line width #2")
//     {
//         const double sigmaLeft = 1.5; // [nm] this is very wide compared to the width of the wavelength window
//         const double sigmaRight = 0.5; // [nm]
//         std::vector<double> spec = CreateAsymmetricGaussian(gaussianCenter, sigmaLeft, sigmaRight, idealGaussian.m_wavelength);
//         memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));
// 
//         AsymmetricGaussianLineShape result;
//         auto ret = FitInstrumentLineShape(idealGaussian, result);
// 
//         REQUIRE(ret == ILF_RETURN_CODE::SUCCESS);
//         REQUIRE(fabs(result.sigmaLeft - sigmaLeft) < 0.005);
//         REQUIRE(fabs(result.sigmaRight - sigmaRight) < 0.005);
//     }
// }


// -------- Fitting symmetrical SuperGaussian --------
TEST_CASE("FitInstrumentLineShape (Symmetric Super Gaussian): Simple Gaussian input returns correct fitted function", "[InstrumentLineShape][SuperGaussian]")
{
    CSpectrum idealGaussian;
    idealGaussian.m_wavelength = CreatePixelToWavelengthMapping(-2.0, +2.0, 43); // 4nm range with 43 sample pointss
    idealGaussian.m_length = (long)idealGaussian.m_wavelength.size();

    SECTION("Narrow line width")
    {
        const double sigma = 0.2; // [nm]
        std::vector<double> spec = CreateGaussian(sigma, idealGaussian.m_wavelength);
        memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));

        SuperGaussianLineShape result;
        auto ret = FitInstrumentLineShape(idealGaussian, result);

        REQUIRE(ret == ILF_RETURN_CODE::SUCCESS);
        REQUIRE(fabs(result.sigma - sigma) < 0.005);
        REQUIRE(fabs(result.P - 2) < 0.0005);
    }
}

TEST_CASE("FitInstrumentLineShape (Symmetric Super Gaussian): Simple not centered Gaussian input returns correct fitted function", "[InstrumentLineShape][SuperGaussian]")
{
    const double gaussianCenter = 300.0;
    CSpectrum idealGaussian;
    idealGaussian.m_wavelength = CreatePixelToWavelengthMapping(gaussianCenter - 2.0, gaussianCenter + 2.0, 43); // 4nm range with 43 sample pointss
    idealGaussian.m_length = (long)idealGaussian.m_wavelength.size();

    SECTION("Narrow line width")
    {
        const double sigma = 0.2; // [nm]
        std::vector<double> spec = CreateGaussian(gaussianCenter, sigma, idealGaussian.m_wavelength);
        memcpy(&idealGaussian.m_data, spec.data(), spec.size() * sizeof(double));

        SuperGaussianLineShape result;
        auto ret = FitInstrumentLineShape(idealGaussian, result);

        REQUIRE(ret == ILF_RETURN_CODE::SUCCESS);
        REQUIRE(fabs(result.sigma - sigma) < 0.005);
        REQUIRE(fabs(result.P - 2) < 0.0005);
    }
}


// -------- Fitting to real data--------
TEST_CASE("FitInstrumentLineShape: Real SLF input returns reasonable fitted function", "[InstrumentLineShape]")
{
    std::vector<double> xData = { -1.823155769, -1.740163507, -1.657182786, -1.574213609,
        -1.491255984,     -1.408309915,     -1.325375407,     -1.242452466,
        -1.159541097,     -1.076641306,    -0.993753097,    -0.910876477,
        -0.82801145,    -0.745158022,    -0.662316199,    -0.579485984,
        -0.496667385,    -0.413860406,    -0.331065052,    -0.248281329,
        -0.165509243,    -0.082748798,    0.0,    0.082737146,
        0.165462634,    0.248176459,    0.330878616,    0.4135691,
        0.496247904,    0.578915024,    0.661570455,    0.74421419,
        0.826846226,    0.909466555,    0.992075174,    1.074672076,
        1.157257257,    1.23983071 ,    1.322392432,    1.404942415,
        1.487480656,    1.570007148,    1.652521887,    1.735024866,
        1.817516081 };

    std::vector<double> yData = { 7.897011196, 8.702061524, 9.416269597, 12.20731957,
        16.3484734 , 19.98842858, 22.42863949, 25.44209636,
        28.80952478, 30.98034142, 33.33284258, 37.57423613,
        45.10161333, 60.24407748, 95.47834242, 211.5371543,
        671.311734 , 1918.519509, 4204.884369, 7004.03183 ,
        8959.771604, 9851.169062, 10000      , 9649.502388,
        8045.297081, 5455.735232, 3228.467958, 1912.639823,
        1102.154622, 691.2187178, 484.1829538, 334.8163592,
        206.0427641, 108.4719174, 57.91350376, 37.55230869,
        25.63004586, 20.04794592, 13.96464733, 11.2331147 ,
        8.567364387, 5.785711892, 2.822374886, 2.220936508,
        0 };

    CSpectrum measuredSpectrum;
    measuredSpectrum.m_wavelength = xData;
    measuredSpectrum.m_length = (long)xData.size();
    memcpy(&measuredSpectrum.m_data, yData.data(), yData.size() * sizeof(double));

    SECTION("Asmmetric Gaussian")
    {
        AsymmetricGaussianLineShape result;
        auto ret = FitInstrumentLineShape(measuredSpectrum, result);

        REQUIRE(ret == ILF_RETURN_CODE::SUCCESS);
        REQUIRE(fabs(result.sigmaLeft - 0.20991) < 0.005);
        REQUIRE(fabs(result.sigmaRight - 0.25473) < 0.005);
    }

    SECTION("Symmetric Gaussian")
    {
        GaussianLineShape result;
        auto ret = FitInstrumentLineShape(measuredSpectrum, result);

        REQUIRE(ret == ILF_RETURN_CODE::SUCCESS);
        REQUIRE(fabs(result.sigma - 0.23295) < 0.005);
    }

    SECTION("Symmetric Super Gaussian")
    {
        SuperGaussianLineShape result;
        auto ret = FitInstrumentLineShape(measuredSpectrum, result);

        REQUIRE(ret == ILF_RETURN_CODE::SUCCESS);
        REQUIRE(fabs(result.sigma - 0.249351) < 0.005);
        REQUIRE(fabs(result.P - 2.37497) < 0.0005);
    }
}