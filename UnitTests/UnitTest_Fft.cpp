#include "catch.hpp"
#include <vector>
#include <SpectralEvaluation/Math/FFT.h>
#include <SpectralEvaluation/VectorUtils.h>

using namespace novac;

TEST_CASE("Fft_Real zero valued input, returns all-zero valued result", "[FFT][Fft_Real]")
{
    std::vector<double> inputValues(64, 0.0);
    std::vector<std::complex<double>> outputValues(64, std::complex<double>{1.0, 1.0}); // sets all values to one initially

    Fft_Real(inputValues, outputValues);

    SECTION("Real component is all zeros")
    {
        const std::vector<double> realOutputValues = Real(outputValues);

        REQUIRE(std::abs(Min(realOutputValues)) < std::numeric_limits<float>::epsilon());
        REQUIRE(std::abs(Max(realOutputValues)) < std::numeric_limits<float>::epsilon());
    }

    SECTION("Imaginary component is all zeros")
    {
        const std::vector<double> imagOutputValues = Imag(outputValues);

        REQUIRE(std::abs(Min(imagOutputValues)) < std::numeric_limits<float>::epsilon());
        REQUIRE(std::abs(Max(imagOutputValues)) < std::numeric_limits<float>::epsilon());
    }
}

TEST_CASE("Fft_Real constant valued input, all energy ends up in DC-component", "[FFT][Fft_Real]")
{
    size_t fftLength = 64;
    std::vector<double> inputValues(fftLength, 5.0);
    std::vector<std::complex<double>> outputValues(64, std::complex<double>{1.0, 1.0}); // sets all values to one initially

    Fft_Real(inputValues, outputValues);

    SECTION("Imaginary component is all zeros")
    {
        const std::vector<double> imagOutputValues = Imag(outputValues);

        REQUIRE(std::abs(Min(imagOutputValues)) < std::numeric_limits<float>::epsilon());
        REQUIRE(std::abs(Max(imagOutputValues)) < std::numeric_limits<float>::epsilon());
    }

    SECTION("Real component has all energy in first element")
    {
        const double expectedValue = fftLength * inputValues[0]; // this is the total energy
        const std::vector<double> realOutputValues = Real(outputValues);

        REQUIRE(std::abs(expectedValue - realOutputValues[0]) < std::numeric_limits<float>::epsilon());
        REQUIRE(std::abs(expectedValue - Sum(realOutputValues)) < std::numeric_limits<float>::epsilon());
    }
}

TEST_CASE("Fft_Real constant valued input, length not power of two, all energy ends up in DC-component", "[FFT][Fft_Real]")
{
    size_t fftLength = 26; // this is an even number, but not a power of two
    std::vector<double> inputValues(fftLength, 5.0);
    std::vector<std::complex<double>> outputValues(64, std::complex<double>{1.0, 1.0}); // sets all values to one initially

    Fft_Real(inputValues, outputValues);

    SECTION("Imaginary component is all zeros")
    {
        const std::vector<double> imagOutputValues = Imag(outputValues);

        REQUIRE(std::abs(Min(imagOutputValues)) < std::numeric_limits<float>::epsilon());
        REQUIRE(std::abs(Max(imagOutputValues)) < std::numeric_limits<float>::epsilon());
    }

    SECTION("Real component has all energy in first element")
    {
        const double expectedValue = fftLength * inputValues[0]; // this is the total energy
        const std::vector<double> realOutputValues = Real(outputValues);

        REQUIRE(std::abs(expectedValue - realOutputValues[0]) < std::numeric_limits<float>::epsilon());
        REQUIRE(std::abs(expectedValue - Sum(realOutputValues)) < std::numeric_limits<float>::epsilon());
    }
}

TEST_CASE("InverseFft restores result of Fft", "[FFT][Fft_Real]")
{
    size_t fftLength = 63;
    std::vector<std::complex<double>> intermediate(fftLength, std::complex<double>{1.0, 1.0}); // sets all values to one initially
    std::vector<std::complex<double>> output(fftLength, std::complex<double>{1.0, 1.0}); // sets all values to one initially

    SECTION("Constant input values")
    {
        const double inputConstantValue = 5.0;
        std::vector<std::complex<double>> originalRealData(fftLength, std::complex<double>{inputConstantValue, 0.0});

        // First do the forward FFT
        Fft(originalRealData, intermediate, true);

        // Now to the Inverse FFT to (hopefully) restore the original values
        Fft(intermediate, output, false);

        // The result should be equal to the input (but without the normalization factor)
        const std::vector<double> realOutput = Real(output);
        REQUIRE(std::abs(Min(realOutput) - fftLength * inputConstantValue) < std::numeric_limits<float>::epsilon());
        REQUIRE(std::abs(Max(realOutput) - fftLength * inputConstantValue) < std::numeric_limits<float>::epsilon());

        const std::vector<double> imagOutput = Imag(output);
        REQUIRE(std::abs(Min(imagOutput)) < std::numeric_limits<float>::epsilon());
        REQUIRE(std::abs(Max(imagOutput)) < std::numeric_limits<float>::epsilon());
    }
}

TEST_CASE("InverseFft restores result of FftReal", "[FFT][Fft_Real]")
{
    size_t fftLength = 28;
    std::vector<std::complex<double>> intermediate(fftLength, std::complex<double>{1.0, 1.0}); // sets all values to one initially
    std::vector<std::complex<double>> output(fftLength, std::complex<double>{1.0, 1.0}); // sets all values to one initially

    SECTION("Constant input values")
    {
        const double inputConstantValue = 5.0;
        std::vector<double> originalRealData(fftLength, inputConstantValue);

        // First do the forward FFT
        Fft_Real(originalRealData, intermediate);

        // Now to the Inverse FFT to (hopefully) restore the original values
        Fft(intermediate, output, false);

        // The result should be equal to the input (but without the normalization factor)
        const std::vector<double> realOutput = Real(output);
        REQUIRE(std::abs(Min(realOutput) - fftLength * inputConstantValue) < std::numeric_limits<float>::epsilon());
        REQUIRE(std::abs(Max(realOutput) - fftLength * inputConstantValue) < std::numeric_limits<float>::epsilon());

        const std::vector<double> imagOutput = Imag(output);
        REQUIRE(std::abs(Min(imagOutput)) < std::numeric_limits<float>::epsilon());
        REQUIRE(std::abs(Max(imagOutput)) < std::numeric_limits<float>::epsilon());
    }
}
