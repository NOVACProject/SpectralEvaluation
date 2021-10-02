#include "catch.hpp"
#include <vector>
#include <numeric>
#include <SpectralEvaluation/Math/FFT.h>
#include <SpectralEvaluation/VectorUtils.h>

using namespace novac;

TEST_CASE("Fft_Real zero valued input, returns all-zero valued result", "[Math][FFT][Fft_Real]")
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

TEST_CASE("Fft_Real constant valued input, all energy ends up in DC-component", "[Math][FFT][Fft_Real]")
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

TEST_CASE("Fft_Real returns same result as Fft with real valued input", "[Math][FFT][Fft_Real]")
{
    size_t fftLength = 62;
    std::vector<std::complex<double>> outputFromFftReal(fftLength, std::complex<double>{1.0, 1.0});
    std::vector<std::complex<double>> outputFromFft(fftLength, std::complex<double>{1.0, 1.0});

    SECTION("Constant values")
    {
        const double constantInputValue = 3.1;
        const double expectedValue = constantInputValue * fftLength;
        std::vector<double> realValuedInput(fftLength, constantInputValue);
        std::vector<std::complex<double>> complexValuedInput(fftLength, std::complex<double>{constantInputValue, 0.0});

        Fft_Real(realValuedInput, outputFromFftReal);
        Fft(complexValuedInput, outputFromFft);

        // The imaginary components should here be zero in both cases
        REQUIRE(std::abs(Min(Imag(outputFromFftReal))) < std::numeric_limits<float>::epsilon());
        REQUIRE(std::abs(Max(Imag(outputFromFftReal))) < std::numeric_limits<float>::epsilon());
        REQUIRE(std::abs(Min(Imag(outputFromFft))) < std::numeric_limits<float>::epsilon());
        REQUIRE(std::abs(Max(Imag(outputFromFft))) < std::numeric_limits<float>::epsilon());

        // The real component should have the energy collected in the first element
        const auto realOutputFromFftReal = Real(outputFromFftReal);
        const auto realOutputFromFft = Real(outputFromFft);
        REQUIRE(std::abs(expectedValue - realOutputFromFftReal[0]) < std::numeric_limits<float>::epsilon());
        REQUIRE(std::abs(expectedValue - Sum(realOutputFromFftReal)) < std::numeric_limits<float>::epsilon());
        REQUIRE(std::abs(expectedValue - realOutputFromFft[0]) < std::numeric_limits<float>::epsilon());
        REQUIRE(std::abs(expectedValue - Sum(realOutputFromFft)) < std::numeric_limits<float>::epsilon());
    }

    SECTION("Linearly Increasing Values")
    {
        std::vector<double> realValuedInput(fftLength);
        std::iota(begin(realValuedInput), end(realValuedInput), 5.0);
        std::vector<std::complex<double>> complexValuedInput(fftLength);
        std::generate(begin(complexValuedInput), end(complexValuedInput), [n = 0] () mutable {
            return std::complex<double>{5.0 + n++, 0.0};
        });

        Fft_Real(realValuedInput, outputFromFftReal);
        Fft(complexValuedInput, outputFromFft);

        // The imaginary components should be identical
        REQUIRE(std::abs(SumOfSquaredDifferences(Imag(outputFromFftReal), Imag(outputFromFft))) < std::numeric_limits<float>::epsilon());

        // The real components should be identical
        REQUIRE(std::abs(SumOfSquaredDifferences(Real(outputFromFftReal), Real(outputFromFft))) < std::numeric_limits<float>::epsilon());
    }
}

TEST_CASE("Fft_Real constant valued input, length not power of two, all energy ends up in DC-component", "[Math][FFT][Fft_Real]")
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

TEST_CASE("InverseFft restores result of Fft", "[Math][FFT][Fft_Real]")
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

TEST_CASE("InverseFft restores result of FftReal", "[Math][FFT][Fft_Real]")
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
