#include "catch.hpp"
#include <SpectralEvaluation/Calibration/ReferenceSpectrumConvolution.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/VectorUtils.h>

std::vector<double> CreateGaussian(double sigma, int length = 19)
{
    const int halfLength = (length / 2);
    const int midIdx = halfLength + 1;
    const double s2 = sigma * sigma;

    std::vector<double> slf(length);
    for (int dI = 0; dI < halfLength; ++dI)
    {
        const double value = std::exp(-(dI * dI) / (2.0 * s2));
        slf[midIdx + dI] = value;
        slf[midIdx - dI] = value;
    }

    return slf;
}

std::vector<double> CreateGaussian(double sigma, const std::vector<double>& x)
{
    const double s2 = sigma * sigma;

    std::vector<double> slf(x.size());
    for (size_t i = 0; i < x.size(); ++i)
    {
        slf[i] = std::exp(-(x[i] * x[i]) / (2.0 * s2));
    }

    return slf;
}

std::vector<double> CreateGaussian(double center, double sigma, const std::vector<double>& x)
{
    const double s2 = sigma * sigma;

    std::vector<double> slf(x.size());
    for (size_t i = 0; i < x.size(); ++i)
    {
        const double diff = x[i] - center;
        slf[i] = std::exp(-(diff * diff) / (2.0 * s2));
    }

    return slf;
}



std::vector<double> CreatePixelToWavelengthMapping(double start, double stop, int size = 2048)
{
    const double inv_size = 1.0 / (double)(size - 1);
    std::vector<double> wavelength(size, 0.0);

    for (int ii = 0; ii < size; ++ii)
    {
        wavelength[ii] = start + (stop - start) * ii * inv_size;
    }

    return wavelength;
}

TEST_CASE("CalculateFwhm returns correct value")
{
    const double sigmaInPixels = 0.5;
    Evaluation::CCrossSectionData slf;
    slf.m_waveLength = CreatePixelToWavelengthMapping(-2.0, +2.0, 43);
    slf.m_crossSection = CreateGaussian(sigmaInPixels, slf.m_waveLength);

    const double expectedFwhm = 2.3548 * sigmaInPixels;

    const double result = CalculateFhwm(slf);
    REQUIRE(fabs(result - expectedFwhm) < 0.005);
}

TEST_CASE("Convolve returns expected output - simple input.", "[Convolve]")
{
    Evaluation::CCrossSectionData highResReference;
    highResReference.m_waveLength = CreatePixelToWavelengthMapping(200.0, 300.0, 1000); // 0.1 nm resolution
    highResReference.m_crossSection = std::vector<double>(1000, 0.0); // all-zeros
    highResReference.m_crossSection[500] = 1.0; // one spike
    REQUIRE(highResReference.m_crossSection.size() == highResReference.m_waveLength.size());

    Evaluation::CCrossSectionData slf;
    const double slfSigma = 0.7;
    slf.m_waveLength = CreatePixelToWavelengthMapping(-2.0, +2.0, 41); // 0.1 nm resolution
    slf.m_crossSection = CreateGaussian(slfSigma, slf.m_waveLength);

    std::vector<double> result;

    SECTION("Result has correct length")
    {
        bool retVal = Convolve(slf, highResReference, result);

        REQUIRE(retVal == true);
        REQUIRE(result.size() == highResReference.m_waveLength.size());
    }

    SECTION("Output is centered correctly")
    {
        bool retVal = Convolve(slf, highResReference, result);

        REQUIRE(retVal == true);
        REQUIRE(fabs(result[500] - 1.0) < 0.03); // the spike must be at the correct point
    }

    SECTION("Output has correct width (in nm)")
    {
        bool retVal = Convolve(slf, highResReference, result);
        REQUIRE(retVal == true);

        Evaluation::CCrossSectionData resultingCrossSection;
        resultingCrossSection.m_crossSection = result;
        resultingCrossSection.m_waveLength = highResReference.m_waveLength;

        const double fwhm = CalculateFhwm(resultingCrossSection);

        const double expectedFwhm = slfSigma * 2.3548;

        REQUIRE((fwhm / expectedFwhm) < 1.04);
        REQUIRE((fwhm / expectedFwhm) > 0.94);
    }

    SECTION("Same output even if slf is not normalized")
    {
        // Increase the amplitude of the SLF.
        for (double& v : slf.m_crossSection)
        {
            v *= 123.0;
        }

        bool retVal = Convolve(slf, highResReference, result);

        REQUIRE(retVal == true);
        REQUIRE(fabs(result[500] - 1.0) < 0.03); // the spike must be at the correct point and still have amplitude = 1
    }
}



TEST_CASE("ConvolveReference returns expected output - simple input", "[ConvolveReference]")
{
    std::vector<double> wavelMapping = CreatePixelToWavelengthMapping(278.0, 423.0, 2048);

    Evaluation::CCrossSectionData highResReference;
    highResReference.m_waveLength = CreatePixelToWavelengthMapping(wavelMapping.front(), wavelMapping.back(), 8192);
    highResReference.m_crossSection = std::vector<double>(8192, 0.0); // all-zeros

    // Add one spike to the high-res reference so that we can study what happens to it.
    const int idxOfSpikeInHighResReference = 500;
    highResReference.m_crossSection[idxOfSpikeInHighResReference] = 1.0; // one spike

    Evaluation::CCrossSectionData slf;
    const double slfSigma = 0.7;
    slf.m_waveLength = CreatePixelToWavelengthMapping(-2.0, +2.0, 41); // 0.1 nm resolution
    slf.m_crossSection = CreateGaussian(slfSigma, slf.m_waveLength);

    std::vector<double> result;

    SECTION("Result has correct length")
    {
        bool retVal = ConvolveReference(wavelMapping, slf, highResReference, result);

        REQUIRE(retVal == true);
        REQUIRE(result.size() == wavelMapping.size());
    }

    SECTION("Output is centered correctly")
    {
        // Find the index of the same wavelength in wavelMapping
        const double wavelengthOfSpike = highResReference.m_waveLength[idxOfSpikeInHighResReference];
        const double idxInWavelMappingF = FindValue(wavelMapping, wavelengthOfSpike, 0U, wavelMapping.size());
        const size_t idxInWavelMapping = (size_t)std::round(idxInWavelMappingF);

        bool retVal = ConvolveReference(wavelMapping, slf, highResReference, result);

        const double expectedAmplitude = 0.01;
        REQUIRE(retVal == true);
        REQUIRE(fabs(result[idxInWavelMapping] - expectedAmplitude) < 0.01); // the spike must be at the correct point
    }

    /* SECTION("Output has correct width (in nm)")
    {
        bool retVal = ConvolveReference(wavelMapping, slf, highResReference, result);
        REQUIRE(retVal == true);

        Evaluation::CCrossSectionData resultingCrossSection;
        resultingCrossSection.m_crossSection = result;
        resultingCrossSection.m_waveLength = highResReference.m_waveLength;

        const double fwhm = CalculateFhwm(resultingCrossSection);

        const double expectedFwhm = slfSigma * 2.3548;

        REQUIRE((fwhm / expectedFwhm) < 1.04);
        REQUIRE((fwhm / expectedFwhm) > 0.96);
    } */

    SECTION("Same output even if slf is not normalized")
    {
        // Find the index of the same wavelength in wavelMapping
        const double wavelengthOfSpike = highResReference.m_waveLength[idxOfSpikeInHighResReference];
        const double idxInWavelMappingF = FindValue(wavelMapping, wavelengthOfSpike, 0U, wavelMapping.size());
        const size_t idxInWavelMapping = (size_t)std::round(idxInWavelMappingF);

        // Increase the amplitude of the SLF.
        for (double& v : slf.m_crossSection)
        {
            v *= 123.0;
        }

        bool retVal = ConvolveReference(wavelMapping, slf, highResReference, result);

        const double expectedAmplitude = 0.01;
        REQUIRE(retVal == true);
        REQUIRE(fabs(result[idxInWavelMapping] - expectedAmplitude) < 0.01); // the spike must be at the correct point and still have amplitude = 1
    }

    // TODO: what happens if the reference doesn't cover the same wavelength range as the spectrometer ??
}


TEST_CASE("ConvolveReference (FFT) returns expected output - simple input", "[ConvolveReference][FFT]")
{
    std::vector<double> wavelMapping = CreatePixelToWavelengthMapping(278.0, 423.0, 2048);

    Evaluation::CCrossSectionData highResReference;
    highResReference.m_waveLength = CreatePixelToWavelengthMapping(wavelMapping.front(), wavelMapping.back(), 8192);
    highResReference.m_crossSection = std::vector<double>(8192, 0.0); // all-zeros

    // Add one spike to the high-res reference so that we can study what happens to it.
    const int idxOfSpikeInHighResReference = 500;
    highResReference.m_crossSection[idxOfSpikeInHighResReference] = 1.0; // one spike

    Evaluation::CCrossSectionData slf;
    const double slfSigma = 0.7;
    slf.m_waveLength = CreatePixelToWavelengthMapping(-2.0, +2.0, 41); // 0.1 nm resolution
    slf.m_crossSection = CreateGaussian(slfSigma, slf.m_waveLength);

    SECTION("Result has correct length")
    {
        std::vector<double> result;
        bool retVal = ConvolveReference(wavelMapping, slf, highResReference, result, WavelengthConversion::None, ConvolutionMethod::Fft);

        REQUIRE(retVal == true);
        REQUIRE(result.size() == wavelMapping.size());
    }

    SECTION("Result is identical to output from direct convolution")
    {
        std::vector<double> resultFromFft;
        ConvolveReference(wavelMapping, slf, highResReference, resultFromFft, WavelengthConversion::None, ConvolutionMethod::Fft);

        std::vector<double> resultFromDirectConv;
        ConvolveReference(wavelMapping, slf, highResReference, resultFromDirectConv, WavelengthConversion::None, ConvolutionMethod::Direct);

        // The two vectors should be identical
        REQUIRE(resultFromFft.size() == resultFromDirectConv.size());
        REQUIRE(std::abs(SumOfSquaredDifferences(resultFromFft, resultFromDirectConv)) < std::numeric_limits<float>::epsilon());
    }
}


TEST_CASE("ConvolveReference returns expected output - reference lambda_max < spectrometer lambda_max ", "[ConvolveReference]")
{
    std::vector<double> wavelMapping = CreatePixelToWavelengthMapping(278.0, 423.0, 2048);

    Evaluation::CCrossSectionData highResReference;
    highResReference.m_waveLength = CreatePixelToWavelengthMapping(wavelMapping.front(), 350.0, 8192);
    highResReference.m_crossSection = std::vector<double>(8192, 1.0); // all-ones

    Evaluation::CCrossSectionData slf;
    const double slfSigma = 0.7;
    slf.m_waveLength = CreatePixelToWavelengthMapping(-2.0, +2.0, 41); // 0.1 nm resolution
    slf.m_crossSection = CreateGaussian(slfSigma, slf.m_waveLength);

    std::vector<double> result;

    SECTION("Result has correct length")
    {
        bool retVal = ConvolveReference(wavelMapping, slf, highResReference, result);

        REQUIRE(retVal == true);
        REQUIRE(result.size() == wavelMapping.size());
    }

    SECTION("Output is zero after reference ends.")
    {
        // Find the index of the same wavelength in wavelMapping
        const double idxInWavelMappingF = FindValue(wavelMapping, highResReference.m_waveLength.back(), 0U, wavelMapping.size());
        const size_t idxInWavelMapping = (size_t)std::round(idxInWavelMappingF);

        bool retVal = ConvolveReference(wavelMapping, slf, highResReference, result);

        REQUIRE(retVal == true);
        REQUIRE(fabs(result[idxInWavelMapping + 50]) < 0.01); // the result must be zero after the refer
    }
}

