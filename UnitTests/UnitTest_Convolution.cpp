#include "catch.hpp"
#include <Spectra/ReferenceSpectrumConvolution.h>
#include <VectorUtils.h>

std::vector<double> CreateGaussian(double sigmaInPixels, int length = 19)
{
    const int halfLength    = (length / 2);
    const int midIdx        = halfLength + 1;
    const double s2         = sigmaInPixels * sigmaInPixels;

    std::vector<double> slf(length);
    for (int dI = 0; dI < halfLength; ++dI)
    {
        const double value = std::exp(-(dI * dI) / (2.0 * s2));
        slf[midIdx + dI] = value;
        slf[midIdx - dI] = value;
    }

    return slf;
}

std::vector<double> CreateGaussian(double sigmaInPixels, const std::vector<double>& x)
{
    const double s2 = sigmaInPixels * sigmaInPixels;

    std::vector<double> slf(x.size());
    for(size_t i = 0; i < x.size(); ++i)
    {
        slf[i] = std::exp(-(x[i] * x[i]) / (2.0 * s2));
    }

    return slf;
}


std::vector<double> CreatePixelToWavelengthMapping(double start, double stop, int size = 2048)
{
    const double inv_size   = 1.0 / (double)(size - 1);
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
    SimpleSpectrum slf;
    slf.wavelength  = CreatePixelToWavelengthMapping(-2.0, +2.0, 43);
    slf.value       = CreateGaussian(sigmaInPixels, slf.wavelength);

    const double expectedFwhm = 2.3548 * sigmaInPixels;

    const double result = CalculateFhwm(slf);
    REQUIRE( fabs(result -expectedFwhm) < 0.005);
}

TEST_CASE("Convolve returns expected output - simple input.", "[Convolve]")
{
    SimpleSpectrum highResReference;
    highResReference.wavelength = CreatePixelToWavelengthMapping(200.0, 300.0, 1000); // 0.1 nm resolution
    highResReference.value      = std::vector<double>(1000, 0.0 ); // all-zeros
    highResReference.value[500] = 1.0; // one spike
    REQUIRE(highResReference.value.size() == highResReference.wavelength.size());

    SimpleSpectrum slf;
    const double slfSigma   = 0.7;
    slf.wavelength          = CreatePixelToWavelengthMapping(-2.0, +2.0, 41); // 0.1 nm resolution
    slf.value               = CreateGaussian(slfSigma, slf.wavelength);

    std::vector<double> result;

    SECTION("Result has correct length")
    {
        bool retVal = Convolve(slf, highResReference, result);
    
        REQUIRE(retVal == true);
        REQUIRE(result.size() == highResReference.wavelength.size());
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

        SimpleSpectrum resultingCrossSection;
        resultingCrossSection.value         = result;
        resultingCrossSection.wavelength    = highResReference.wavelength;

        const double fwhm = CalculateFhwm(resultingCrossSection);

        const double expectedFwhm = slfSigma * 2.3548;

        REQUIRE( (fwhm / expectedFwhm) < 1.04);
        REQUIRE( (fwhm / expectedFwhm) > 0.96);
    }

    SECTION("Same output even if slf is not normalized")
    {
        // Increase the amplitude of the SLF.
        for (double& v : slf.value)
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

    SimpleSpectrum highResReference;
    highResReference.wavelength = CreatePixelToWavelengthMapping(wavelMapping.front(), wavelMapping.back(), 8192);
    highResReference.value      = std::vector<double>(8192, 0.0); // all-zeros

    // Add one spike to the high-res reference so that we can study what happens to it.
    const int idxOfSpikeInHighResReference = 500;
    highResReference.value[idxOfSpikeInHighResReference] = 1.0; // one spike

    SimpleSpectrum slf;
    const double slfSigma   = 0.7;
    slf.wavelength          = CreatePixelToWavelengthMapping(-2.0, +2.0, 41); // 0.1 nm resolution
    slf.value               = CreateGaussian(slfSigma, slf.wavelength);

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
        const double wavelengthOfSpike = highResReference.wavelength[idxOfSpikeInHighResReference];
        const double idxInWavelMappingF = FindValue(wavelMapping, wavelengthOfSpike, 0U, wavelMapping.size());
        const size_t idxInWavelMapping = (size_t)std::round(idxInWavelMappingF);

        bool retVal = ConvolveReference(wavelMapping, slf, highResReference, result);

        REQUIRE(retVal == true);
        REQUIRE(fabs(result[idxInWavelMapping] - 1.0) < 0.01); // the spike must be at the correct point
    }

    /* SECTION("Output has correct width (in nm)")
    {
        bool retVal = ConvolveReference(wavelMapping, slf, highResReference, result);
        REQUIRE(retVal == true);

        SimpleSpectrum resultingCrossSection;
        resultingCrossSection.value = result;
        resultingCrossSection.wavelength = highResReference.wavelength;

        const double fwhm = CalculateFhwm(resultingCrossSection);

        const double expectedFwhm = slfSigma * 2.3548;

        REQUIRE((fwhm / expectedFwhm) < 1.04);
        REQUIRE((fwhm / expectedFwhm) > 0.96);
    } */

    SECTION("Same output even if slf is not normalized")
    {
        // Find the index of the same wavelength in wavelMapping
        const double wavelengthOfSpike = highResReference.wavelength[idxOfSpikeInHighResReference];
        const double idxInWavelMappingF = FindValue(wavelMapping, wavelengthOfSpike, 0U, wavelMapping.size());
        const size_t idxInWavelMapping = (size_t)std::round(idxInWavelMappingF);

        // Increase the amplitude of the SLF.
        for (double& v : slf.value)
        {
            v *= 123.0;
        }

        bool retVal = ConvolveReference(wavelMapping, slf, highResReference, result);

        REQUIRE(retVal == true);
        REQUIRE(fabs(result[idxInWavelMapping] - 1.0) < 0.01); // the spike must be at the correct point and still have amplitude = 1
    }

    // TODO: what happens if the reference doesn't cover the same wavelength range as the spectrometer ??
}


TEST_CASE("ConvolveReference returns expected output - reference lambda_max < spectrometer lambda_max ", "[ConvolveReference]")
{
    std::vector<double> wavelMapping = CreatePixelToWavelengthMapping(278.0, 423.0, 2048);

    SimpleSpectrum highResReference;
    highResReference.wavelength = CreatePixelToWavelengthMapping(wavelMapping.front(), 350.0, 8192);
    highResReference.value = std::vector<double>(8192, 1.0); // all-ones

    SimpleSpectrum slf;
    const double slfSigma = 0.7;
    slf.wavelength = CreatePixelToWavelengthMapping(-2.0, +2.0, 41); // 0.1 nm resolution
    slf.value = CreateGaussian(slfSigma, slf.wavelength);

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
        const double idxInWavelMappingF = FindValue(wavelMapping, highResReference.wavelength.back(), 0U, wavelMapping.size());
        const size_t idxInWavelMapping = (size_t)std::round(idxInWavelMappingF);

        bool retVal = ConvolveReference(wavelMapping, slf, highResReference, result);

        REQUIRE(retVal == true);
        REQUIRE(fabs(result[idxInWavelMapping + 50]) < 0.01); // the result must be zero after the refer
    }
}

