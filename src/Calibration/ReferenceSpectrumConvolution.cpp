#include <SpectralEvaluation/Calibration/ReferenceSpectrumConvolution.h>
#include <SpectralEvaluation/Calibration/InstrumentLineShapeEstimation.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Fit/CubicSplineFunction.h>
#include <SpectralEvaluation/File/File.h>
#include <SpectralEvaluation/VectorUtils.h>
#include <SpectralEvaluation/Spectra/Grid.h>
#include <SpectralEvaluation/Air.h>
#include <SpectralEvaluation/Math/FFT.h>
#include <iostream>
#include <assert.h>

namespace novac
{

double CalculateFhwm(const CCrossSectionData& slf)
{
    // 1. Find the maximum value
    size_t idxOfMax = 0U;
    const double maxValue = Max(slf.m_crossSection, idxOfMax);

    // 2. Find the (fractional) index to the left and to the right where the amplitude has fallen to half.
    const double leftIdx = FindValue(slf.m_crossSection, maxValue * 0.5, 0U, idxOfMax - 1);
    const double rightIdx = FindValue(slf.m_crossSection, maxValue * 0.5, idxOfMax + 1, slf.m_crossSection.size());

    // 3. Get the x-axis values at the respective indices
    const double leftX = GetAt(slf.m_waveLength, leftIdx);
    const double rightX = GetAt(slf.m_waveLength, rightIdx);

    // 4. We now have the FWHM
    return (rightX - leftX);
}

/* Performs a convolution between the input and core and stores the result in 'result' */
void ConvolutionCore(const std::vector<double>& input, const std::vector<double>& core, std::vector<double>& result)
{
    const size_t inputSize = input.size();
    const size_t coreSize = core.size();

    result.resize(inputSize + coreSize - 1, 0.0);

    // The actual convolution. Here a dead-simple raw convolution calculation. This can be made faster if required.
    // The convolution is separated into three parts, the beginning, the middle and the end, for efficiently reasons since
    // some math in the inner loop can be avoided in the major loop this way.
    //  Get the pointers to the data, this avoids range-checking for each and every value.
    const double* inputPtr = input.data();
    const double* corePtr = core.data();
    double* resultPtr = result.data();

    // The first part, until the entire core fits into the input.
    for (std::int64_t n = coreSize / 2; n < (std::int64_t)coreSize - 1; ++n)
    {
        resultPtr[n] = 0;

        const std::int64_t kmin = n - (coreSize - 1);
        const std::int64_t kmax = n;

        for (std::int64_t k = std::max(kmin, (std::int64_t)0); k <= kmax; k++)
        {
            resultPtr[n] += inputPtr[k] * corePtr[std::max(n - k, (std::int64_t)0)];
        }
    }

    // The middle, this is normally the longest part of the convolution.
    for (size_t n = coreSize - 1; n < inputSize - 1; ++n)
    {
        resultPtr[n] = 0;

        const size_t kmin = n - (coreSize - 1);
        const size_t kmax = n;

        for (size_t k = kmin; k <= kmax; k++)
        {
            resultPtr[n] += inputPtr[k] * corePtr[n - k];
        }
    }

    // The end
    for (size_t n = inputSize - 1; n < inputSize + coreSize - 1; ++n)
    {
        resultPtr[n] = 0;

        const size_t kmin = n - (coreSize - 1);
        const size_t kmax = n;

        for (size_t k = kmin; k <= std::min(kmax, inputSize - 1); k++)
        {
            resultPtr[n] += inputPtr[k] * corePtr[n - k];
        }
    }
}

/* Performs a convolution using FFT between the input and core and stores the result in 'result'.
    The result will have length equal to (input.size() + core.size() - 1)
*/
void ConvolutionCoreFft(const std::vector<double>& input, const std::vector<double>& core, std::vector<double>& result)
{
    const size_t outputSize = input.size() + core.size() - 1;
    const size_t fftSize = (outputSize % 2 == 0) ? outputSize : outputSize + 1; // the FFT_real used below only accepts vectors of even length

    // The convolution is performed by taking the fft of both the input and the core, multiplying their (complex) outputs and taking the inverse fft.

    // Do the fft of the input
    std::vector<std::complex<double>> dftOfInput(fftSize);
    {
        // Create a padded input, where the original input has been padded with zeros to the correct length.
        std::vector<double> paddedInput(fftSize, 0.0);
        memcpy(paddedInput.data(), input.data(), input.size() * sizeof(double));
        novac::Fft_Real(paddedInput, dftOfInput);
    }

    // Do the fft of the core
    std::vector<std::complex<double>> dftOfCore(fftSize);
    {
        // Create a padded core, where the original core has been padded with zeros to the correct length.
        std::vector<double> paddedCore(fftSize, 0.0);
        memcpy(paddedCore.data(), core.data(), core.size() * sizeof(double));
        novac::Fft_Real(paddedCore, dftOfCore);
    }

    // Multiply
    std::vector<std::complex<double>> product(fftSize);
    for (size_t ii = 0; ii < fftSize; ++ii)
    {
        product[ii] = dftOfInput[ii] * dftOfCore[ii];
    }

    // Inverse transform
    const bool forwardTransform = false;
    std::vector<std::complex<double>> complexResult(fftSize);
    novac::Fft(product, complexResult, forwardTransform);

    // Extract the result and return, remember to scale with the length of the fft (since the Ifft above doesn't do that).
    result = novac::Abs(complexResult);
    Mult(result, 1.0 / (double)fftSize);

    if (outputSize != fftSize)
    {
        result.resize(outputSize);
    }
}

bool ConvolveReference(
    const std::string& pixelToWavelengthMappingFile,
    const std::string& slfFile,
    const std::string& highResReferenceFile,
    CCrossSectionData& result,
    WavelengthConversion conversion)
{
    try
    {
        CCrossSectionData pixelToWavelengthMapping;
        if (!ReadCrossSectionFile(pixelToWavelengthMappingFile, pixelToWavelengthMapping, true))
        {
            return false;
        }

        CCrossSectionData slf;
        if (!ReadCrossSectionFile(slfFile, slf))
        {
            return false;
        }

        CCrossSectionData highResReference;
        if (!ReadCrossSectionFile(highResReferenceFile, highResReference))
        {
            return false;
        }

        ConvolveReference(pixelToWavelengthMapping.m_waveLength, slf, highResReference, result.m_crossSection, conversion);

        // save the wavelengths as well
        result.m_waveLength = pixelToWavelengthMapping.m_waveLength;

        return true;
    }
    catch (std::exception& e)
    {
        std::cout << "Exception in ConvolveReference: " << e.what() << std::endl;
        return false;
    }
}

void CreateGaussian(double gaussianSigma, double deltaLambda, CCrossSectionData& gaussian)
{
    const double fwhm = gaussianSigma * 2.0 * std::sqrt(2.0 * std::log(2.0));
    const size_t slfSize = std::max(size_t(35), 1 + 2 * ((size_t)(0.0001 + std::round(4.0 * fwhm / deltaLambda)) / 2));
    gaussian.m_waveLength.resize(slfSize);
    gaussian.m_crossSection.resize(slfSize);

    const size_t halfLength = (slfSize / 2);
    const size_t midIdx = halfLength;
    const double s2_inv = 1.0 / (2.0 * gaussianSigma * gaussianSigma);

    for (size_t dI = 0; dI <= halfLength; ++dI)
    {
        const double x = dI * deltaLambda;
        gaussian.m_waveLength[midIdx + dI] = x;
        gaussian.m_waveLength[midIdx - dI] = -x;

        const double value = std::exp(-(x * x) * s2_inv);
        gaussian.m_crossSection[midIdx + dI] = value;
        gaussian.m_crossSection[midIdx - dI] = value;
    }
}

bool ConvolveReferenceGaussian(
    const std::string& pixelToWavelengthMappingFile,
    double gaussianSigma,
    const std::string& highResReferenceFile,
    CCrossSectionData& result,
    WavelengthConversion conversion)
{
    try
    {
        CCrossSectionData pixelToWavelengthMapping;
        if (!ReadCrossSectionFile(pixelToWavelengthMappingFile, pixelToWavelengthMapping, true))
        {
            return false;
        }

        CCrossSectionData highResReference;
        if (!ReadCrossSectionFile(highResReferenceFile, highResReference))
        {
            return false;
        }

        // Generate a gaussian with high enough resolution
        CCrossSectionData slf;
        const double slfDeltaLambda = std::min(0.25 * Resolution(highResReference.m_waveLength), gaussianSigma * 0.125);
        CreateGaussian(gaussianSigma, slfDeltaLambda, slf);

        ConvolveReference(pixelToWavelengthMapping.m_waveLength, slf, highResReference, result.m_crossSection, conversion);

        // save the wavelengths as well
        result.m_waveLength = pixelToWavelengthMapping.m_waveLength;

        return true;
    }
    catch (std::exception& e)
    {
        std::cout << "Exception in ConvolveReferenceGaussian: " << e.what() << std::endl;
        return false;
    }
}

void Convert(const CCrossSectionData& highResReference, WavelengthConversion conversion, CCrossSectionData& result)
{
    result = highResReference;

    if (conversion == WavelengthConversion::VacuumToAir)
    {
        for (size_t ii = 0; ii < highResReference.m_waveLength.size(); ++ii)
        {
            result.m_waveLength[ii] = NmVacuumToNmAir(highResReference.m_waveLength[ii]);
        }
    }
}

bool Convolve(
    const CCrossSectionData& slf,
    const CCrossSectionData& highResReference,
    std::vector<double>& result,
    WavelengthConversion conversion)
{
    if (slf.m_waveLength.size() != slf.m_crossSection.size())
    {
        std::cout << " Error in call to 'Convolve', the SLF must have as many values as wavelength values." << std::endl;
        return false;
    }
    if (highResReference.m_waveLength.size() != highResReference.m_crossSection.size())
    {
        std::cout << " Error in call to 'Convolve', the reference must have as many values as wavelength values." << std::endl;
        return false;
    }

    // If desired, convert the high-res reference from vacuum to air
    CCrossSectionData convertedHighResReference;
    Convert(highResReference, conversion, convertedHighResReference);

    // We need to make sure we work on the correct resolution, the highest possible to get the most accurate results.
    const double highestResolution = std::min(Resolution(convertedHighResReference.m_waveLength), Resolution(slf.m_waveLength));

    // Resample the highResReference to be on a uniform grid.
    // Study the SLF and the reference to see which has the highest resolution and take that.
    UniformGrid highResGrid;
    highResGrid.minValue = convertedHighResReference.m_waveLength.front();
    highResGrid.maxValue = convertedHighResReference.m_waveLength.back();
    highResGrid.length = (size_t)((highResGrid.maxValue - highResGrid.minValue) / highestResolution);

    std::vector<double> uniformHighResReference;
    Resample(convertedHighResReference, highResGrid.Resolution(), uniformHighResReference);

    // We also need to resample the slit-function to be on the same wavelength-grid as the high-res reference.
    std::vector<double> uniformSlf;
    Resample(slf, highResGrid.Resolution(), uniformSlf);

    // To preserve the energy, we need to normalize the slit-function to the range [0->1]
    std::vector<double> normalizedUniformSlf;
    Normalize(uniformSlf, normalizedUniformSlf);
    assert(normalizedUniformSlf.size() == uniformSlf.size());

    const size_t refSize = uniformHighResReference.size();
    const size_t coreSize = normalizedUniformSlf.size();

    // Do the actual convolution
    std::vector<double> intermediate;
    ConvolutionCore(uniformHighResReference, normalizedUniformSlf, intermediate);

    // Cut down the result to the same size as the output is supposed to be
    CCrossSectionData resultSpec;
    resultSpec.m_crossSection = std::vector<double>(begin(intermediate) + coreSize / 2, begin(intermediate) + refSize + coreSize / 2);
    resultSpec.m_waveLength = std::vector<double>(resultSpec.m_crossSection.size());
    for (size_t ii = 0; ii < resultSpec.m_waveLength.size(); ++ii)
    {
        resultSpec.m_waveLength[ii] = highResGrid.minValue + (highResGrid.maxValue - highResGrid.minValue) * (double)ii / (resultSpec.m_waveLength.size() - 1);
    }

    Resample(resultSpec, Resolution(convertedHighResReference.m_waveLength), result);

    return true;
}

// Small helper method, returns the remainder after as many
//  factors of 2, 3 and 5 as possible have been removed.
size_t GetNonReduciblePrime(size_t number)
{
    while (number % 2 == 0)
    {
        number /= 2;
    }
    while (number % 3 == 0)
    {
        number /= 3;
    }
    while (number % 5 == 0)
    {
        number /= 5;
    }
    return number;
}

void ConvolveReference(
    const std::vector<double>& pixelToWavelengthMapping,
    const CCrossSectionData& slf,
    const CCrossSectionData& highResReference,
    std::vector<double>& result,
    WavelengthConversion conversion,
    ConvolutionMethod method,
    double fwhmOfInstrumentLineShape,
    bool normalizeSlf)
{
    if (slf.m_waveLength.size() != slf.m_crossSection.size())
    {
        throw std::invalid_argument(" Error in call to 'ConvolveReference', the SLF must have as many values as wavelength values.");
    }
    if (highResReference.m_waveLength.size() != highResReference.m_crossSection.size())
    {
        throw std::invalid_argument(" Error in call to 'ConvolveReference', the reference must have as many values as wavelength values.");
    }

    // If desired, convert the high-res reference from vacuum to air
    CCrossSectionData convertedHighResReference;
    Convert(highResReference, conversion, convertedHighResReference);

    // We need to make sure we work on the correct resolution, the highest possible to get the most accurate results.
    if (fwhmOfInstrumentLineShape < std::numeric_limits<float>::epsilon())
    {
        fwhmOfInstrumentLineShape = GetFwhm(slf);
    }
    if (fwhmOfInstrumentLineShape < std::numeric_limits<float>::epsilon())
    {
        throw std::invalid_argument(" Error in call to 'ConvolveReference', the estimated fwhm of the instrument line shape is zero.");
    }
    const double minimumAllowedResolution = 0.02 * fwhmOfInstrumentLineShape; // do use at least 50 points per FWHM of the SLF
    const double maximumAllowedResolution = 0.01 * fwhmOfInstrumentLineShape; // do not use more than 100 points per FWHM of the SLF
    const double resolutionOfReference = Resolution(convertedHighResReference.m_waveLength);
    const double highestResolution = std::max(std::min(resolutionOfReference, maximumAllowedResolution), minimumAllowedResolution);

    // Resample the convertedHighResReference to be on a uniform grid.
    // Study the SLF and the reference to see which has the highest resolution and take that.
    UniformGrid convolutionGrid;
    convolutionGrid.minValue = convertedHighResReference.m_waveLength.front();
    convolutionGrid.maxValue = convertedHighResReference.m_waveLength.back();

    if (method == ConvolutionMethod::Fft)
    {
        // For the sake of efficiency of the fft below, find the length which contain as many factors of two, tree and five as possible
        //  while still being between minimumAllowedResolution and maximumAllowedResolution.
        // This is done by finding the smallest length for which the remainder after removing all 2, 3 and 5 prime factors is minimal (usually one).
        const size_t maximumLength = 2 * (size_t)((convolutionGrid.maxValue - convolutionGrid.minValue) / maximumAllowedResolution);

        size_t testLengthForConvolutionGrid = 2 * (size_t)((convolutionGrid.maxValue - convolutionGrid.minValue) / highestResolution);
        size_t optimumConvolutionGridLength = testLengthForConvolutionGrid;
        convolutionGrid.length = testLengthForConvolutionGrid;
        size_t slfLength = 1 + (size_t)(std::round((slf.m_waveLength.back() - slf.m_waveLength.front()) / convolutionGrid.Resolution()));
        size_t fftLength = testLengthForConvolutionGrid + slfLength - 1; // take into account that what we're calculating the fft of isn't just the input length but the input + slf - 1.
        size_t optimumRemainder = GetNonReduciblePrime(fftLength);

        while (testLengthForConvolutionGrid < maximumLength)
        {
            ++testLengthForConvolutionGrid;

            convolutionGrid.length = testLengthForConvolutionGrid;
            slfLength = 1 + (size_t)(std::round((slf.m_waveLength.back() - slf.m_waveLength.front()) / convolutionGrid.Resolution()));
            fftLength = testLengthForConvolutionGrid + slfLength - 1;

            const size_t remainder = GetNonReduciblePrime(fftLength);

            if (remainder < optimumRemainder)
            {
                optimumRemainder = remainder;
                optimumConvolutionGridLength = testLengthForConvolutionGrid;
            }
        }
        convolutionGrid.length = optimumConvolutionGridLength;

        // std::cout << "Convolution is using length: " << optimumConvolutionGridLength << " with smallest remainder: " << optimumRemainder << std::endl;
    }
    else
    {
        // Use the highest resolution
        convolutionGrid.length = 2 * (size_t)((convolutionGrid.maxValue - convolutionGrid.minValue) / highestResolution);
    }

    std::vector<double> uniformHighResReference;
    Resample(convertedHighResReference, convolutionGrid.Resolution(), uniformHighResReference);
    assert(uniformHighResReference.size() == convolutionGrid.length);

    // We also need to resample the slit-function to be on the same wavelength-grid as the high-res reference.
    std::vector<double> resampledSlf;
    Resample(slf, convolutionGrid.Resolution(), resampledSlf);

    // To preserve the energy, we need to normalize the slit-function to the range [0->1]
    std::vector<double> normalizedSlf;
    if (normalizeSlf)
    {
        NormalizeArea(resampledSlf, normalizedSlf);
        assert(normalizedSlf.size() == resampledSlf.size());
    }
    else
    {
        normalizedSlf = resampledSlf;
    }

    const size_t refSize = uniformHighResReference.size();
    const size_t coreSize = normalizedSlf.size();

    // Do the actual convolution
    std::vector<double> intermediate;
    if (method == ConvolutionMethod::Direct)
    {
        ConvolutionCore(uniformHighResReference, normalizedSlf, intermediate);
    }
    else if (method == ConvolutionMethod::Fft)
    {
        ConvolutionCoreFft(uniformHighResReference, normalizedSlf, intermediate);
    }

    CCrossSectionData resultSpec;

    // Cut down the result to the same size as the output is supposed to be
    resultSpec.m_crossSection = std::vector<double>(begin(intermediate) + coreSize / 2, begin(intermediate) + refSize + coreSize / 2);

    resultSpec.m_waveLength = std::vector<double>(resultSpec.m_crossSection.size());
    convolutionGrid.Generate(resultSpec.m_waveLength);

    Resample(resultSpec, pixelToWavelengthMapping, result);
}

bool ConvolveReference_Fast(
    const std::vector<double>& pixelToWavelengthMapping,
    const CCrossSectionData& slf,
    const CCrossSectionData& highResReference,
    std::vector<double>& result,
    WavelengthConversion conversion)
{
    if (slf.m_waveLength.size() != slf.m_crossSection.size())
    {
        std::cout << " Error in call to 'ConvolveReference', the SLF must have as many values as wavelength values." << std::endl;
        return false;
    }
    if (highResReference.m_waveLength.size() != highResReference.m_crossSection.size())
    {
        std::cout << " Error in call to 'ConvolveReference', the reference must have as many values as wavelength values." << std::endl;
        return false;
    }

    // If desired, convert the high-res reference from vacuum to air
    CCrossSectionData convertedHighResReference;
    Convert(highResReference, conversion, convertedHighResReference);

    // We need to make sure we work on the correct resolution, the highest possible to get the most accurate results.
    const double highestResolution = std::min(Resolution(convertedHighResReference.m_waveLength), Resolution(slf.m_waveLength));

    // Resample the convertedHighResReference to be on a uniform grid.
    // Study the SLF and the reference to see which has the highest resolution and take that.
    UniformGrid convolutionGrid;
    convolutionGrid.minValue = convertedHighResReference.m_waveLength.front();
    convolutionGrid.maxValue = convertedHighResReference.m_waveLength.back();
    convolutionGrid.length = 2 * (size_t)((convolutionGrid.maxValue - convolutionGrid.minValue) / highestResolution);

    std::vector<double> uniformHighResReference;
    Resample(convertedHighResReference, convolutionGrid.Resolution(), uniformHighResReference);
    assert(uniformHighResReference.size() == convolutionGrid.length);

    // We also need to resample the slit-function to be on the same wavelength-grid as the high-res reference.
    std::vector<double> resampledSlf;
    Resample(slf, convolutionGrid.Resolution(), resampledSlf);

    // To preserve the energy, we need to normalize the slit-function to the range [0->1]
    std::vector<double> normalizedSlf;
    NormalizeArea(resampledSlf, normalizedSlf);
    assert(normalizedSlf.size() == resampledSlf.size());

    const size_t refSize = uniformHighResReference.size();
    const size_t coreSize = normalizedSlf.size();

    // Do the actual convolution
    std::vector<double> intermediate;
    ConvolutionCore(uniformHighResReference, normalizedSlf, intermediate);

    // Cut down the result to the same size as the output is supposed to be
    CCrossSectionData resultSpec;
    resultSpec.m_crossSection = std::vector<double>(begin(intermediate) + coreSize / 2, begin(intermediate) + refSize + coreSize / 2);
    resultSpec.m_waveLength = std::vector<double>(resultSpec.m_crossSection.size());
    convolutionGrid.Generate(resultSpec.m_waveLength);

    Resample(resultSpec, pixelToWavelengthMapping, result);

    return true;
}


}

