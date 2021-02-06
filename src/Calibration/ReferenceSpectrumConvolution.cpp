#include <SpectralEvaluation/Calibration/ReferenceSpectrumConvolution.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Fit/CubicSplineFunction.h>
#include <SpectralEvaluation/File/File.h>
#include <SpectralEvaluation/VectorUtils.h>
#include <SpectralEvaluation/Spectra/Grid.h>
#include <SpectralEvaluation/Air.h>
#include <SpectralEvaluation/Math/FFT.h>
#include <iostream>
#include <assert.h>

namespace Evaluation
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

// Resamples the provided slf to the provided wavelength grid.
void Resample(const CCrossSectionData& slf, const std::vector<double>& wavelength, std::vector<double>& resampledSlf)
{
    const double xMin = slf.m_waveLength.front();
    const double xMax = slf.m_waveLength.back();

    std::vector<double> xCopy(begin(slf.m_waveLength), end(slf.m_waveLength)); // a non-const local copy
    std::vector<double> yCopy(begin(slf.m_crossSection), end(slf.m_crossSection)); // a non-const local copy

    MathFit::CVector slfX(xCopy.data(), (int)xCopy.size(), 1, false);
    MathFit::CVector slfY(yCopy.data(), (int)yCopy.size(), 1, false);

    // Create a spline from the slit-function.
    MathFit::CCubicSplineFunction spline(slfX, slfY);

    // do the resampling...
    resampledSlf.resize(wavelength.size());
    for (size_t ii = 0; ii < wavelength.size(); ++ii)
    {
        if (wavelength[ii] >= xMin && wavelength[ii] <= xMax)
        {
            resampledSlf[ii] = spline.GetValue(wavelength[ii]);
        }
        else
        {
            resampledSlf[ii] = 0.0;
        }
    }
}

/* Performs a convolution between the input and core and stores the result in 'result' */
void ConvolutionCore(const std::vector<double>& input, const std::vector<double>& core, std::vector<double>& result)
{
    const size_t refSize = input.size();
    const size_t coreSize = core.size();

    result.resize(refSize + coreSize - 1, 0.0);

    // The actual convolution. Here a dead-simple raw convolution calculation. This can be made faster using FFT if required.
    //  Get the pointers to the data, this avoids range-checking for each and every value (at least in debug mode)
    const double* inputPtr = input.data();
    const double* corePtr = core.data();
    double* resultPtr = result.data();
    for (size_t n = 0; n < refSize + coreSize - 1; ++n)
    {
        result[n] = 0;

        const size_t kmin = (n >= coreSize - 1) ? n - (coreSize - 1) : 0;
        const size_t kmax = (n < refSize - 1) ? n : refSize - 1;

        for (size_t k = kmin; k <= kmax; k++)
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
    CCrossSectionData pixelToWavelengthMapping;
    if (!FileIo::ReadCrossSectionFile(pixelToWavelengthMappingFile, pixelToWavelengthMapping, true))
    {
        return false;
    }

    CCrossSectionData slf;
    if (!FileIo::ReadCrossSectionFile(slfFile, slf))
    {
        return false;
    }

    CCrossSectionData highResReference;
    if (!FileIo::ReadCrossSectionFile(highResReferenceFile, highResReference))
    {
        return false;
    }

    if (!ConvolveReference(pixelToWavelengthMapping.m_waveLength, slf, highResReference, result.m_crossSection, conversion))
    {
        return false;
    }

    // save the wavelengths as well
    result.m_waveLength = pixelToWavelengthMapping.m_waveLength;

    return true;
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
    CCrossSectionData pixelToWavelengthMapping;
    if (!FileIo::ReadCrossSectionFile(pixelToWavelengthMappingFile, pixelToWavelengthMapping, true))
    {
        return false;
    }

    CCrossSectionData highResReference;
    if (!FileIo::ReadCrossSectionFile(highResReferenceFile, highResReference))
    {
        return false;
    }

    // Generate a gaussian with high enough resolution
    CCrossSectionData slf;
    const double slfDeltaLambda = std::min(0.25 * Resolution(highResReference.m_waveLength), gaussianSigma * 0.125);
    CreateGaussian(gaussianSigma, slfDeltaLambda, slf);

    if (!ConvolveReference(pixelToWavelengthMapping.m_waveLength, slf, highResReference, result.m_crossSection, conversion))
    {
        return false;
    }

    // save the wavelengths as well
    result.m_waveLength = pixelToWavelengthMapping.m_waveLength;

    return true;
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

bool ConvolveReference(
    const std::vector<double>& pixelToWavelengthMapping,
    const CCrossSectionData& slf,
    const CCrossSectionData& highResReference,
    std::vector<double>& result,
    WavelengthConversion conversion,
    ConvolutionMethod method)
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
    const double resolutionOfSlf = Resolution(slf.m_waveLength);
    const double maximumAllowedResolution = 0.2 * resolutionOfSlf; // do not use more than 5 points per pixel of the SLF
    const double resolutionOfReference = Resolution(convertedHighResReference.m_waveLength);
    const double highestResolution = std::max(std::min(resolutionOfReference, resolutionOfSlf), maximumAllowedResolution);

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

    return true;
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

