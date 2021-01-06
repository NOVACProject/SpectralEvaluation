#include <SpectralEvaluation/Math/FFT.h>

#pragma warning(push)
#pragma warning(disable:4244)
#pragma warning(disable:4267)

// run the FFT using double (not float as default)
#define kiss_fft_scalar double

#include "kissfft/kiss_fft.h"
#include "kissfft/kiss_fft.c"
#include "kissfft/kiss_fftr.h"
#include "kissfft/kiss_fftr.c"

namespace novac
{

template<class T>
std::vector<T> RealT(const std::vector<std::complex<T>>& complexData)
{
    std::vector<T> result(complexData.size());
    for (size_t ii = 0; ii < complexData.size(); ++ii)
    {
        result[ii] = complexData[ii].real();
    }
    return result;
}

template<class T>
std::vector<T> ImagT(const std::vector<std::complex<T>>& complexData)
{
    std::vector<T> result(complexData.size());
    for (size_t ii = 0; ii < complexData.size(); ++ii)
    {
        result[ii] = complexData[ii].imag();
    }
    return result;
}

std::vector<float> Real(const std::vector<std::complex<float>>& complexData)
{
    return RealT(complexData);
}
std::vector<double> Real(const std::vector<std::complex<double>>& complexData)
{
    return RealT(complexData);
}

std::vector<float> Imag(const std::vector<std::complex<float>>& complexData)
{
    return ImagT(complexData);
}
std::vector<double> Imag(const std::vector<std::complex<double>>& complexData)
{
    return ImagT(complexData);
}


// The FFT implementation here makes use of kiss fft https://github.com/mborgerding/kissfft

void Fft(const std::vector<std::complex<double>>& input, std::vector<std::complex<double>>& result, bool forward)
{
    const int isInverse = (forward) ? 0 : 1;
    const size_t length = input.size();

    kiss_fft_cfg cfg = kiss_fft_alloc(static_cast<int>(length), isInverse, nullptr, nullptr);
    if (cfg != nullptr)
    {
        // Copy the input data
        std::vector<kiss_fft_cpx> cx_in(length);
        for (size_t ii = 0; ii < length; ++ii)
        {
            cx_in[ii].r = static_cast<kiss_fft_scalar>(input[ii].real());
            cx_in[ii].i = static_cast<kiss_fft_scalar>(input[ii].imag());
        }

        // Perform the FFT
        std::vector<kiss_fft_cpx> cx_out(length);
        kiss_fft(cfg, cx_in.data(), cx_out.data());

        // Copy out the result
        result.resize(length);
        for (size_t ii = 0; ii < length; ++ii)
        {
            result[ii] = std::complex<double>(static_cast<double>(cx_out[ii].r), static_cast<double>(cx_out[ii].i));
        }

        kiss_fft_free(cfg);
    }
}

void Fft_Real(const std::vector<double>& input, std::vector<std::complex<double>>& result, bool forward)
{
    const int isInverse = (forward) ? 0 : 1;
    const size_t length = input.size();

    kiss_fftr_cfg cfg = kiss_fftr_alloc(static_cast<int>(length), isInverse, nullptr, nullptr);
    if (cfg != nullptr)
    {
        // Perform the FFT
        kiss_fft_cpx defaultValue{};
        std::vector<kiss_fft_cpx> cx_out(length, defaultValue);
        kiss_fftr(cfg, input.data(), cx_out.data());

        // Copy out the result, the output is symmetrical around the mid point and kiss_fftr only fills in the first half of the values
        //  the second loop here is to fill in the second half (with changed sign on the imaginary component).
        result.resize(length);
        const size_t midLength = length / 2;
        for (size_t ii = 0; ii < midLength; ++ii)
        {
            result[ii] = std::complex<double>(static_cast<double>(cx_out[ii].r), static_cast<double>(cx_out[ii].i));
        }
        for (size_t ii = midLength; ii < length; ++ii)
        {
            result[ii] = std::complex<double>(static_cast<double>(cx_out[length - ii].r), static_cast<double>(-cx_out[length - ii].i));
        }

        kiss_fft_free(cfg);
    }
}

}

#pragma warning(pop)

