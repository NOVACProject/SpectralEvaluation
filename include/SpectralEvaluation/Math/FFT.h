#pragma once

#include <vector>
#include <complex>

// ---------------------------------------------------------------------
// ----------------- CALCULATING THE FOURIER TRANSFORM -----------------
// ---------------------------------------------------------------------

namespace novac
{

/** Returns the real value component of the provided vector of complex values */
std::vector<float> Real(const std::vector<std::complex<float>>& complexData);
std::vector<double> Real(const std::vector<std::complex<double>>& complexData);

/** Returns the imaginary value component of the provided vector of complex values */
std::vector<float> Imag(const std::vector<std::complex<float>>& complexData);
std::vector<double> Imag(const std::vector<std::complex<double>>& complexData);

/** Calculates the Fourier Transform of the provided complex-valued input sequence
    and stores the result in the given output parameter 'result'.
    @param input The data to calculate the fourier transform of.
    @param result Will on scucessful return be filled with the calculated fourier transform data of the input.
        This will be resized to input.size() if required.
    @param forward If true then the forward Fourier transform is calculated, if false then the inverse Fourier transform is calculated. 
    The length of the input does not have to be an even power of two, but the algorithm is faster if it is. */
void Fft(const std::vector<std::complex<double>>& input, std::vector<std::complex<double>>& result, bool forward = true);

/** Calculates the Fourier Transform of the provided real-valued input sequence
    and stores the result in the given output parameter 'result'.
    This is considerably faster than using the more general method 'Fft' when the input is real valued.
    @param input The data to calculate the fourier transform of.
    @param result Will on scucessful return be filled with the calculated fourier transform data of the input.
        This will be resized to input.size() if required.
    @param forward If true then the forward Fourier transform is calculated, if false then the inverse Fourier transform is calculated.
    The length of the input MUST be an even number. */
void Fft_Real(const std::vector<double>& input, std::vector<std::complex<double>>& result, bool forward = true);


}