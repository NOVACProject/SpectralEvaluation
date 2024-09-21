#pragma once

#include <vector>

namespace novac
{
    class CSpectrum;
    class IScanSpectrumSource;

    /** Calculates the average or sum of the spectra with the given indices in the given scan.
        @param scan The source of the spectra from one scan.
        @param indices The (zero based) indices of the spectra in the scan to extract and average.
        @param average If true then the result is averaged, if false then the result is added.
        @return the number of spectra averaged/added. */
    int AverageSpectra(IScanSpectrumSource& scan, const std::vector<int>& indices, CSpectrum& result, bool average = true);

}