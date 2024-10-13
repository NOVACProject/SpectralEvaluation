#include <SpectralEvaluation/Math/SpectrumMath.h>
#include <SpectralEvaluation/Spectra/IScanSpectrumSource.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>

namespace novac
{

int AverageSpectra(IScanSpectrumSource& scan, const std::vector<int>& indices, CSpectrum& result, bool average)
{
    novac::LogContext context; // TODO: Get from input

    if (indices.size() == 0)
    {
        return 0;
    }

    if (scan.GetSpectrum(context, indices[0], result))
    {
        return 0;
    }
    int nofAveragedSpectra = result.NumSpectra();

    for (size_t ii = 1; ii < indices.size(); ++ii)
    {
        CSpectrum tmpSpec;
        if (0 == scan.GetSpectrum(context, indices[ii], tmpSpec))
        {
            result.Add(tmpSpec);
            nofAveragedSpectra += tmpSpec.NumSpectra();
        }
    }

    if (average)
    {
        result.Div((double)nofAveragedSpectra);
    }

    return nofAveragedSpectra;
}
}
