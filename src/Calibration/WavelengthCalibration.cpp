#include <algorithm>
#include <SpectralEvaluation/Calibration/WavelengthCalibration.h>
#include <SpectralEvaluation/Calibration/ReferenceSpectrumConvolution.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Evaluation/WavelengthFit.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/VectorUtils.h>
#include <SpectralEvaluation/File/File.h>

namespace Calibration
{

bool MercuryCalibration(
    const CSpectrum& measuredMercurySpectrum,
    int polynomialOrder,
    double minimumWavelength,
    double maximumWavelength,
    std::vector<SpectrumDataPoint>& /*foundPeaks*/,
    std::vector<double> pixelToWavelengthPolynomial)
{
    if (measuredMercurySpectrum.m_length < 50)
    {
        return false;
    }
    else if (polynomialOrder < 0 || polynomialOrder > 5)
    {
        return false;
    }
    else if (maximumWavelength < minimumWavelength)
    {
        return false;
    }

    // List of known mercury lines, in nm(air)
    const std::vector<double> knownMercuryLines = { 284.7675, 302.1498, 312.5668, 313.1548, 313.1839, 365.0153, 365.4836, 366.3279, 398.3931, 404.6563, 435.8328 };



    return false;
}


bool EstimateWavelengthToPixelMapping(
    const WavelengthCalibrationSetup& /*calibrationSetup*/,
    const SpectrometerCalibration& /*initialCalibration*/,
    const CSpectrum& /*measuredspectrum*/,
    SpectrometerCalibration& /*result*/)
{
 


    return true;
}
}
