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
    const WavelengthCalibrationSetup& calibrationSetup,
    const SpectrometerCalibration& initialCalibration,
    const CSpectrum& measuredspectrum,
    SpectrometerCalibration& /*result*/)
{
    ::Evaluation::CCrossSectionData measuredCopy;
    measuredCopy.m_waveLength = initialCalibration.wavelengthToPixelMapping;
    measuredCopy.m_crossSection = std::vector<double>(measuredspectrum.m_data, measuredspectrum.m_data + measuredspectrum.m_length);

    // TODO: we need a way to pass in the wavelength range as parameter here
    double lambdaLow = 310.0;
    double lambdaHigh = 320.0;
    ::Evaluation::CCrossSectionData localSolarAtlas{ calibrationSetup.solarAtlas, lambdaLow, lambdaHigh };

    // Extract the measured spectrum in this local region, saves processing time...
    ::Evaluation::CCrossSectionData localMeasured{ measuredCopy, lambdaLow, lambdaHigh };

    // Convolve the solar atlas and the references
    ::Evaluation::CCrossSectionData convolvedSolarAtlas;
    if (!ConvolveReference(localMeasured.m_waveLength, initialCalibration.slf, localSolarAtlas, convolvedSolarAtlas.m_crossSection))
    {
        return false;
    }
    convolvedSolarAtlas.m_waveLength = localMeasured.m_waveLength;

    std::vector<::Evaluation::CCrossSectionData> convolvedReference;
    convolvedReference.reserve(calibrationSetup.crossSections.size());
    for (size_t refIdx = 0; refIdx < calibrationSetup.crossSections.size(); ++refIdx)
    {
        ::Evaluation::CCrossSectionData reference;
        if (!ConvolveReference(localMeasured.m_waveLength, initialCalibration.slf, calibrationSetup.crossSections[refIdx], reference.m_crossSection))
        {
            return false;
        }
        reference.m_waveLength = localMeasured.m_waveLength;

        convolvedReference.push_back(reference);
    }

    // Setup a DOAS fit here where we fit the kurucz + high-res references against the measured spectrum.
    ::Evaluation::CFitWindow doasFitWindow;
    doasFitWindow.fitLow = 0; // TODO: setup properly
    doasFitWindow.fitHigh = localMeasured.GetSize(); // TODO: setup properly
    doasFitWindow.fitType = ::Evaluation::FIT_POLY;
    doasFitWindow.polyOrder = 2;
    doasFitWindow.includeIntensitySpacePolyominal = true;
    doasFitWindow.ringCalculation = ::Evaluation::RING_CALCULATION_OPTION::CALCULATE_RING_X2;
    doasFitWindow.specLength = (int)localMeasured.m_crossSection.size();
    doasFitWindow.shiftSky = true;
    for (size_t ii = 0; ii < convolvedReference.size(); ++ii)
    {
        doasFitWindow.ref[ii] = ::Evaluation::CReferenceFile(convolvedReference[ii]);
    }
    doasFitWindow.nRef = (int)convolvedReference.size();

    /*
    CWavelengthFit doasFit{ doasFitWindow };
    doasFit.SetSkySpectrum(convolvedSolarAtlas);
    doasFit.Evaluate(localMeasured.m_crossSection.data(), localMeasured.m_crossSection.size());
    */

    return true;
}
}
