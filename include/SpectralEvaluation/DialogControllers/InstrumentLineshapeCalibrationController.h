#pragma once

#include <map>
#include <memory>
#include <vector>
#include <string>
#include <SpectralEvaluation/Spectra/SpectrumUtils.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Calibration/InstrumentCalibration.h>

class InstrumentLineshapeCalibrationController
{
public:
    InstrumentLineshapeCalibrationController();

    enum class LineShapeFunction
    {
        None,
        Gaussian,
        SuperGauss
    };

    /** Input: The full path to the mercury spectrum to retreive the line shape from. */
    std::string m_inputSpectrumPath;

    /** Input: The full path to the dark to use (if any) */
    std::string m_darkSpectrumPath;

    /** Input: set to true if we should read the wavelength calibration
        from file (if available). Set to false to attempt to auto-determine
        the wavelength calibration from the mercury emission lines. */
    bool m_readWavelengthCalibrationFromFile = true;

    /** Output: this is our list of found peaks in the m_inputSpectrum. */
    std::vector<novac::SpectrumDataPoint> m_peaksFound;

    /** Output: this is our list of peaks which have been found in the m_inputSpectrum
        but have been rejected for some reason */
    std::vector<novac::SpectrumDataPoint> m_rejectedPeaks;

    /** This is our best estimate of the instrument calibration.
        The pixel-to-wavelength mapping can either be read from file or by
        performing a pixel-to-wavelength calibration based on the location of found mercury lines. */
    std::unique_ptr<novac::InstrumentCalibration> m_resultingCalibration;

    /** Output: The read in mercury spectrum data */
    std::vector<double> m_inputSpectrum;

    /** Output: This is set to true (when reading the input spectrum)
        if the read in spectrum already contains a pixel-to-wavelength calibration.
        Otherwise false. */
    bool m_inputSpectrumContainsWavelength = false;

    /** Output: this is set to true if we are to perform a wavelength calibration ourselves
        (i.e. if m_readWavelengthCalibrationFromFile==false or m_inputSpectrumContainsWavelength==false)
        and the wavelength calibration succeeded. */
    bool m_wavelengthCalibrationSucceeded = false;

    /** Locates peaks in the spectrum.
        This updates m_inputSpectrum and m_peaksFound. */
    void Update();

    /** Fits a function to the peak with the provided index.
        This updates m_fittedLineShape and m_sampledLineShapeFunction */
    void FitFunctionToLineShape(size_t peakIdx, LineShapeFunction function);

    /** Extracts the peak with the provided index from the read in spectrum file.
        The peak will be normalized and have its baseline subtracted. */
    std::unique_ptr<novac::CSpectrum> ExtractSelectedMercuryPeak(size_t peakIdx, double* baseline = nullptr, size_t* startIdx = nullptr) const;

    /** Returns a textual summary of parameters / properties of the last fitted function
        (i.e. the last call to FitFunctionToLineShape).
        Returns an empty set if no function has been fitted. */
    std::vector<std::pair<std::string, std::string>> GetFittedFunctionDescription() const;

    /** Saves the resulting instrument line shape information as a .std file.
        Notice that both 'Update' and 'FitFunctionToLineShape' must have been called prior to this. */
    void SaveResultAsStd(size_t peakIdx, const std::string& filename);

    /** Saves the resulting pixel-to-wavelength mapping information as a .clb file. */
    void SaveResultAsClb(const std::string& filename);

    /** Saves the resulting instrument line shape information as a .slf file. */
    void SaveResultAsSlf(size_t peakIdx, const std::string& filename);

private:
    void ClearFittedLineShape();

    /** Subtracts the baseline from the provided spectrum.
     *  @param spectrum The spectrum to modify.
     *  @return The baseline which was subtracted */
    static double SubtractBaseline(novac::CSpectrum& spectrum);

    /** The meta-data regarding the last read in mercury spectrum. */
    novac::CSpectrumInfo m_inputspectrumInformation;

};