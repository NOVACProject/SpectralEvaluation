#pragma once

#include <memory>
#include <vector>
#include <string>

// ---------------------------------------------------------------------------------------------------------------
// --- This header contains CrossSectionSpectrumGenerator which helps with convolving cross section references ---
// ---------------------------------------------------------------------------------------------------------------

namespace novac
{
class CSpectrum;
class CCrossSectionData;
struct WavelengthRange;

class ICrossSectionSpectrumGenerator
{
public:

    /// <summary>
    /// Returns the wavelength range over which subsequent calls to 'GetCrossSection' with the provided
    /// pixel-to-wavelength calibration will be valid. This is determined both by the provided range and by the
    /// range of the included high resolution absorption cross section.
    /// </summary>
    virtual WavelengthRange GetSpectrumRange(const std::vector<double>& wavelengthCalibration) = 0;

    /// <summary>
    /// Creates a cross section reference spectrum using the provided pixel-to-wavelength mapping and measured instrument line shape.
    /// This will determine the fwhm of the provided instrument line shape and use this value to determine the convolution grid.
    /// </summary>
    /// <param name="pixelToWavelengthMapping">The wavelength (in nm air) for each pixel on the detector.</param>
    /// <param name="measuredInstrumentLineShape">A measurement of the instrument line shape</param>
    /// <returns>The high resolution absorption cross section spectrum convolved with the measured slf and resample to the provided grid.</returns>
    virtual std::unique_ptr<CSpectrum> GetCrossSection(
        const std::vector<double>& wavelengthCalibration,
        const novac::CCrossSectionData& measuredInstrumentLineShape) = 0;

    /// <summary>
    /// Creates a cross section reference spectrum using the provided pixel-to-wavelength mapping and measured instrument line shape.
    /// This will determine the fwhm of the provided instrument line shape and use this value to determine the convolution grid.
    /// </summary>
    /// <param name="pixelToWavelengthMapping">The wavelength (in nm air) for each pixel on the detector.</param>
    /// <param name="measuredInstrumentLineShape">A measurement of the instrument line shape</param>
    /// <returns>The high resolution absorption cross section spectrum convolved with the measured slf and resample to the provided grid.</returns>
    virtual std::unique_ptr<CSpectrum> GetCrossSection(
        const std::vector<double>& wavelengthCalibration,
        const novac::CCrossSectionData& measuredInstrumentLineShape,
        double fwhmOfInstrumentLineShape,
        bool normalize) = 0;
};


class CrossSectionSpectrumGenerator : public ICrossSectionSpectrumGenerator
{
public:
    CrossSectionSpectrumGenerator(const std::string& highResolutionCrossSection)
        : m_crossSectionFile(highResolutionCrossSection)
    {
    }

    /// <summary>
    /// Returns the wavelength range over which subsequent calls to 'GetCrossSection' with the provided
    /// pixel-to-wavelength calibration will be valid. This is determined both by the provided range and by the
    /// range of the included high resolution absorption cross section.
    /// </summary>
    virtual WavelengthRange GetSpectrumRange(const std::vector<double>& wavelengthCalibration) override;

    /// <summary>
    /// Creates a cross section reference spectrum using the provided pixel-to-wavelength mapping and measured instrument line shape.
    /// This will determine the fwhm of the provided instrument line shape and use this value to determine the convolution grid.
    /// </summary>
    /// <param name="pixelToWavelengthMapping">The wavelength (in nm air) for each pixel on the detector.</param>
    /// <param name="measuredInstrumentLineShape">A measurement of the instrument line shape</param>
    /// <returns>The high resolution absorption cross section spectrum convolved with the measured slf and resample to the provided grid.</returns>
    /// <throws>std::invalid_argument if the measuredInstrumentLineShape does not have a valid wavelength calibration 
    ///     or the fwhm could not be determined from the instrument line shape. </throws>
    virtual std::unique_ptr<CSpectrum> GetCrossSection(
        const std::vector<double>& wavelengthCalibration,
        const novac::CCrossSectionData& measuredInstrumentLineShape) override;

    /// <summary>
    /// Creates a cross section reference spectrum using the provided pixel-to-wavelength mapping and measured instrument line shape.
    /// This will determine the fwhm of the provided instrument line shape and use this value to determine the convolution grid.
    /// </summary>
    /// <param name="pixelToWavelengthMapping">The wavelength (in nm air) for each pixel on the detector.</param>
    /// <param name="measuredInstrumentLineShape">A measurement of the instrument line shape</param>
    /// <returns>The high resolution absorption cross section spectrum convolved with the measured slf and resample to the provided grid.</returns>
    /// <throws>std::invalid_argument if the measuredInstrumentLineShape does not have a valid wavelength calibration 
    ///     or the fwhm was not provided and could not be determined from the instrument line shape. </throws>
    virtual std::unique_ptr<CSpectrum> GetCrossSection(
        const std::vector<double>& wavelengthCalibration,
        const novac::CCrossSectionData& measuredInstrumentLineShape,
        double fwhmOfInstrumentLineShape,
        bool normalize) override;

private:

    const std::string m_crossSectionFile;

    std::unique_ptr<novac::CCrossSectionData> m_highResolutionCrossSection;

    void ReadCrossSection();
};

}