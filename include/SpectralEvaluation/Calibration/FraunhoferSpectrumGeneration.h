#pragma once

#include <memory>
#include <vector>
#include <string>
#include <utility>

// ---------------------------------------------------------------------------------------------------------------
// ---- This header contains FraunhoferSpectrumGeneration which helps with convolving Fraunhofer references  -----
// ---------------------------------------------------------------------------------------------------------------

namespace novac
{

class CSpectrum;
class CCrossSectionData;

class IFraunhoferSpectrumGenerator
{
public:
    /// <summary>
    /// Creates a Fraunhofer reference spectrum using the provided pixel-to-wavelength mapping and measured instrument line shape.
    /// This will determine the fwhm of the provided instrument line shape and use this value to determine the convolution grid.
    /// </summary>
    /// <param name="pixelToWavelengthMapping">The wavelength (in nm air) for each pixel on the detector.</param>
    /// <param name="measuredInstrumentLineShape">A measurement of the instrument line shape</param>
    /// <returns>The high resolution solar spectrum convolved with the measured slf and resample to the provided grid.</returns>
    virtual std::unique_ptr<CSpectrum> GetFraunhoferSpectrum(
        const std::vector<double>& wavelengthCalibration,
        const novac::CCrossSectionData& measuredInstrumentLineShape) = 0;

    /// <summary>
    /// Creates a Fraunhofer reference spectrum using the provided pixel-to-wavelength mapping and measured instrument line shape.
    /// This will determine the fwhm of the provided instrument line shape and use this value to determine the convolution grid.
    /// </summary>
    /// <param name="pixelToWavelengthMapping">The wavelength (in nm air) for each pixel on the detector.</param>
    /// <param name="measuredInstrumentLineShape">A measurement of the instrument line shape</param>
    /// <returns>The high resolution solar spectrum convolved with the measured slf and resample to the provided grid.</returns>
    virtual std::unique_ptr<CSpectrum> GetFraunhoferSpectrum(
        const std::vector<double>& wavelengthCalibration,
        const novac::CCrossSectionData& measuredInstrumentLineShape,
        bool normalize) = 0;

    /// <summary>
    /// Creates a Fraunhofer reference spectrum using the provided pixel-to-wavelength mapping and differential instrument line shape.
    /// </summary>
    /// <param name="pixelToWavelengthMapping">The wavelength (in nm air) for each pixel on the detector.</param>
    /// <param name="measuredInstrumentLineShape">A measurement of the instrument line shape.</param>
    /// <param name="fwhmOfInstrumentLineShape>The Full Width at Half Maximum of the provided instrument line shape.</param>
    /// <returns>The high resolution solar spectrum convolved with the measured slf and resample to the provided grid.</returns>
    virtual std::unique_ptr<CSpectrum> GetDifferentialFraunhoferSpectrum(
        const std::vector<double>& wavelengthCalibration,
        const novac::CCrossSectionData& measuredInstrumentLineShape,
        double fwhmOfInstrumentLineShape) = 0;
};

/// <summary>
/// This is a helper class for generating a Fraunhofer spectrum from a high resolved
/// solar spectrum, a likewise high resolved ozone spectrum and a given instrument setup.
/// Notice that this class will read in the high-resolved solar spectrum when needed (calling GetFraunhoferSpectrum)
/// and will keep it in memory to save loading time. If memory is a consern, then make sure that this object gets destructed when no longer needed.
/// </summary>
class FraunhoferSpectrumGeneration : public IFraunhoferSpectrumGenerator
{
public:
    /// <summary>
    /// Sets up the generation parameters
    /// </summary>
    /// <param name="highResolutionSolarAtlas">The full path to the high resolved solar atlas. This must be in nm air.</param>
    /// <param name="highResolutionCrossSections">The full path to a set of high resolved molecular cross section together with the total column for them.
    ///     These must have x-axis unit of nm air and y-axis unit of molecules / cm2</param>
    FraunhoferSpectrumGeneration(const std::string& highResolutionSolarAtlas, const std::vector<std::pair<std::string, double>>& highResolutionCrossSections)
        : solarAtlasFile(highResolutionSolarAtlas), crossSectionsToInclude(highResolutionCrossSections.size())
    {
        for (size_t ii = 0; ii < highResolutionCrossSections.size(); ++ii)
        {
            crossSectionsToInclude[ii].path = highResolutionCrossSections[ii].first;
            crossSectionsToInclude[ii].totalColumn = highResolutionCrossSections[ii].second;
        }
    }

    virtual std::unique_ptr<CSpectrum> GetFraunhoferSpectrum(
        const std::vector<double>& pixelToWavelengthMapping,
        const novac::CCrossSectionData& measuredInstrumentLineShape) override;

    virtual std::unique_ptr<CSpectrum> GetFraunhoferSpectrum(
        const std::vector<double>& pixelToWavelengthMapping,
        const novac::CCrossSectionData& measuredInstrumentLineShape,
        bool normalize) override;

    virtual std::unique_ptr<CSpectrum> GetDifferentialFraunhoferSpectrum(
        const std::vector<double>& wavelengthCalibration,
        const novac::CCrossSectionData& measuredInstrumentLineShape,
        double fwhmOfInstrumentLineShape) override;

#ifdef USE_DOAS_FIT

    /// <summary>
    /// Creates a Fraunhofer reference spectrum using the provided pixel-to-wavelength mapping and measured instrument line shape.
    ///     A DOAS fit is applied in order to determine the total columns of each high resolution cross section in this setup.
    ///     I.e., the total columns passed in to the constructor are ignored.
    /// This will throw an std::invalid_argument if no high resolution cross section was provided at creation (indicating an invalid setup)
    /// </summary>
    /// <param name="measuredInstrumentLineShape">A measurement of the instrument line shape</param>
    /// <param name="measuredSpectrum">A measured spectrum. </param>
    /// <returns>The high resolution solar spectrum convolved with the measured slf and resample to the provided grid.</returns>
    std::unique_ptr<CSpectrum> GetFraunhoferSpectrumMatching(
        const std::vector<double>& pixelToWavelengthMapping,
        const novac::CSpectrum& measuredSpectrum,
        const novac::CCrossSectionData& measuredInstrumentLineShape);

#endif // USE_DOAS_FIT

private:
    /// <summary>
    /// The path and filename of the solar atlas file to use.
    /// </summary>
    const std::string solarAtlasFile;

    struct AbsorbingCrossSection
    {
        AbsorbingCrossSection() = default;

        AbsorbingCrossSection(const std::pair<std::string, double>& value)
            : path(value.first), totalColumn(value.second)
        {
        }

        std::string path;
        double totalColumn = 0.0;
        std::unique_ptr<novac::CCrossSectionData> crossSectionData;
    };

    /// <summary>
    /// The path and total column of the high resolved absorption cross section files to include.
    /// </summary>
    std::vector<AbsorbingCrossSection> crossSectionsToInclude;

    /// <summary>
    /// The read in high resolution solar cross section, saved in order to reduce file-io time.
    /// </summary>
    std::unique_ptr<novac::CCrossSectionData> solarCrossSection;

    /// <summary>
    /// This creates a Fraunhofer spectrum using the provided Absorbing cross sections instead of using the member
    /// </summary>
    /// <param name="pixelToWavelengthMapping"></param>
    /// <param name="measuredInstrumentLineShape"></param>
    /// <param name="crossSectionsToInclude"></param>
    /// <returns></returns>
    std::unique_ptr<CSpectrum> GetFraunhoferSpectrum(
        const std::vector<double>& pixelToWavelengthMapping,
        const novac::CCrossSectionData& measuredInstrumentLineShape,
        std::vector<AbsorbingCrossSection>& crossSectionsToInclude,
        double fwhmOfInstrumentLineShape,
        bool normalize);

    std::unique_ptr<CSpectrum> GetDifferentialFraunhoferSpectrum(
        const std::vector<double>& pixelToWavelengthMapping,
        const novac::CCrossSectionData& measuredInstrumentLineShape,
        std::vector<AbsorbingCrossSection>& crossSectionsToInclude,
        double fwhmOfInstrumentLineShape);

    void ReadSolarCrossSection();
};

}

