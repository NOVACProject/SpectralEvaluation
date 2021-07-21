#include <SpectralEvaluation/Spectra/SpectrumUtils.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/File/TXTFile.h>
#include "catch.hpp"

namespace novac
{
std::string GetMercuryFileName()
{
#ifdef _MSC_VER
    return std::string("../TestData/Hg_D2J2200_all.Master.Sample.txt");
#else
    return std::string("TestData/Hg_D2J2200_all.Master.Sample.txt");
#endif // _MSC_VER
}

// -------- Doing a calibration from a mercury spectrum --------
TEST_CASE("D2J2200 Mercury Spectrum", "[SpectrumUtils][InstrumentCalibration][IntegrationTest][MercurySpectrum]")
{
    CSpectrum spectrum;
    CTXTFile::ReadSpectrum(spectrum, GetMercuryFileName());

    SECTION("Find peak locates all mercury peaks")
    {
        const double minimumIntensity = 500.0;
        std::vector<SpectrumDataPoint> result;

        FindPeaks(spectrum, minimumIntensity, result);

        REQUIRE(9 == result.size());
    }

    // SECTION("Do mercury calibration locates all relevant mercury peaks and creates wavelength calibration")
    // {
    //     std::vector<SpectrumDataPoint> foundPeaks;
    //     std::vector<double> polynomial;
    //     MercuryCalibration(spectrum, 3, foundPeaks, polynomial);
    // 
    // }
}

}