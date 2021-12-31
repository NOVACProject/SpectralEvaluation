#include "catch.hpp"
#include <string.h>
#include <SpectralEvaluation/Spectra/SpectrumUtils.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/File/STDFile.h>
#include <sstream>
#include "TestData.h"

std::vector<double> CreateGaussian(double center, double sigma, const std::vector<double>& x); // located elsewhere

namespace novac
{
    static bool ContainsPointAtPixel(const std::vector<SpectrumDataPoint>& points, double expectedPixel, double margin)
    {
        std::stringstream peakList;
        for (const auto& pt : points)
        {
            if (std::abs(pt.pixel - expectedPixel) < margin)
            {
                return true;
            }
            peakList << " " << pt.pixel << " (with intensity: " << pt.intensity << ", base: [" << pt.leftPixel << ", " << pt.rightPixel << "])" << std::endl;
        }

        std::cout << "Failed to locate peak/valley at pixel: " << expectedPixel << " (+- " << margin << "). Found points are: " << std::endl;
        std::cout << peakList.str() << std::endl;

        return false;
    }

    // -------- Calculating the derivative --------
    TEST_CASE("Derivative", "[SpectrumUtils]")
    {
        CSpectrum inputSpectrum;
        std::vector<double> result;

        SECTION("Flat Line - returns all zeros")
        {
            inputSpectrum.m_length = 31;
            for (int ii = 0; ii < inputSpectrum.m_length; ++ii)
            {
                inputSpectrum.m_data[ii] = 1.0;
            }

            bool ret = novac::Derivative(inputSpectrum.m_data, inputSpectrum.m_length, result);

            REQUIRE(ret == true);
            REQUIRE(inputSpectrum.m_length == result.size());
            for (size_t ii = 0; ii < result.size(); ++ii)
            {
                REQUIRE(std::abs(result[ii] < 1e-9));
            }
        }
    }

    TEST_CASE("FindPeaks in synthetic spectra", "[SpectrumUtils][FindPeaks]")
    {
        CSpectrum inputSpectrum;
        inputSpectrum.m_length = 31;
        inputSpectrum.m_wavelength.resize(inputSpectrum.m_length);
        for (int ii = 0; ii < inputSpectrum.m_length; ++ii)
        {
            inputSpectrum.m_wavelength[ii] = 310.0 + ii * 0.3;
        }

        std::vector<novac::SpectrumDataPoint> result;

        SECTION("Flat Line - returns no peaks")
        {
            for (int ii = 0; ii < inputSpectrum.m_length; ++ii)
            {
                inputSpectrum.m_data[ii] = 1.0;
            }

            novac::FindPeaks(inputSpectrum, 0.0, result);

            REQUIRE(0 == result.size());
        }

        SECTION("Narrow ApproximateGaussian 1 - locates peak")
        {
            const double sigma = 0.2; // [nm]
            const double center = inputSpectrum.m_wavelength[inputSpectrum.m_length / 2];
            std::vector<double> spec = CreateGaussian(center, sigma, inputSpectrum.m_wavelength);
            memcpy(&inputSpectrum.m_data, spec.data(), spec.size() * sizeof(double));

            novac::FindPeaks(inputSpectrum, 0.0, result);

            REQUIRE(1 == result.size());
            REQUIRE(result[0].pixel == inputSpectrum.m_length / 2);
            REQUIRE(result[0].wavelength == center);
            REQUIRE(result[0].intensity == 1.0);
        }

        SECTION("Narrow ApproximateGaussian 2 - locates peak")
        {
            const double sigma = 0.2; // [nm]
            const double center = inputSpectrum.m_wavelength[0] + 0.42 * (inputSpectrum.m_wavelength[inputSpectrum.m_length - 1] - inputSpectrum.m_wavelength[0]);
            std::vector<double> spec = CreateGaussian(center, sigma, inputSpectrum.m_wavelength);
            memcpy(&inputSpectrum.m_data, spec.data(), spec.size() * sizeof(double));

            novac::FindPeaks(inputSpectrum, 0.0, result);

            REQUIRE(1 == result.size());
            REQUIRE(std::abs(result[0].wavelength - center) < 0.01);
        }
    }

    TEST_CASE("FindPeaks in measured spectra spectra from MAYP11440", "[SpectrumUtils][FindPeaks][MayaPro]")
    {
        // Prepare by reading in the spectrum
        CSpectrum inputSpectrum;
        CSTDFile::ReadSpectrum(inputSpectrum, TestData::GetSkySpectrumName_MAYP11440());
        {
            CSpectrum darkSpectrum;
            CSTDFile::ReadSpectrum(darkSpectrum, TestData::GetDarkSpectrumName_MAYP11440());
            inputSpectrum.Sub(darkSpectrum);
        }

        std::vector<novac::SpectrumDataPoint> result;

        SECTION("Very high threshold - does not return any peaks")
        {
            const double threshold = 1e5; // above any points in the spectrum.
            novac::FindPeaks(inputSpectrum, threshold, result);

            REQUIRE(0 == result.size());
        }

        SECTION("Returns found peaks sorted with increasing pixel.")
        {
            const double goodThreshold = 1000;
            novac::FindPeaks(inputSpectrum, goodThreshold, result);

            REQUIRE(result.size() > 10);
            for (size_t idx = 1; idx < result.size(); ++idx)
            {
                REQUIRE(result[idx].pixel > result[idx - 1].pixel);
            }
        }

        SECTION("Threshold 1000 - Returns only peaks with intensity above threshold.")
        {
            const double goodThreshold = 1000;
            novac::FindPeaks(inputSpectrum, goodThreshold, result);

            REQUIRE(result.size() > 0);
            for (const auto& peak : result)
            {
                REQUIRE(peak.intensity >= goodThreshold);
            }
        }

        SECTION("Threshold 30000 - Returns only peaks with intensity above threshold.")
        {
            const double goodThreshold = 30000;
            novac::FindPeaks(inputSpectrum, goodThreshold, result);

            REQUIRE(result.size() > 0);
            for (const auto& peak : result)
            {
                REQUIRE(peak.intensity >= goodThreshold);
            }
        }

        SECTION("Finds expected peaks in measured spectrum")
        {
            const double goodThreshold = 1000;
            novac::FindPeaks(inputSpectrum, goodThreshold, result);

            REQUIRE(ContainsPointAtPixel(result, 1719.0, 2.0));
            REQUIRE(ContainsPointAtPixel(result, 1691.0, 2.0));
            REQUIRE(ContainsPointAtPixel(result, 1748.2, 1.0));
            REQUIRE(ContainsPointAtPixel(result, 1849.0, 1.0));
            REQUIRE(ContainsPointAtPixel(result, 1922.0, 1.0));
            REQUIRE(ContainsPointAtPixel(result, 1984.0, 1.5));
            REQUIRE(ContainsPointAtPixel(result, 2017.0, 2.0));
            REQUIRE(ContainsPointAtPixel(result, 666.7, 2.0));
            REQUIRE(ContainsPointAtPixel(result, 686.8, 2.0));
        }
    }

    TEST_CASE("FindValleys in measured spectra spectra from MAYP11440", "[SpectrumUtils][FindValleys][MayaPro]")
    {
        // Prepare by reading in the spectrum
        CSpectrum inputSpectrum;
        CSTDFile::ReadSpectrum(inputSpectrum, TestData::GetSkySpectrumName_MAYP11440());
        {
            CSpectrum darkSpectrum;
            CSTDFile::ReadSpectrum(darkSpectrum, TestData::GetDarkSpectrumName_MAYP11440());
            inputSpectrum.Sub(darkSpectrum);
        }

        std::vector<novac::SpectrumDataPoint> result;

        SECTION("Very high threshold - does not return any valleys")
        {
            const double threshold = 1e5; // above any points in the spectrum.
            novac::FindValleys(inputSpectrum, threshold, result);

            REQUIRE(0 == result.size());
        }

        SECTION("Returns found valleys sorted with increasing pixel.")
        {
            const double goodThreshold = 1000;
            novac::FindValleys(inputSpectrum, goodThreshold, result);

            REQUIRE(result.size() > 10);
            for (size_t idx = 1; idx < result.size(); ++idx)
            {
                REQUIRE(result[idx].pixel > result[idx - 1].pixel);
            }
        }

        SECTION("No overlap between found valleys.")
        {
            const double goodThreshold = 1000;
            novac::FindValleys(inputSpectrum, goodThreshold, result);

            REQUIRE(result.size() > 10);
            for (size_t idx = 1; idx < result.size(); ++idx)
            {
                REQUIRE(result[idx].leftPixel > result[idx - 1].rightPixel);
            }
        }

        SECTION("Threshold 1000 - Returns only valleys with intensity above threshold.")
        {
            const double goodThreshold = 1000;
            novac::FindValleys(inputSpectrum, goodThreshold, result);

            REQUIRE(result.size() > 0);
            for (const auto& valley : result)
            {
                REQUIRE(valley.intensity >= goodThreshold);
            }
        }

        SECTION("Threshold 30000 - Returns only valleys with intensity above threshold.")
        {
            const double goodThreshold = 30000;
            novac::FindValleys(inputSpectrum, goodThreshold, result);

            REQUIRE(result.size() > 0);
            for (const auto& valley : result)
            {
                REQUIRE(valley.intensity >= goodThreshold);
            }
        }

        SECTION("Finds expected valleys in measured spectrum")
        {
            const double goodThreshold = 1000;
            novac::FindValleys(inputSpectrum, goodThreshold, result);

            REQUIRE(ContainsPointAtPixel(result, 1972.0, 2.0));
            REQUIRE(ContainsPointAtPixel(result, 2002.0, 2.0));
            REQUIRE(ContainsPointAtPixel(result, 1905.0, 1.0));
            REQUIRE(ContainsPointAtPixel(result, 1765.0, 1.0));

            REQUIRE(ContainsPointAtPixel(result, 1677.0, 1.0));
            REQUIRE(ContainsPointAtPixel(result, 1700.5, 1.0));
            REQUIRE(ContainsPointAtPixel(result, 1730.5, 1.0));

            REQUIRE(ContainsPointAtPixel(result, 1303.0, 2.0));
            REQUIRE(ContainsPointAtPixel(result, 705.0, 2.0));
            REQUIRE(ContainsPointAtPixel(result, 675.0, 2.0));
            REQUIRE(ContainsPointAtPixel(result, 569.0, 2.0));
        }
    }
}
