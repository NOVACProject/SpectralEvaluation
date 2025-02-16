#include "catch.hpp"
#include <SpectralEvaluation/Math/SpectrumMath.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>

// Simple Mock implementation of an IScanSpectrumSource
class MockSpectrumSource : public novac::IScanSpectrumSource
{
public:
    std::vector<novac::CSpectrum> m_spectra;
    int m_internalSpectrumCounter = 0;

    virtual int GetSpectrum(novac::LogContext context, int specNumber, novac::CSpectrum& spec) override
    {
        if (specNumber >= 0 && specNumber < (int)m_spectra.size())
        {
            spec = m_spectra[specNumber];
            return 0;
        }

        return 1;
    }

    virtual int GetSpectrumNumInFile() const override
    {
        return static_cast<int>(m_spectra.size());
    }

    virtual void ResetCounter() override
    {
        m_internalSpectrumCounter = 0;
    }

    virtual int GetNextMeasuredSpectrum(novac::LogContext context, novac::CSpectrum& spec) override
    {
        if ((size_t)m_internalSpectrumCounter < m_spectra.size())
        {
            spec = m_spectra[m_internalSpectrumCounter++];
            return 0;
        }
        return 1;
    }

    virtual int GetSky(novac::CSpectrum& /*spec*/) const override {
        return 1; // TODO: Implement
    };

    virtual int GetDark(novac::CSpectrum& /*result*/) const override {
        return 1; // TODO: Implement
    };

    virtual int GetOffset(novac::CSpectrum& /*spec*/) const override {
        return 1; // TODO: Implement
    };

    virtual int GetDarkCurrent(novac::CSpectrum& /*spec*/) const override {
        return 1; // TODO: Implement
    };

    virtual std::string GetFileName() const override {
        return "";
    };

    /** Retrieves the time when the first spectrum was collected */
    virtual novac::CDateTime GetScanStartTime() const override {
        return m_spectra.front().m_info.m_startTime;
    }

    /** Retrieves the time when the last spectrum was collected */
    virtual novac::CDateTime GetScanStopTime() const override {
        return m_spectra.back().m_info.m_startTime;
    }

    /** Retrieves the serial number of the device which collected this scan. */
    virtual std::string GetDeviceSerial() const override {
        return m_spectra.front().m_info.m_device;
    }
};

// TODO: Move
novac::CSpectrum GenerateSpectrum()
{
    novac::CSpectrum spectrum;
    spectrum.m_length = 123;
    for (int idx = 0; idx < 123; ++idx)
    {
        spectrum.m_data[idx] = (double)idx;
    }
    spectrum.m_wavelength.resize(spectrum.m_length);
    std::generate(begin(spectrum.m_wavelength), end(spectrum.m_wavelength), [n = 0]() mutable {
        return 300 + n++;
    });

    spectrum.m_info.m_numSpec = 15;
    spectrum.m_info.m_startTime.SetToNow();
    spectrum.m_info.m_startTime.Decrement(5); // 5 seconds
    spectrum.m_info.m_stopTime.SetToNow();

    return spectrum;
}

TEST_CASE("AverageSpectra", "[AverageSpectra]")
{
    SECTION("No indices provided returns zero.")
    {
        novac::CSpectrum result;
        MockSpectrumSource spectrumSource;

        int numberOfSpectraAveraged = novac::AverageSpectra(spectrumSource, std::vector<int>{}, result);

        REQUIRE(0 == numberOfSpectraAveraged);
        REQUIRE(0 == result.m_info.m_numSpec);
    }

    SECTION("One index provided - returns one and returns contents of set spectrum.")
    {
        novac::CSpectrum result;
        MockSpectrumSource spectrumSource;
        novac::CSpectrum originalSpectrum = GenerateSpectrum();
        spectrumSource.m_spectra.push_back(originalSpectrum);

        int numberOfSpectraAveraged = novac::AverageSpectra(spectrumSource, std::vector<int>{0}, result);

        REQUIRE(originalSpectrum.m_info.m_numSpec == numberOfSpectraAveraged);
        REQUIRE(originalSpectrum.m_info.m_numSpec == result.m_info.m_numSpec);
    }

    SECTION("Three indices provided - returns one and returns contents of set spectrum.")
    {
        novac::CSpectrum result;
        MockSpectrumSource spectrumSource;
        novac::CSpectrum originalSpectrum = GenerateSpectrum();
        spectrumSource.m_spectra.push_back(originalSpectrum);
        spectrumSource.m_spectra.push_back(originalSpectrum);
        spectrumSource.m_spectra.push_back(originalSpectrum);
        spectrumSource.m_spectra.push_back(originalSpectrum);
        spectrumSource.m_spectra.push_back(originalSpectrum);

        int numberOfSpectraAveraged = novac::AverageSpectra(spectrumSource, std::vector<int>{0, 2, 4}, result);

        REQUIRE(3 * originalSpectrum.m_info.m_numSpec == numberOfSpectraAveraged);
        REQUIRE(3 * originalSpectrum.m_info.m_numSpec == result.m_info.m_numSpec);
    }
}

