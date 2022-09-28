#include <SpectralEvaluation/File/SpectrumIO.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/File/MKPack.h>
#include "catch.hpp"
#include "TestData.h"

namespace novac
{
    TEST_CASE("SpectrumIO CountSpectra", "[SpectrumIO][IntegrationTests]")
    {
        CSpectrumIO sut;
        int number = sut.CountSpectra(TestData::GetMeasuredSpectrumName_I2J8549());

        REQUIRE(53 == number);
    }

    TEST_CASE("MKZYhdr Alignment is correct", "[MKZYhdr]")
    {
        MKZYhdr sut;

        REQUIRE((size_t)&sut.flag - (size_t)&sut.ident == 51);
        REQUIRE((size_t)&sut.date - (size_t)&sut.ident == 52);
        REQUIRE((size_t)&sut.starttime - (size_t)&sut.ident == 56);
        REQUIRE((size_t)&sut.stoptime - (size_t)&sut.ident == 60);
        REQUIRE((size_t)&sut.lon - (size_t)&sut.ident == 72);
        REQUIRE((size_t)&sut.ADC - (size_t)&sut.ident == 98);

        // TODO: Add more tests here, to make sure that the values are correctly aligned further on...
    }


    TEST_CASE("SpectrumIO ReadSpectrum can parse MKZY header of first spectrum from file", "[SpectrumIO][MKZY][IntegrationTests]")
    {
        CSpectrumIO sut;
        CSpectrum ignoredSpectrum;
        MKZYhdr header;
        int headerSize = 0;
        int returnValue = sut.ReadSpectrum(TestData::GetMeasuredSpectrumName_I2J8549(), 0, ignoredSpectrum, (char*)&header, sizeof(header), &headerSize);

        REQUIRE(returnValue);

        REQUIRE(header.ident[0] == 'M');
        REQUIRE(header.ident[1] == 'K');
        REQUIRE(header.ident[2] == 'Z');
        REQUIRE(header.ident[3] == 'Y');

        REQUIRE(header.hdrsize == 114);
        REQUIRE(header.pixels == 2048);
        REQUIRE(header.viewangle == 0);
        REQUIRE(header.scans == 15);
        REQUIRE(header.exptime == -704);
        REQUIRE(header.channel == 0);
        REQUIRE(header.flag == 0);
    }

    TEST_CASE("SpectrumIO ReadSpectrum can read first spectrum from file", "[SpectrumIO][ReadSpectrum][IntegrationTests]")
    {
        CSpectrumIO sut;
        CSpectrum spectrum;
        int returnValue = sut.ReadSpectrum(TestData::GetMeasuredSpectrumName_I2J8549(), 0, spectrum);

        REQUIRE(returnValue);

        REQUIRE(spectrum.m_length == 2048);
        REQUIRE(spectrum.m_info.m_name == "sky");
        REQUIRE(spectrum.m_info.m_device == "I2J8549");

        REQUIRE(spectrum.m_info.m_numSpec == 15);
        REQUIRE(spectrum.m_info.m_exposureTime == 704);
        REQUIRE(spectrum.m_info.m_channel == 0);
        REQUIRE(spectrum.m_info.m_flag == 0);
        REQUIRE(spectrum.m_info.m_startChannel == 0);

        REQUIRE(spectrum.m_info.m_startTime.year == 2017);
        REQUIRE(spectrum.m_info.m_startTime.month == 2);
        REQUIRE(spectrum.m_info.m_startTime.day == 16);
        REQUIRE(spectrum.m_info.m_startTime.hour == 12);
        REQUIRE(spectrum.m_info.m_startTime.minute == 30);
        REQUIRE(spectrum.m_info.m_startTime.second == 44);

        REQUIRE(spectrum.m_info.m_stopTime.year == 2017);
        REQUIRE(spectrum.m_info.m_stopTime.month == 2);
        REQUIRE(spectrum.m_info.m_stopTime.day == 16);
        REQUIRE(spectrum.m_info.m_stopTime.hour == 12);
        REQUIRE(spectrum.m_info.m_stopTime.minute == 30);
        REQUIRE(spectrum.m_info.m_stopTime.second == 56);

        REQUIRE(spectrum.m_info.m_scanIndex == 0);

        REQUIRE(spectrum.m_info.m_compass == 34.0);

    }
}