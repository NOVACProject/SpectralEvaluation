#include <SpectralEvaluation/File/SpectrumIO.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/File/MKPack.h>
#include "catch.hpp"
#include "TestData.h"

namespace novac
{

TEST_CASE("SpectrumIO CountSpectra", "[SpectrumIO][MKZY][IntegrationTests]")
{
    CSpectrumIO sut;

    SECTION("I2J8549")
    {
        // Act
        int number = sut.CountSpectra(TestData::GetMeasuredSpectrumName_I2J8549());

        // Assert
        REQUIRE(53 == number);
        REQUIRE(sut.m_lastError == FileError::SpectrumNotFound);
    }

    SECTION("I2J8549")
    {
        // Act
        int number = sut.CountSpectra(TestData::GetMeasuredSpectrumName_2002126M1());

        // Assert
        REQUIRE(52 == number); // one spectrum here could not be read
        REQUIRE(sut.m_lastError == FileError::SpectrumNotFound);
    }

    SECTION("Not existing file")
    {
        // Act
        int number = sut.CountSpectra("/some-not/existing/file.txt");

        // Assert
        REQUIRE(0 == number);
        REQUIRE(sut.m_lastError == FileError::CouldNotOpenfile);
    }
}

TEST_CASE("SpectrumIO FindSpectrumNumber", "[SpectrumIO][MKZY][IntegrationTests]")
{
    CSpectrumIO sut;
    CSpectrum readSpectrum;

    SECTION("Last spectrum in file")
    {
        // Arrange
        FILE* f = fopen(TestData::GetMeasuredSpectrumName_I2J8549().c_str(), "rb");
        REQUIRE(f != nullptr);

        // Act
        const bool success = sut.FindSpectrumNumber(f, 52);

        // Assert
        REQUIRE(success);

        // read the next spectrum, this should now be spectrum #52
        REQUIRE(sut.ReadNextSpectrum(f, readSpectrum));
        REQUIRE(readSpectrum.m_info.m_scanIndex == 52);
    }

    SECTION("One after last spectrum in file")
    {
        // Arrange
        FILE* f = fopen(TestData::GetMeasuredSpectrumName_I2J8549().c_str(), "rb");
        REQUIRE(f != nullptr);

        // Act
        const bool success = sut.FindSpectrumNumber(f, 53);

        // Assert
        REQUIRE(success == false);
    }

    SECTION("First spectrum after some spectra have been read")
    {
        // Arrange
        FILE* f = fopen(TestData::GetMeasuredSpectrumName_I2J8549().c_str(), "rb");
        REQUIRE(f != nullptr);

        const int numberOfSpectraToRead = 21;
        for (int k = 0; k < numberOfSpectraToRead; ++k)
        {
            REQUIRE(sut.ReadNextSpectrum(f, readSpectrum));
            REQUIRE(readSpectrum.m_info.m_scanIndex == k);
        }

        // Act
        const bool success = sut.FindSpectrumNumber(f, 3);

        // Assert
        REQUIRE(success);

        // read the next spectrum, this should now be spectrum #3
        REQUIRE(sut.ReadNextSpectrum(f, readSpectrum));
        REQUIRE(readSpectrum.m_info.m_scanIndex == 3);
    }
}

TEST_CASE("SpectrumIO ScanSpectrumFile", "[SpectrumIO][MKZY][IntegrationTests]")
{
    CSpectrumIO sut;
    const std::vector<std::string> strings = { std::string("sky"), std::string("zenith"), std::string("dark") };
    std::vector<int> indices(3, -1);

    SECTION("I2J8549")
    {
        // Act
        int number = sut.ScanSpectrumFile(TestData::GetMeasuredSpectrumName_I2J8549(), strings, indices);

        // Assert
        REQUIRE(53 == number);
        REQUIRE(sut.m_lastError == FileError::SpectrumNotFound);
        REQUIRE(indices[0] == 0); // first spectrum
        REQUIRE(indices[1] == -1); // not found
        REQUIRE(indices[2] == 1); // second spectrum
    }

    SECTION("I2J8549")
    {
        // Act
        int number = sut.ScanSpectrumFile(TestData::GetMeasuredSpectrumName_2002126M1(), strings, indices);

        // Assert
        REQUIRE(52 == number); // There are in fact 52 spectra in the file, but one cannot be read..
        REQUIRE(sut.m_lastError == FileError::SpectrumNotFound);
        REQUIRE(indices[0] == 0); // first spectrum
        REQUIRE(indices[1] == -1); // not found
        REQUIRE(indices[2] == 1); // second spectrum
    }

    SECTION("Not existing file")
    {
        // Act
        int number = sut.ScanSpectrumFile("/some-not/existing/file.txt", strings, indices);

        // Assert
        REQUIRE(0 == number);
        REQUIRE(sut.m_lastError == FileError::CouldNotOpenfile);
        REQUIRE(indices[0] == -1);
        REQUIRE(indices[1] == -1);
        REQUIRE(indices[2] == -1);
    }
}

TEST_CASE("MKZYhdr Alignment is correct", "[SpectrumIO][MKZY]")
{
    MKZYhdr sut;

    REQUIRE(sizeof(MKZYhdr) == 114);
    REQUIRE((size_t)&sut.flag - (size_t)&sut.ident == 51);
    REQUIRE((size_t)&sut.date - (size_t)&sut.ident == 52);
    REQUIRE((size_t)&sut.starttime - (size_t)&sut.ident == 56);
    REQUIRE((size_t)&sut.stoptime - (size_t)&sut.ident == 60);
    REQUIRE((size_t)&sut.lon - (size_t)&sut.ident == 72);
    REQUIRE((size_t)&sut.ADC - (size_t)&sut.ident == 98);
}

TEST_CASE("SpectrumIO ReadSpectrum - Not existing file", "[SpectrumIO][MKZY][IntegrationTests]")
{
    CSpectrumIO sut;
    CSpectrum ignoredSpectrum;
    MKZYhdr header;
    int headerSize = 0;
    const int spectrumNumberToRead = 0;

    // Act
    bool readSuccessfully = sut.ReadSpectrum(
        "/some-not/existing/file.txt",
        spectrumNumberToRead,
        ignoredSpectrum,
        (char*)&header,
        sizeof(header),
        &headerSize);

    // Assert
    REQUIRE(false == readSuccessfully);
    REQUIRE(sut.m_lastError == FileError::CouldNotOpenfile);
}

TEST_CASE("SpectrumIO ReadSpectrum can parse MKZY header of first spectrum in file", "[SpectrumIO][MKZY][IntegrationTests]")
{
    CSpectrumIO sut;
    CSpectrum ignoredSpectrum;
    MKZYhdr header;
    int headerSize = 0;
    const std::string filename = TestData::GetMeasuredSpectrumName_I2J8549();
    const int spectrumNumberToRead = 0;

    // Act
    bool readSuccessfully = sut.ReadSpectrum(
        filename,
        spectrumNumberToRead,
        ignoredSpectrum,
        (char*)&header,
        sizeof(header),
        &headerSize);

    // Assert
    REQUIRE(readSuccessfully);

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

TEST_CASE("SpectrumIO ReadSpectrum and file with checksum errors", "[SpectrumIO][MKZY][IntegrationTests]")
{
    CSpectrumIO sut;
    CSpectrum ignoredSpectrum;
    const std::string filename = TestData::GetMeasuredSpectrumName_2002126M1();

    MKZYhdr header;
    int headerSize = 0;

    SECTION("Spectrum before malformed can be read")
    {
        // Act
        bool readSuccessfully = sut.ReadSpectrum(
            filename,
            30,
            ignoredSpectrum,
            (char*)&header,
            sizeof(header),
            &headerSize);

        // Assert
        REQUIRE(readSuccessfully);
        REQUIRE(sut.m_lastError == FileError::NoError);
    }

    SECTION("Malformed spectrum cannot be read.")
    {
        // Act
        bool readSuccessfully = sut.ReadSpectrum(
            filename,
            31,
            ignoredSpectrum,
            (char*)&header,
            sizeof(header),
            &headerSize);

        // Assert
        REQUIRE(readSuccessfully == false);
        REQUIRE(sut.m_lastError == FileError::ChecksumMismatch);
    }

    SECTION("Spectrum after malformed can be read")
    {
        // Act
        bool readSuccessfully = sut.ReadSpectrum(
            filename,
            32,
            ignoredSpectrum,
            (char*)&header,
            sizeof(header),
            &headerSize);

        // Assert
        REQUIRE(readSuccessfully);
        REQUIRE(sut.m_lastError == FileError::NoError);
    }
}

TEST_CASE("SpectrumIO ReadSpectrum regular pak file.", "[SpectrumIO][MKZY][ReadSpectrum][IntegrationTests]")
{
    CSpectrumIO sut;
    CSpectrum spectrum;
    const std::string filename = TestData::GetMeasuredSpectrumName_I2J8549();

    SECTION("First spectrum in file")
    {
        // Act
        bool readSuccessfully = sut.ReadSpectrum(filename, 0, spectrum);

        // Assert
        REQUIRE(readSuccessfully);
        REQUIRE(sut.m_lastError == FileError::NoError);

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

    SECTION("Last spectrum in file")
    {
        // Act
        bool successfullyRead = sut.ReadSpectrum(filename, 52, spectrum);

        // Assert
        REQUIRE(successfullyRead);
        REQUIRE(sut.m_lastError == FileError::NoError);

        REQUIRE(spectrum.m_length == 2048);
        REQUIRE(spectrum.m_info.m_name == "scan");
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
        REQUIRE(spectrum.m_info.m_startTime.minute == 45);
        REQUIRE(spectrum.m_info.m_startTime.second == 6);

        REQUIRE(spectrum.m_info.m_stopTime.year == 2017);
        REQUIRE(spectrum.m_info.m_stopTime.month == 2);
        REQUIRE(spectrum.m_info.m_stopTime.day == 16);
        REQUIRE(spectrum.m_info.m_stopTime.hour == 12);
        REQUIRE(spectrum.m_info.m_stopTime.minute == 45);
        REQUIRE(spectrum.m_info.m_stopTime.second == 18);
        REQUIRE(spectrum.m_info.m_scanIndex == 52);
        REQUIRE(spectrum.m_info.m_compass == 34.0);
    }

    SECTION("One after last spectrum in file results in error")
    {
        // Act
        bool readSuccessfully = sut.ReadSpectrum(filename, 53, spectrum);

        // Assert
        REQUIRE(readSuccessfully == false);
        REQUIRE(sut.m_lastError == FileError::SpectrumNotFound);
    }
}


TEST_CASE("SpectrumIO ReadNextSpectrum regular pak file.", "[SpectrumIO][MKZY][ReadSpectrum][IntegrationTests]")
{
    CSpectrumIO sut;
    CSpectrum spectrum;
    const std::string filename = TestData::GetMeasuredSpectrumName_I2J8549();

    FILE* f = fopen(filename.c_str(), "rb");
    REQUIRE(f != nullptr);

    SECTION("First spectrum in file")
    {
        // Act
        bool readSuccessfully = sut.ReadNextSpectrum(f, spectrum);

        // Assert
        REQUIRE(readSuccessfully);
        REQUIRE(sut.m_lastError == FileError::NoError);

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

    SECTION("Reads all spectra in file")
    {
        const int numberOfSpectraInFile = 53;

        // Act
        for (int idx = 0; idx < numberOfSpectraInFile; ++idx)
        {
            bool readSuccessfully = sut.ReadNextSpectrum(f, spectrum);

            REQUIRE(readSuccessfully);
            REQUIRE(sut.m_lastError == FileError::NoError);
            REQUIRE(spectrum.m_length == 2048);
        }
    }

    SECTION("Returns error after last spectrum in file")
    {
        // Arrange
        const int numberOfSpectraInFile = 53;
        for (int idx = 0; idx < numberOfSpectraInFile; ++idx)
        {
            bool readSuccessfully = sut.ReadNextSpectrum(f, spectrum);

            REQUIRE(readSuccessfully);
        }

        // Act
        bool readSuccessfully = sut.ReadNextSpectrum(f, spectrum);

        // Assert
        REQUIRE(readSuccessfully == false);
        REQUIRE(sut.m_lastError == FileError::NoError);
    }
}
}