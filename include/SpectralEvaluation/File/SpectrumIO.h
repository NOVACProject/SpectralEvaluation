#pragma once

#include "MKPack.h"
#include <string>

namespace novac
{
class CSpectrum;
class CDateTime;

/** <b>CSpectrumIO</b> is a class for reading and storing Spectra */

class CSpectrumIO
{
public:
    CSpectrumIO();
    ~CSpectrumIO();

    static const int ERROR_NO_ERROR = 0;
    static const int ERROR_EOF = 1;
    static const int ERROR_COULD_NOT_OPEN_FILE = 2;
    static const int ERROR_CHECKSUM_MISMATCH = 3;
    static const int ERROR_SPECTRUM_TOO_LARGE = 4;  // the size of the uncompressed spectrum is too large to handle
    static const int ERROR_SPECTRUM_NOT_FOUND = 5;  // the given spectrum index was not found (end of file)
    static const int ERROR_DECOMPRESS = 6;  // an error occurred while decompressing the spectrum
    static const int ERROR_SPECTRUM_NOT_COMPLETE = 7;  // the whole spectrum was not saved
    static const int ERROR_COULD_NOT_CHANGE_POS = 8;  // the position in the file could not be changed

    /** Reads spectrum number 'spectrumNumber' in the provided spectrum file.
        Spectrum files in the MKPack format can contain an unlimited amount of spectra.
        @param fileName - The name of the spectrum file.
        @param spectrumNumber - The spectrum number in the file.
        @param spec - Will on successful return contain the desired spectrum.
        @return true if all is ok. */
    bool ReadSpectrum(const std::string& fileName, const int spectrumNumber, CSpectrum& spec, char* headerBuffer = NULL, int headerBufferSize = 0, int* headerSize = NULL);

    /** Reads the next spectrum in the provided spectrum file.
            The spectrum file (which must be in the .pak format) must be opened for reading
            in binary mode. File will not be closed by this routine.
        @param f - The opened spectrum file.
        @param spec - Will on successful return contain the desired spectrum.
        @return true if all is ok. */
    bool ReadNextSpectrum(FILE* f, CSpectrum& spec);

    /** Rewinds the gien file to the beginning and forwards the current position
            in the file to the beginning of spectrum number 'spectrumNumber' (zero-based index).
            Return true if all is ok, return false if the file is corrupt in some
                way or the spectrum number 'spectrumNumber' does not exist in this file. */
    bool FindSpectrumNumber(FILE* f, int spectrumNumber);

    /** Reads the next spectrum in the provided spectrum file.
            The spectrum file (which must be in the .pak format) must be opened for reading
            in binary mode. File will not be closed by this routine.
        @param f - The opened spectrum file.
        @param spec - Will on successful return contain the desired spectrum.
        @param headerBuffer - if this is not null it will on successfull return be filled
            with the full header of the spectrum in binary format. Useful if the header in the .pak
            file is of a newer version than the programs headerversion
        @param headerBufferSize - the size of the headerBuffer.
        @param headerSize - will on successfull return be the size of the binary header (in bytes)
        @return true if all is ok. */
    bool ReadNextSpectrum(FILE* f, CSpectrum& spec, int& headerSize, char* headerBuffer = NULL, int headerBufferSize = 0);

    /** Adds a new spectrum to the given pak-file. If the file does not exist, it will be created.
        @param fileName - The name of the spectrum file in which to store the spectrum.
        @param spec - The spectrum to save.
        @param headerBuffer - if not null then this will be written as header instead of the header in the CSpectrum.
            useful if the spectrum was read from a newer version of the .pak-file than this
            program can handle. Copying the header directly ensures that no data is lost.
        @param headersize - The size of the headerBuffer. only useful if headerBuffer is not null.
        */
    int AddSpectrumToFile(const std::string& fileName, const CSpectrum& spec, const char* headerBuffer = nullptr, int headerSize = 0, bool overwrite = false);

    /** Opens The spectrum file and counts how many spectra there are in the file.
            @param fileName - the name and path of the .pak-file to open
            @return - The number of spectra in the spectrum file. */
    int CountSpectra(const std::string& fileName);

    /** Opens the spectrum file and searches for the occurence of certain spectra inside.
        The function can e.g. try to localize the spectrum with the name 'zenith' inside the
        spectrum-file. 'specNamesToLookFor[i]' will then be 'zenith' and on return 'indices[i]'
        will be the (zero-based) index into the file which has this name i.e.
        'ReadSpectrum(fileName, indices[i], spec)' will return the 'zenith' spectrum.

        @param specNamesToLookFor - this function can optionally look for the occurrence
            of spectra with given names, e.g. 'offset, 'sky', 'dar' etc while scanning the .pak-file.
            The strings to look for are given in 'specnNamesToLookFor'

        @param numSpecNames - the length of the array 'specNamesToLookFor' (must be equal to
            the length of 'incides' as well)

        @param incides - if spectra with certain names are to be looked for then
            'indices' will on return the (zero-based) index

        @return - The number of spectra in the spectrum file */
    std::uint32_t ScanSpectrumFile(const std::string& fileName, const std::string* specNamesToLookFor, int numSpecNames, int* indices);

    /** If any error occurs in the reading of the file, this int is set to
        any of the errors defined above. */
    int m_lastError;

private:
    /** A buffer for reading data */
    unsigned char buffer[16384];

    /** A buffer for temporary output */
    long outbuf[16384];

    /** The maximum value read */
    double maxv;

    /** The last read header */
    struct MKZYhdr MKZY;

    /** ?? */
    long multisize;

    /** Reads a spectrum header from the supplied file. The result
            will be saved to the member-variable 'MKZY'. If a CSpectrum
            is provided, the header information will also be saved in the spectrum.
            @param headerBuffer - if this is not null it will on successfull return be filled
                    with the full header of the spectrum in binary format. Useful if the header in the .pak
                    file is of a newer version than the programs headerversion
            @param headerBufferSize - the size of the headerBuffer.
            @param headerSize - will on successfull return be the size of the binary header (in bytes)
            @return 0 - on success
            @return 1 - ...*/
    int ReadNextSpectrumHeader(FILE* f, int& headerSize, CSpectrum* spec = NULL, char* headerBuffer = NULL, int headerBufferSize = 0);

    /** Converts a time from std::uint32_t to CDateTime */
    void ParseTime(const std::uint32_t t, CDateTime& time) const;

    /** Converts a time from CDateTime to std::uint32_t */
    void WriteTime(std::uint32_t& t, const CDateTime& time) const;

    /** Converts a date from std::uint32_t to CDateTime */
    void ParseDate(const std::uint32_t d, CDateTime& day) const;

    /** Converts a date from CDateTime to std::uint32_t */
    void WriteDate(std::uint32_t& d, const CDateTime& day) const;
};
}