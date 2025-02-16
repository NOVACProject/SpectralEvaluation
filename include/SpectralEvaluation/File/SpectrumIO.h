#pragma once

#include <string>
#include <vector>
#include <SpectralEvaluation/NovacEnums.h>

namespace novac
{
class CSpectrum;
class CDateTime;
struct MKZYhdr;


/** CSpectrumIO is a class for reading and storing spectra
    using the compressed .pak file format used in the Novac project. */
class CSpectrumIO
{
public:

    // The maximum length of a read spectrum here, any longer spectra on disk will result in FileError::SpectrumToLarge
    static const int MaxOutputSpectrumLength = 4096;

    /** Reads spectrum number 'spectrumNumber' in the provided spectrum file.
        Spectrum files in the MKPack format can contain an unlimited amount of spectra.
        @param fileName - The name of the spectrum file.
        @param spectrumNumber - The spectrum number in the file.
        @param spec - Will on successful return contain the desired spectrum.
        @return true if all is ok. */
    bool ReadSpectrum(const std::string& fileName, int spectrumNumber, CSpectrum& spec, char* headerBuffer = nullptr, int headerBufferSize = 0, int* headerSize = nullptr);

    /** Reads the next spectrum in the provided spectrum file.
            The spectrum file (which must be in the .pak format) must be opened for reading
            in binary mode. File will not be closed by this routine.
        @param f - The opened spectrum file.
        @param spec - Will on successful return contain the desired spectrum.
        @return true if all is ok. */
    bool ReadNextSpectrum(FILE* f, CSpectrum& spec);

    /** Rewinds the given file to the beginning and forwards the current position
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
    bool ReadNextSpectrum(FILE* f, CSpectrum& spec, int& headerSize, char* headerBuffer = nullptr, int headerBufferSize = 0);

    /** Adds a new spectrum to the given pak-file. If the file does not exist, it will be created.
        @param fileName - The name of the spectrum file in which to store the spectrum.
        @param spec - The spectrum to save.
        @param headerBuffer - if not null then this will be written as header instead of the header in the CSpectrum.
            useful if the spectrum was read from a newer version of the .pak-file than this
            program can handle. Copying the header directly ensures that no data is lost.
        @param headersize - The size of the headerBuffer. only useful if headerBuffer is not null. */
    int AddSpectrumToFile(const std::string& fileName, const CSpectrum& spec, const char* headerBuffer = nullptr, int headerSize = 0, bool overwrite = false);

    /** Opens The spectrum file and counts how many spectra there are in the file.
        @param fileName - the name and path of the .pak-file to open
        @return - The number of spectra in the spectrum file.
        This sets m_lastError to 'ERROR_SPECTRUM_NOT_FOUND' upon successful completion.
        This sets m_lastError to 'ERROR_COULD_NOT_OPEN_FILE' if the file could not be opened.*/
    int CountSpectra(const std::string& fileName);

    /** Opens the spectrum file and searches for the occurence of certain spectra inside.
        The function can e.g. try to localize the spectrum with the name 'zenith' inside the
        spectrum-file. 'specNamesToLookFor[i]' will then be 'zenith' and on return 'indices[i]'
        will be the (zero-based) index into the file which has this name i.e.
        'ReadSpectrum(fileName, indices[i], spec)' will return the 'zenith' spectrum.
        @param specNamesToLookFor - this function can optionally look for the occurrence
            of spectra with given names, e.g. 'offset, 'sky', 'dar' etc while scanning the .pak-file.
            The strings to look for are given in 'specnNamesToLookFor'
        @param incides - if spectra with certain names are to be looked for then
            'indices' will on return the (zero-based) index
        @return - The number of spectra in the spectrum file.
        This sets m_lastError to 'ERROR_SPECTRUM_NOT_FOUND' upon successful completion.
        This sets m_lastError to 'ERROR_COULD_NOT_OPEN_FILE' if the file could not be opened.*/
    std::uint32_t ScanSpectrumFile(const std::string& fileName, const std::vector<std::string>& specNamesToLookFor, std::vector<int>& indices);

    /** If any error occurs in the reading of the file, this is set to any of the errors defined above. */
    FileError m_lastError = FileError::NoError;

private:

    enum class HeaderReadingStatusCode
    {
        Success,
        CouldNotReadData,
        InvalidHeader
    };

    /** Reads a spectrum header from the supplied file.
        @param MKZYHeader Will on successful return contain the read in result.
        @param spec If not null then the header information will also be saved in the spectrum (not the spectral data).
        @param headerBuffer - if this is not null it will on successfull return be filled
                with the full header of the spectrum in binary format. Useful if the header in the .pak
                file is of a newer version than the programs headerversion
        @param headerBufferSize - the size of the headerBuffer.
        @param headerSize - will on successfull return be the size of the binary header (in bytes) */
    HeaderReadingStatusCode ReadNextSpectrumHeader(FILE* f, struct MKZYhdr& MKZYHeader, int& headerSize, CSpectrum* spec = nullptr, char* headerBuffer = nullptr, int headerBufferSize = 0);

    // Moves the current position in the given file to the beginning of the next spectrum file.
    // This assumes that we have read the header and wants to skip the reading of the spectral data.
    // @param numberOfBytesToSkip is the (compressed) size of the spectrum to skip.
    bool GotoNextSpectrum(FILE* f, std::uint16_t numberOfBytesToSkip) const;
};
}