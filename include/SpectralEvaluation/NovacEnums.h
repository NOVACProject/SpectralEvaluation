#pragma once

#include <string>

namespace novac
{

// The list of instrument types available in the Novac Network
enum class NovacInstrumentType
{
    Gothenburg,
    Heidelberg
};

// Returns a standard name for the given instrument type.
std::string ToString(NovacInstrumentType type);

// Various types of errors which could occur when reading a spectrum file.
enum class FileError
{
    NoError,
    EndOfFile,
    CouldNotOpenfile,
    ChecksumMismatch,
    SpectrumTooLarge,   // the size of the uncompressed spectrum is too large to handle
    SpectrumNotFound,   // the given spectrum index was not found (end of file)
    DecompressionError, // an error occurred while decompressing the spectrum
    SpectrumNotComplete,        // the whole spectrum was not saved
    CouldNotChangeFilePosition, // the position in the file could not be changed
};

// Returns the given error code nicely formatted
std::string ToString(FileError err);

// The various kinds of measurement modes that we have.
// This determines the intent behind the measurement.
enum class MeasurementMode
{
    Unknown,
    Flux,
    Windspeed,
    Stratosphere,
    DirectSun,
    Composition,
    Lunar,
    Troposphere,
    MaxDoas
};

// Returns a standard name for the given measurement mode.
std::string ToString(MeasurementMode mode);

}  // namespace novac
