#include <SpectralEvaluation/NovacEnums.h>

namespace novac
{

std::string ToString(NovacInstrumentType type)
{
    switch (type)
    {
    case novac::NovacInstrumentType::Heidelberg:
        return "heidelberg";
    default:
        return "gothenburg";
    }
}

std::string ToString(FileError err)
{
    switch (err)
    {
    case FileError::NoError: return "No error";
    case FileError::EndOfFile: return "End of File";
    case FileError::CouldNotOpenfile: return "Could not open file";
    case FileError::ChecksumMismatch: return "Checksum mismatch";
    case FileError::SpectrumTooLarge: return "Spectrum too large to be opened";
    case FileError::SpectrumNotFound: return "Requested spectrum not found";
    case FileError::DecompressionError: return "Error decompressing spectrum";
    case FileError::SpectrumNotComplete: return "Spectrum not complete in file";
    case FileError::CouldNotChangeFilePosition: return "Could not change position in file";
    default:
        return "Unknown error";
    }

}

std::string ToString(MeasurementMode mode)
{
    switch (mode)
    {
    case novac::MeasurementMode::Unknown:
        return "unknown";
    case novac::MeasurementMode::Flux:
        return "flux";
    case novac::MeasurementMode::Windspeed:
        return "wind";
    case novac::MeasurementMode::Stratosphere:
        return "stratospheric";
    case novac::MeasurementMode::DirectSun:
        return "direct_sun";
    case novac::MeasurementMode::Composition:
        return "composition";
    case novac::MeasurementMode::Lunar:
        return "lunar";
    case novac::MeasurementMode::Troposphere:
        return "tropospheric";
    case novac::MeasurementMode::MaxDoas:
        return "max_doas";
    default:
        return "unknown";
    }
}

} // namespace novac
