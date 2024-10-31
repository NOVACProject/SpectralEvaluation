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
