#pragma once

namespace novac
{

// The list of instrument types available in the Novac Network
enum class NovacInstrumentType
{
    Gothenburg,
    Heidelberg
};

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

}  // namespace novac
