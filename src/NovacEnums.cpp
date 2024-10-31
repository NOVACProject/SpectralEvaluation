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
