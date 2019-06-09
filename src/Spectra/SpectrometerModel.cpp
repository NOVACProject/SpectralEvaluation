#include <SpectralEvaluation/Spectra/SpectrometerModel.h>
#include <SpectralEvaluation/Utils.h>

CSpectrometerModel::CSpectrometerModel()
{
}

CSpectrometerModel::~CSpectrometerModel()
{
}

/** Retrieves the maximum intensity for the supplied spectrometer model */
double	CSpectrometerModel::GetMaxIntensity(const std::string& modelNumber)
{
    return GetMaxIntensity(GetModel(modelNumber));
}

double  CSpectrometerModel::GetMaxIntensity(const SPECTROMETER_MODEL modelNumber)
{
    switch (modelNumber) {
    case S2000:		return 4095;
    case USB2000: return 4095;
    case USB4000:	return 65535;
    case HR2000:	return 4095;
    case HR4000:	return 16535;
    case QE65000:	return 65535;
    case MAYAPRO:	return 65535;
    default:	return 4095;
    }
}

bool CSpectrometerModel::ToString(SPECTROMETER_MODEL model, std::string &str)
{
    if (S2000 == model) {
        str = "S2000";
        return true;
    }
    if (USB2000 == model) {
        str = "USB2000";
        return true;
    }
    if (USB4000 == model) {
        str = "USB4000";
        return true;
    }
    if (HR2000 == model) {
        str = "HR2000";
        return true;
    }
    if (HR4000 == model) {
        str = "HR4000";
        return true;
    }
    if (QE65000 == model) {
        str = "QE65000";
        return true;
    }
    if (MAYAPRO == model) {
        str = "MAYAPRO";
        return true;
    }

    str = "Unknown";
    return false;
}

SPECTROMETER_MODEL CSpectrometerModel::GetModel(const std::string &str)
{
    if (EqualsIgnoringCase(str, "S2000"))
        return S2000;
    if (EqualsIgnoringCase(str, "USB2000"))
        return USB2000;
    if (EqualsIgnoringCase(str, "HR2000"))
        return HR2000;
    if (EqualsIgnoringCase(str, "HR4000"))
        return HR4000;
    if (EqualsIgnoringCase(str, "USB4000"))
        return USB4000;
    if (EqualsIgnoringCase(str, "QE65000"))
        return QE65000;
    if (EqualsIgnoringCase(str, "MAYAPRO"))
        return MAYAPRO;

    // not defined
    return UNKNOWN_SPECTROMETER;
}

int	CSpectrometerModel::GetNumSpectrometerModels()
{
    return NUM_CONF_SPEC_MODELS;
}
