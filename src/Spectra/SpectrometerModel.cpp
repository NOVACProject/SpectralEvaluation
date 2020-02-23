#include <SpectralEvaluation/Spectra/SpectrometerModel.h>
#include <SpectralEvaluation/Utils.h>
#include <SpectralEvaluation/Spectra/SpectrumInfo.h>

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
    case AVASPEC:   return 65535;
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
    if (AVASPEC == model) {
        str = "AVASPEC";
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
    if (EqualsIgnoringCase(str, "AVASPEC"))
        return AVASPEC;

    // not defined
    return UNKNOWN_SPECTROMETER;
}

int	CSpectrometerModel::GetNumSpectrometerModels()
{
    return NUM_CONF_SPEC_MODELS;
}

CSpectrometerDatabase::CSpectrometerDatabase()
{
    SpectrometerModel s2000 = SpectrometerModel_S2000();
    modelDb.push_back(s2000);

    SpectrometerModel USB2000 = SpectrometerModel_USB2000();
    modelDb.push_back(USB2000);

    SpectrometerModel USB4000 = SpectrometerModel_USB4000();
    modelDb.push_back(USB4000);

    SpectrometerModel HR2000 = SpectrometerModel_HR2000();
    modelDb.push_back(HR2000);

    SpectrometerModel HR4000 = SpectrometerModel_HR4000();
    modelDb.push_back(HR4000);

    SpectrometerModel QE65000 = SpectrometerModel_QE65000();
    modelDb.push_back(QE65000);

    SpectrometerModel MAYAPRO = SpectrometerModel_MAYAPRO();
    modelDb.push_back(MAYAPRO);

    SpectrometerModel AVASPEC = SpectrometerModel_AVASPEC();
    modelDb.push_back(AVASPEC);
}

SpectrometerModel CSpectrometerDatabase::GuessModelFromSerial(const std::string& deviceSerialNumber)
{
    if (Contains(deviceSerialNumber, "D2J"))
    {
        return CSpectrometerDatabase::SpectrometerModel_S2000();
    }
    else if (Contains(deviceSerialNumber, "I2J"))
    {
        return CSpectrometerDatabase::SpectrometerModel_S2000();
    }
    else if (Contains(deviceSerialNumber, "USB2"))
    {
        return CSpectrometerDatabase::SpectrometerModel_USB2000();
    }
    else if (Contains(deviceSerialNumber, "USB4C"))
    {
        return CSpectrometerDatabase::SpectrometerModel_USB4000();
    }
    else if (Contains(deviceSerialNumber, "HR2"))
    {
        return CSpectrometerDatabase::SpectrometerModel_HR2000();
    }
    else if (Contains(deviceSerialNumber, "HR4"))
    {
        return CSpectrometerDatabase::SpectrometerModel_HR4000();
    }
    else if (Contains(deviceSerialNumber, "QE"))
    {
        return CSpectrometerDatabase::SpectrometerModel_QE65000();
    }
    else if (Contains(deviceSerialNumber, "MAYAPRO"))
    {
        return CSpectrometerDatabase::SpectrometerModel_MAYAPRO();
    }
    else if (Contains(deviceSerialNumber, "M1") && deviceSerialNumber.size() == 9)
    {
        return CSpectrometerDatabase::SpectrometerModel_AVASPEC();
    }
    else
    {
        return CSpectrometerDatabase::SpectrometerModel_Unknown();
    }
}

SpectrometerModel CSpectrometerDatabase::GetModel(const std::string& modelname)
{
    for (const SpectrometerModel& entry : modelDb)
    {
        if (EqualsIgnoringCase(modelname, entry.modelName))
        {
            return entry;
        }
    }

    return unknown;
}

SpectrometerModel CSpectrometerDatabase::GetModel(int modelIndex)
{
    if ((unsigned int)modelIndex >= modelDb.size())
    {
        return unknown;
    }
    return modelDb[modelIndex];
}

int CSpectrometerDatabase::GetModelIndex(const std::string& modelname)
{
    for (size_t ii = 0; ii < modelDb.size(); ++ii)
    {
        if (EqualsIgnoringCase(modelname, modelDb[ii].modelName))
        {
            return (int)ii;
        }
    }

    return -1;
}

std::vector<std::string> CSpectrometerDatabase::ListModels() const
{
    std::vector<std::string> items;
    items.reserve(modelDb.size());

    for (const SpectrometerModel& entry : modelDb)
    {
        items.push_back(entry.modelName);
    }

    return items;
}

double FullDynamicRangeForSpectrum(const CSpectrumInfo& info)
{
    SpectrometerModel model = CSpectrometerDatabase::GetInstance().GuessModelFromSerial(info.m_device);

    if (model.averagesSpectra)
    {
        return model.maximumIntensity;
    }
    else if (info.m_numSpec > 0)
    {
        return model.maximumIntensity * info.m_numSpec;
    }
    else
    {
        // the number of co-added spectra can sometimes be wrongly set to zero. Handle this case as well...
        return (long)(model.maximumIntensity * (info.m_peakIntensity / model.maximumIntensity));
    }
}