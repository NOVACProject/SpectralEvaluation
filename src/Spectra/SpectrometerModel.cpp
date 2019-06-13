#include <SpectralEvaluation/Spectra/SpectrometerModel.h>
#include <SpectralEvaluation/Utils.h>

CSpectrometerDatabase::CSpectrometerDatabase()
{
    SpectrometerModel s2000{ "S2000", 4095 };
    modelDb.push_back(s2000);

    SpectrometerModel USB2000{ "USB2000", 4095 };
    modelDb.push_back(USB2000);

    SpectrometerModel USB4000{ "USB4000", 65535 };
    modelDb.push_back(USB4000);

    SpectrometerModel HR2000{ "HR2000", 4095 };
    modelDb.push_back(HR2000);

    SpectrometerModel HR4000{ "HR4000", 16535 };
    modelDb.push_back(HR4000);

    SpectrometerModel QE65000{ "QE65000", 65535 };
    modelDb.push_back(QE65000);

    SpectrometerModel MAYAPRO{ "MAYAPRO", 65535 };
    modelDb.push_back(MAYAPRO);
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
