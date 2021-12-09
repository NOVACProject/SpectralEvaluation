#include <SpectralEvaluation/Spectra/SpectrometerModel.h>
#include <SpectralEvaluation/StringUtils.h>
#include <algorithm>
#include <SpectralEvaluation/Spectra/SpectrumInfo.h>

namespace novac
{
    SpectrometerModel::PixelRange::PixelRange()
        : from(0), to(0)
    {
    }

    SpectrometerModel::PixelRange::PixelRange(int low, int high)
        : from(std::min(low, high)), to(std::max(low, high))
    {
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
        else if (Contains(deviceSerialNumber, "MAYAPRO") || Contains(deviceSerialNumber, "MAYP"))
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

    bool CSpectrometerDatabase::Exists(const std::string& modelName) const
    {
        for (const SpectrometerModel& entry : modelDb)
        {
            if (EqualsIgnoringCase(modelName, entry.modelName))
            {
                return true;
            }
        }
        return false;
    }

    bool CSpectrometerDatabase::AddModel(const SpectrometerModel& newModel)
    {
        if (Exists(newModel.modelName))
        {
            return false;
        }

        modelDb.push_back(newModel);

        return true;
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

}
