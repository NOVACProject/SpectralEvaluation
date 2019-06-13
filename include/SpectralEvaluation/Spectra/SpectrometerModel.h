#pragma once

#include <string>
#include <vector>

struct SpectrometerModel
{
    SpectrometerModel()
        : modelName("S2000"), maximumIntensity(4095), numberOfChannels(1)
    {
    }

    SpectrometerModel(const std::string& name, double maxIntensity, bool custom = true)
        :modelName(name), maximumIntensity(maxIntensity), numberOfChannels(1), isCustom(custom)
    {
    }

    std::string modelName = "S2000";
    double maximumIntensity = 4096;
    int numberOfChannels = 1;
    bool isCustom = true;

    bool IsUnknown() const { return modelName.size() == 0; }
};

class CSpectrometerDatabase
{
public:
    static CSpectrometerDatabase& GetInstance()
    {
        static CSpectrometerDatabase singletonInstance;
        return singletonInstance;
    }

    /** @return The properties of the provided spectrometer model.
        The model name comparison is not case sensitive.
        If the model cannot be found, then an 'unknown spectrometer' is returned. */
    SpectrometerModel GetModel(const std::string& modelName);

    /** @return The properties of the spectrometer with the provided index into this database.
        If the model cannot be found, then an 'unknown spectrometer' is returned. */
    SpectrometerModel GetModel(int modelIndex);

    /** Attempts to guess the spectrometer model from a provided serialnumber. */
    SpectrometerModel GuessModelFromSerial(const std::string& deviceSerialNumber);

    /** @return The index of the provided spectrometer model in the list returned by 'ListModels()'.
        The model name comparison is not case sensitive.
        @return -1 if the model cannot be found. */
    int GetModelIndex(const std::string& modelName);

    /** @return a list of all configured spectrometer models */
    std::vector<std::string> ListModels() const;

    /** @return true if a model with the given name exists in the database.
        The model name comparison is not case sensitive. */
    bool Exists(const std::string& modelName) const;

    /** Adds a new model to this database.
        @return true if successful */
    bool AddModel(const SpectrometerModel& newModel);

    static SpectrometerModel SpectrometerModel_Unknown() { return SpectrometerModel{ "", 4095, true }; }
    static SpectrometerModel SpectrometerModel_S2000() { return SpectrometerModel{ "S2000", 4095, true }; }
    static SpectrometerModel SpectrometerModel_USB2000() { return SpectrometerModel{ "USB2000", 4095, true }; }
    static SpectrometerModel SpectrometerModel_USB4000() { return SpectrometerModel{ "USB4000", 65535, true }; }
    static SpectrometerModel SpectrometerModel_HR2000() { return SpectrometerModel{ "HR2000", 4095, true }; }
    static SpectrometerModel SpectrometerModel_HR4000() { return SpectrometerModel{ "HR4000", 16535, true }; }
    static SpectrometerModel SpectrometerModel_QE65000() { return SpectrometerModel{ "QE65000", 65535, true }; }
    static SpectrometerModel SpectrometerModel_MAYAPRO() { return SpectrometerModel{ "MAYAPRO", 65535, true }; }

private:
    CSpectrometerDatabase();

    std::vector<SpectrometerModel> modelDb;
    const SpectrometerModel unknown = SpectrometerModel{ "", 4095 };
};
