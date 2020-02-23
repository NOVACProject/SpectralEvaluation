#pragma once

#include <string>
#include <vector>

class CSpectrumInfo;

/** The class <b>CSpectrometerModel</b> contains a collection of the
        the spectrometer types/models that can be connected with the scanning
        instruments. */
enum SPECTROMETER_MODEL {
    S2000,
    USB2000,
    USB4000,
    HR2000,
    HR4000,
    QE65000,
    MAYAPRO,
    AVASPEC,
    UNKNOWN_SPECTROMETER,
    NUM_CONF_SPEC_MODELS // the number of spectrometers that are configured
};

class CSpectrometerModel
{
public:
    CSpectrometerModel(void);
    ~CSpectrometerModel(void);

    /** Retrieves the maximum intensity for the supplied spectrometer model */
    static double GetMaxIntensity(const std::string& modelNumber);
    static double GetMaxIntensity(const SPECTROMETER_MODEL modelNumber);

    /** Converts a SPECTROMETER_MODEL to a string item */
    static bool ToString(SPECTROMETER_MODEL model, std::string &str);

    /** Converts a string item to a SPECTROMETER_MODEL */
    static SPECTROMETER_MODEL GetModel(const std::string &str);

    /** Gets the number of configured spectrometer models */
    static int GetNumSpectrometerModels();
};

struct SpectrometerModel
{
    SpectrometerModel(const std::string& name, SPECTROMETER_MODEL modelIdx, double maxIntensity, bool averages = false)
        :modelName(name), model(modelIdx), maximumIntensity(maxIntensity), averagesSpectra(averages)
    {
    }

    std::string modelName = "S2000";
    double maximumIntensity = 4096;
    int numberOfChannels = 1;
    bool averagesSpectra = false;
    SPECTROMETER_MODEL model = UNKNOWN_SPECTROMETER;
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
    If the model cannot be found, then an 'unknown spectrometer' is returned. */
    SpectrometerModel GetModel(const std::string& modelname);

    /** @return The properties of the spectrometer with the provided index into this database.
    If the model cannot be found, then an 'unknown spectrometer' is returned. */
    SpectrometerModel GetModel(int modelIndex);

    /** Attempts to guess the spectrometer model from a provided serialnumber. */
    SpectrometerModel GuessModelFromSerial(const std::string& deviceSerialNumber);

    /** @return The index of the provided spectrometer model in the list returned by 'ListModels()'.
    @return -1 if the model cannot be found. */
    int GetModelIndex(const std::string& modelname);

    /** @return a list of all configured spectrometer models */
    std::vector<std::string> ListModels() const;

    static SpectrometerModel SpectrometerModel_Unknown() { return SpectrometerModel{ "", UNKNOWN_SPECTROMETER, 4095 }; }
    static SpectrometerModel SpectrometerModel_S2000() { return SpectrometerModel{ "S2000", S2000, 4095 }; }
    static SpectrometerModel SpectrometerModel_USB2000() { return SpectrometerModel{ "USB2000", USB2000, 4095 }; }
    static SpectrometerModel SpectrometerModel_USB4000() { return SpectrometerModel{ "USB4000", USB4000, 65535 }; }
    static SpectrometerModel SpectrometerModel_HR2000() { return SpectrometerModel{ "HR2000", HR2000, 4095 }; }
    static SpectrometerModel SpectrometerModel_HR4000() { return SpectrometerModel{ "HR4000", HR4000, 16535 }; }
    static SpectrometerModel SpectrometerModel_QE65000() { return SpectrometerModel{ "QE65000", QE65000, 65535 }; }
    static SpectrometerModel SpectrometerModel_MAYAPRO() { return SpectrometerModel{ "MAYAPRO", MAYAPRO, 65535 }; }
    static SpectrometerModel SpectrometerModel_AVASPEC() { return SpectrometerModel{ "AVASPEC", AVASPEC, 262143, true }; }

private:
    CSpectrometerDatabase();

    std::vector<SpectrometerModel> modelDb;
    const SpectrometerModel unknown = SpectrometerModel{ "", UNKNOWN_SPECTROMETER, 4095 };
};

/// ----------- Free Functions, useful for spectrometer related questions ----------- 

// Calculates the maximum intensity necessary to display a spectrum with the 
//  provided collection properties, based on the number of spectra measured and
//  the spectrometer model.
double FullDynamicRangeForSpectrum(const CSpectrumInfo& info);
