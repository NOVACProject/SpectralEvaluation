#pragma once

#include <string>
#include <vector>

namespace novac
{

class CSpectrumInfo;

struct SpectrometerModel
{
    struct PixelRange
    {
        PixelRange();
        PixelRange(int low, int high);

        int from;
        int to;
    };

    SpectrometerModel()
        : modelName("S2000"), maximumIntensityForSingleReadout(4095), numberOfChannels(1)
    {
    }

    // Adds a not-custom spectrometer model (for internal use)
    SpectrometerModel(const std::string& name, double maxIntensity, int pixelsOnDetector = 2048, bool averages = false)
        :modelName(name), maximumIntensityForSingleReadout(maxIntensity), numberOfPixels(pixelsOnDetector), numberOfChannels(1), isCustom(false), averagesSpectra(averages)
    {
    }

    SpectrometerModel(const std::string& name, double maxIntensity, bool custom = true, bool averages = false)
        :modelName(name), maximumIntensityForSingleReadout(maxIntensity), numberOfChannels(1), isCustom(custom), averagesSpectra(averages)
    {
    }

    // The given name for this spectrometer model. Name comparisons are not case sensitive.
    std::string modelName = "S2000";

    // The maximum intensity of a single readout of this device, defines the dynamic range of the ADC.
    // NB: To get the maximum intensity of a particular spectrum call FullDynamicRangeForSpectrum
    double maximumIntensityForSingleReadout = 4096;

    // The number of pixels on the detector.
    int numberOfPixels = 2048;

    // The number of optically covered pixels, used to determine an electronic offset.
    PixelRange coveredPixels = PixelRange(20, 200);

    // The (maximum) number of channels.
    int numberOfChannels = 1;

    // Set to true for user defined models, used to separate built-in models from user defined ones.
    bool isCustom = true;

    // @return true if this model is not well defined.
    bool IsUnknown() const { return modelName.size() == 0 || modelName == "UNKNOWN"; }

    /** For models where averagesSpectra is true, the maximum intensity of any spectrum equals the maximumIntensityForSingleReadout.
        For models where averagesSpectra is false, the maximum intensity of any spectrum equals maximumIntensityForSingleReadout * number of spectra */
    bool averagesSpectra = false;

    /** Calculates the maximum possible intensity of the given spectrum,
        based on the number of spectra measured and the current spectrometer model. */
    double FullDynamicRangeForSpectrum(const CSpectrumInfo& info) const;
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
    static SpectrometerModel GuessModelFromSerial(const std::string& deviceSerialNumber);

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

    static SpectrometerModel SpectrometerModel_Unknown() { return SpectrometerModel{ "UNKNOWN", 4095, false, false }; }
    static SpectrometerModel SpectrometerModel_S2000() { return SpectrometerModel{ "S2000", 4095, false, false }; }
    static SpectrometerModel SpectrometerModel_USB2000() { return SpectrometerModel{ "USB2000", 4095, false, false }; }
    static SpectrometerModel SpectrometerModel_USB2000Plus() { return SpectrometerModel{ "USB2000+", 65535, false, false }; }
    static SpectrometerModel SpectrometerModel_USB4000() { return SpectrometerModel{ "USB4000", 65535, false, false }; }
    static SpectrometerModel SpectrometerModel_HR2000() { return SpectrometerModel{ "HR2000", 4095, false, false }; }
    static SpectrometerModel SpectrometerModel_HR4000() { return SpectrometerModel{ "HR4000", 16535, false, false }; }
    static SpectrometerModel SpectrometerModel_QE65000() { return SpectrometerModel{ "QE65000", 65535, 1024, false }; }
    static SpectrometerModel SpectrometerModel_MAYAPRO() { return SpectrometerModel{ "MAYAPRO", 65535, false, false }; }
    static SpectrometerModel SpectrometerModel_AVASPEC() { return SpectrometerModel{ "AVASPEC", 16383, false, false }; }
    static SpectrometerModel SpectrometerModel_FLAME() { return SpectrometerModel{ "FLAME", 65536, false, true }; }

private:
    CSpectrometerDatabase();

    std::vector<SpectrometerModel> modelDb;
    const SpectrometerModel unknown = SpectrometerModel{ "", 4095, false, false };
};

/// ----------- Free Functions, useful for spectrometer related questions ----------- 

/** Calculates the maximum intensity necessary to display a spectrum with the
    provided collection properties, based on the number of spectra measured and
    the spectrometer model. */
double FullDynamicRangeForSpectrum(const CSpectrumInfo& info);

}
