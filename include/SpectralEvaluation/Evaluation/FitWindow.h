#pragma once

#include <SpectralEvaluation/Evaluation/ReferenceFile.h>
#include <SpectralEvaluation/Evaluation/DoasFitEnumDeclarations.h>
#include <vector>

// TODO: Remove
#define MAX_N_REFERENCES 10

namespace novac
{
    enum class RING_CALCULATION_OPTION
    {
        DO_NOT_CALCULATE_RING,  // the user supplies a ring spectrum, don't calculate any
        CALCULATE_RING,         // calculate a ring spectrum and include in the fit
        CALCULATE_RING_X2       // calculate 2 ring spectra and include in the fit (regular ring and ring * lambda^4)
    };

    /** CFitWindow is a class for describing the settings
        used for evaluation of a spectral region. Each CFitWindow object
        contains settings such as the lower and higher edge of the fit-region,
        which references to use, the length of the spectra, how the spectra were
        collected from the spectrometer (i.e. if it's a whole spectrum/partial spectrum)
        etc. */
    class CFitWindow
    {
    public:
        CFitWindow() = default;
        ~CFitWindow() = default;

        CFitWindow& operator=(const CFitWindow& other);
        CFitWindow(const CFitWindow& other);

        CFitWindow& operator=(CFitWindow&& other);
        CFitWindow(CFitWindow&& other);

        /** The lower edge of the fit window (in pixels) */
        int fitLow = 320;

        /** The upper edge of the fit window (in pixels) */
        int fitHigh = 460;

        /** The channel of the spectrometer for which this fit window is valid.
            Almost all spectrometers have only one channel however the Ocean Optics SD2000 may have
            up to eight channels inside one casing. Each such channel is a separate spectrometer, with different
            optical properties, however they all share the same serial number and are therefore handled as one device.
            For a normal spectrometer this is equal to 0
            For the SD2000-spectrometers this can be 0 or 1 */
        int channel = 0;

        /** The reference files to use */
        // TODO: Use std::vector for this, to remove unnecessary limit on number of references
        CReferenceFile ref[MAX_N_REFERENCES];

        /** The number of references to use */
        // TODO: Remove when 'ref' is a std::vector.
        int nRef = 0;

        /** The Fraunhofer-reference spectrum which we can use
            to determine the shift (&squeeze) between the measured
            spectra and the reference-files.
            It is necessary that the wavelength calibration of the 'fraunhoferRef'
            is same as the wavelength calibration of each reference file. */
        CReferenceFile fraunhoferRef;

        /** The order of the polynomial that is fitted in optical depth space. */
        int polyOrder = 5;

        /** Set to true to include a polynomial (currently only 0th degree)
            fitted in intensity space. Used to correct for stray light in spectrometer.
            Notice: This can only be included when fitType is FIT_POLY  */
        bool includeIntensitySpacePolyominal = false;

        /** The option for if a ring spectrum should be calculated and included in the fit.
            Possible options are
                1) don't calculate a ring
                2) calculate a ring from the sky spectrum and include this in the fit
                3) calculate a ring from the sky spectrum and include this and ring * lambda^4 into the fit. */
                // TODO: Resolve how this should be able to get a correct wavelength calibration for the sky spectrum !
        RING_CALCULATION_OPTION ringCalculation = RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING;

        /** The length of the spectra */
        int specLength = 2048;

        /** The pixel number that corresponds to the first data point in the spectrum
            (this is normally 0 but for the OceanOptics spectrometers it is possible
            to read out only a small part of the pixels, e.g. a spectrum that is
            the values in the pixels 321 to 470) */
        int startChannel = 0;

        /** The name of the fit window */
        std::string name = "SO2";

        /** The type of fit. Novac standard is FIT_HP_DIV. */
        FIT_TYPE fitType = FIT_TYPE::FIT_HP_DIV;

        /** Set to zero (0) if the sky-spectrum should have a fixed shift/squeeze of 0.0 / 1.0.
            Set to one (1) if the sky-spectrum should be allowed to shift freely.
            Set to two (2) if the sky-spectrum should have a fixed shift/squeeze set by the skyShift / skySqueeze parameters.
            only useful if fitType is FIT_HP_SUB or FIT_POLY */
        int shiftSky = 0;

        /** The amount of shift to apply to the sky-spectrum.
            Only used if shiftSky == 2 */
        double skyShift = 0.0;
        double skySqueeze = 1.0;

        /** Larger than 1 if the spectra are read out in an interlaced way.
            This parameter works in the same way as the CSpectrumInfo::m_interlaceStep */
        int interlaceStep = 1;

        /** 'UV' is true if the start wavelength of the spectrum is 290 nm or shorter */
        // TODO: replace this with the pixel-range which should be used for the offset-correction (which what this variable is used for)
        int UV = 1;

        /** True if the scan should be twice, once for finding the highest column value.
            The spectrum with the highest column value is then evluated again with
            the shift of all references set to free and the optimum shift value is found.
            The scan is then evaluated again using this shift value. */
        int findOptimalShift = 0;

        /** This is an optional set of child windows which are to be evaluated in conjunction with the current one.
            These are typically used for ratio-evaluations, where the spectrum is evaluated in one region
            for e.g. SO2 and another region for e.g. BrO, but could also in the future be used for
            e.g. studying what happens if the spectrum is evaluated in two different regions. */
        std::vector<CFitWindow> child;

        // ------------------ METHODS ----------------------------------

        /** Clearing the settings do default values */
        void Clear();
    };

    /** Reads all the references specified in the window from disk.
        This requires that all references have an absolute path such that the files can be found.
        @return true if all references could be read sucessfully. */
    bool ReadReferences(CFitWindow& window);

    /** Performs a high-pass filtering on all references which are not labelled as m_isFiltered. */
    void HighPassFilterReferences(CFitWindow& window);

    /** Performs a rescaling of all read in references to the unit of Molecules/cm2.
        It is here assumed that any reference which is NOT filtered is in the unit of Molecules/cm2
        and any reference which IS filtered is in the unit of PPMM. The reference will be scaled accordingly. */
    void ScaleReferencesToMolecCm2(CFitWindow& window);

    /** Adds the provided vector of data as a reference to the current window. No processing of the reference data will be done.
        The reference will be added as free column, fixed squeeze and shift either fixed (if linkShiftToIdx is -1) or as linked to the reference with the provided idx.
        This will increase the number of references included in the window by one. */
    void AddAsReference(CFitWindow& window, const std::vector<double>& referenceData, const std::string& name, int linkShiftToIdx = -1);

    /* Adds the provided vector of data as a 'sky' spectrum to the current window. No processing of the reference data will be done.
        The reference will be added as column fixed to -1, squeeze fixed to one and the shift according to the provided shift option.
        This will increase the number of references inlcuded in the window by one.
        @return the index of the newly inserted reference. */
    int AddAsSky(CFitWindow& window, const std::vector<double>& referenceData, SHIFT_TYPE shiftOption = SHIFT_TYPE::SHIFT_FIX);

}