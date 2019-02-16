#pragma once

#include "ReferenceFile.h"

// TODO: Remove
#define MAX_N_REFERENCES 10

namespace Evaluation
{
    enum FIT_TYPE {
        FIT_HP_DIV,
        FIT_HP_SUB,
        FIT_POLY
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
        CFitWindow();

        ~CFitWindow();

        CFitWindow &operator=(const CFitWindow &w2);
        CFitWindow(const CFitWindow &wnd);

        /** The lower edge of the fit window (in pixels) */
        int fitLow = 320;

        /** The upper edge of the fit window (in pixels) */
        int fitHigh = 460;

        /** The channel which is used in this fit window
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

        /** The order of the polynomial that will also be fitted */
        int polyOrder = 5;

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
        FIT_TYPE fitType = FIT_HP_DIV;

        /** true if the sky-spectrum should be allowed to shift.
                only useful if fitType is FIT_HP_SUB or FIT_POLY */
        bool shiftSky;

        /** Larger than 1 if the spectra are read out in an interlaced way.
            This parameter works in the same way as the CSpectrumInfo::m_interlaceStep */
        int interlaceStep = 1;

        /** 'UV' is true if the start wavelength of the spectrum is 290 nm or shorter */
        // TODO: replace this with the pixel-range which should be used for the offset-correction (which this variable is used for)
        bool UV = true;

        /** True if the scan should be twice, once for finding the highest column value.
            The spectrum with the highest column value is then evluated again with
            the shift of all references set to free and the optimum shift value is found.
            The scan is then evaluated again using this shift value. */
        bool findOptimalShift = false;

        // ------------------ METHODS ----------------------------------

        /** Clearing the settings do default values */
        void Clear();
   };

   /** Reads all the references specified in the window from disk.
        @return true if all references could be read sucessfully. */
   bool ReadReferences(CFitWindow& window);
}