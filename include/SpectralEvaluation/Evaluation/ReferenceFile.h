#pragma once

#include <memory>
#include <string>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>

namespace novac
{
    /* the options for the shift and squeeze */
    enum class SHIFT_TYPE{
        SHIFT_FREE,     /* Include the shift and squeeze in the fit */
        SHIFT_FIX,      /* Set the shift/squeee to a fixed value */
        SHIFT_LINK,     /* Links the shift/squeeze to the value of another reference */
        SHIFT_OPTIMAL,  /* (Only used in MobileDOAS) find an optimal shift/squeeze value from one spectrum in the traverse */
        SHIFT_LIMIT     /* (NovacProgram/NovacPPP) Set the shift/squeeze to free but limit it's maximum value. */
    };


    /** The class <b>CReferenceFile</b> is used to store
        the options for a single reference file that is
        to be included in the DOAS fit.
        Options stored are e.g. the options for the shift
        and squeeze and the path to the file containing the
        reference-spectrum. 
    */
    class CReferenceFile
    {
    public:
        CReferenceFile() = default;
        ~CReferenceFile() = default;

        /** Creates this reference file as an in-memory copy of the 
            provided cross section */
        CReferenceFile(const CCrossSectionData& contents);

        /** Copying this object */
        CReferenceFile& operator=(const CReferenceFile& other);
        CReferenceFile(const CReferenceFile& other);

        /** Moving this object */
        CReferenceFile& operator=(CReferenceFile&& other);
        CReferenceFile(CReferenceFile&& other);

        /** The name of the specie, e.g. SO2 or O3. */
        std::string m_specieName = "";

        /** The path to the reference file.
            This is a ready-to-use cross section on the correct wavelength grid convolved 
            with the instruments slit-function. This may or may not be high-pass filtered
            (that option is given by 'm_isFiltered' below).
            The main option is to use this file to evaluate the spectra. 
            This may be empty, in which case the cross section will be generated using
            m_crossSectionFile, m_slitFunctionFile and m_wavelengthCalibrationFile below. */
        std::string m_path = "";

        /** The path to a high-resolved cross section which can be convolved with 
            'm_slitFunctionFile' and resampled to the wavelength in 'm_wavelengthCalibrationFile'
            to generate a reference file. If m_path is empty then these three must be provided
            in order to generate the reference on the fly. */
        std::string m_crossSectionFile = "";
        std::string m_slitFunctionFile = "";
        std::string m_wavelengthCalibrationFile = "";

        /** The magic gas-factor is the conversion factor 
            between ppmm and mg and is necessary to calculate a flux.
            The factor 2.66 is for SO2, other gases needs other values. */
        double m_gasFactor = 2.66;

        /** The option for the column value. */
        SHIFT_TYPE m_columnOption = SHIFT_TYPE::SHIFT_FREE;

        /** The value for the column value (only used if m_columnOption is not SHIFT_FREE) */
        double m_columnValue = 0.0;

        /** if m_columnOption is SHIFT_LIMIT,
            this is the maximum column value allowed
            and m_columnValue is the minimum column value allowed */
        double m_columnMaxValue = 0.0;

        /** The option for the shift */
        SHIFT_TYPE m_shiftOption = SHIFT_TYPE::SHIFT_FIX;

        /** The value for the shift */
        double m_shiftValue = 0.0;

        /** if m_shiftOption is SHIFT_LIMIT,
            this is the maximum shift value allowed
            and m_shiftValue is the minimum shift value allowed */
        double m_shiftMaxValue = 0.0;

        /** The option for the squeeze */
        SHIFT_TYPE m_squeezeOption = SHIFT_TYPE::SHIFT_FIX;

        /** The value for the squeeze */
        double m_squeezeValue = 1.0;

        /** if m_squeezeOption is SHIFT_LIMIT,
            this is the maximum squeeze value allowed
            and m_squeezeValue is the minimum squeeze value allowed */
        double m_squeezeMaxValue = 1.0;

        /** This is true if this cross section file is high-pass filtered already on disk
            If this is false and the fit-type if HP_SUB or HP_DIV then the reference will be 
            High-pass filtered when they are being read in (FitWindow::ReadReferences). 
            NOTICE that in the NovacProgram are the references given for real-time evaluation always filtered
            and this flag must be set to true before running the evaluation. In NovacPPP are the references given
            always NOT filtered. */
        bool m_isFiltered = false;

        /** Set this to false to not include the reference into the DOAS fit. 
            This can be used to selectively include/exclude references. */
        bool m_include = true;

        /** The actual data of this reference - file. This is equal to
                nullptr if the reference data has not been read in yet, otherwise
                this will contain the data from the reference file. */
        std::unique_ptr<CCrossSectionData> m_data = nullptr;

        // ------------------------ METHODS ---------------------------

        /** Setting the column.
            if(SHIFT_TYPE) is SHIFT_LIMIT then 'value' is the lower limit
            and value2 is the upper limit 
            otherwise value2 is not used */
        void SetColumn(SHIFT_TYPE option, double value, double value2 = 1e16);

        /** Setting the shift
            if(SHIFT_TYPE) is SHIFT_LIMIT then 'value' is the lower limit
            and value2 is the upper limit 
            otherwise value2 is not used */
        void SetShift(SHIFT_TYPE option, double value, double value2 = 1e16);

        /** Setting the squeeze
            if(SHIFT_TYPE) is SHIFT_LIMIT then 'value' is the lower limit
            and value2 is the upper limit 
            otherwise value2 is not used */
        void SetSqueeze(SHIFT_TYPE option, double value, double value2 = 1e16);

        /** Reads the data of this reference file from disk.
            The file is taken from the member variable 'm_path' (which of course must be initialized first)
            and the result is written to 'm_data'. If this fails then m_data will be NULL.
            @return 0 on success. */
        int ReadCrossSectionDataFromFile();

        /** Performs a convolution using the files m_crossSectionFile, m_slitFunctionFile and m_wavelengthCalibrationFile
            and saves the result in m_data.
            @return 0 on success.*/
        int ConvolveReference();
    };
}