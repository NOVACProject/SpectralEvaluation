#pragma once

#include <string>

namespace Evaluation
{
    class CCrossSectionData;

    /* the options for the shift and squeeze */
    const enum SHIFT_TYPE{
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
        CReferenceFile();
        ~CReferenceFile();

        /** Copying this object */
        CReferenceFile &operator=(const CReferenceFile& other);
        CReferenceFile(const CReferenceFile& other);

        /** The name of the specie */
        std::string   m_specieName = "";

        /** The path to the reference file */
        std::string   m_path = "";

        /** The magic gas-factor is the conversion factor 
            between ppmm and mg */
        double        m_gasFactor = 2.66;

        /** The option for the column value. */
        SHIFT_TYPE    m_columnOption = SHIFT_FREE;

        /** The value for the column value (only used if m_columnOption is not SHIFT_FREE) */
        double        m_columnValue = 0.0;

        /** if m_columnOption is SHIFT_LIMIT,
            this is the maximum column value allowed
            and m_columnValue is the minimum column value allowed */
        double        m_columnMaxValue = 0.0;

        /** The option for the shift */
        SHIFT_TYPE    m_shiftOption = SHIFT_FIX;

        /** The value for the shift */
        double        m_shiftValue = 0.0;

        /** if m_shiftOption is SHIFT_LIMIT,
            this is the maximum shift value allowed
            and m_shiftValue is the minimum shift value allowed */
        double        m_shiftMaxValue = 0.0;

        /** The option for the squeeze */
        SHIFT_TYPE    m_squeezeOption = SHIFT_FIX;

        /** The value for the squeeze */
        double        m_squeezeValue = 1.0;

        /** if m_squeezeOption is SHIFT_LIMIT,
            this is the maximum squeeze value allowed
            and m_squeezeValue is the minimum squeeze value allowed */
        double        m_squeezeMaxValue = 1.0;

        /** This is true if this cross section file is filtered already on disk
            If this is false and the fit-tyype if HP_SUB or HP_DIV then the reference will be 
            High-pass filtered when they are being read in (FitWindow::ReadReferences). */
        bool          m_isFiltered = false;

        /** The actual data of this reference - file. This is equal to
                nullptr if the reference data has not been read in yet, otherwise
                this will contain the data from the reference file. 
                (only used in NovacPPP) */
        CCrossSectionData *m_data = nullptr;

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
            and the result is written to 'm_data'. If this fails then m_data will be NULL */
        int ReadCrossSectionDataFromFile();
    };
}