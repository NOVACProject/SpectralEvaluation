#pragma once

#include <vector>
#include <string>

namespace MathFit
{
    class CVector;
}

namespace Evaluation{
    /** 
        The CReferenceData class holds information on the cross sections
        used in the fitting procedure. Each instance of this class holds
        the information of one reference used.
        The references can be either differential or not. */
    class CCrossSectionData
    {
    public:
        CCrossSectionData();
        
        CCrossSectionData &operator=(const CCrossSectionData &xs2);
        
        ~CCrossSectionData();
        
        // ----------------------- DATA -----------------------

        /** An array containing the wavelength information.*/
        std::vector <double> m_waveLength;

        /** An array containing the actual cross-section */
        std::vector <double> m_crossSection;

        // ----------------------- METHODS -----------------------

        /** Sets the cross-section information at the given pixel */
        void SetAt(int index, double wavel, double value);

        /** Sets the cross-section information to the values in the 
            supplied array */
        void Set(double *wavelength, double *crossSection, unsigned long pointNum);

        /** Sets the cross-section information to the values in the 
            supplied array */
        void Set(double *crossSection, unsigned long pointNum);

        /** Sets the cross-section information to the values in the 
            supplied array */
        void Set(MathFit::CVector &crossSection, unsigned long pointNum);

        /** Gets the cross section at the given pixel */
        double GetAt(unsigned int index) const;

        /** Gets the length of this cross section */
        unsigned long GetSize() const;

        /** Gets the wavelength at the given pixel */
        double GetWavelengthAt(unsigned int index) const;

        /** Reads the cross section from a file 
            @return 0 on success
            @return non-zero value on fail */
        int ReadCrossSectionFile(const std::string &fileName);
    };

    /** Performs a high-pass filtering of this cross section file */
    int HighPassFilter(CCrossSectionData& crossSection);
    int HighPassFilter_Ring(CCrossSectionData& crossSection);

    /** Multiplies this cross section with the given constant */
    int Multiply(CCrossSectionData& crossSection, double scalar);

    /** Calculates the logarithm of this cross section file */
    int Log(CCrossSectionData& crossSection);
}
