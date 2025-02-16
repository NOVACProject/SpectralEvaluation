#pragma once

#include <vector>
#include <string>
#include <SpectralEvaluation/Units.h>

namespace MathFit
{
class CVector;
}
namespace novac
{
class CSpectrum;

/**
    The CReferenceData class holds information on the cross sections
    used in the fitting procedure. Each instance of this class holds
    the information of one reference used.
    The references can be either differential or not. */
class CCrossSectionData
{
public:
    CCrossSectionData();

    /** Creates a full copy of another cross section */
    CCrossSectionData(const CCrossSectionData& other);

    /** Creates a copy of a range of another cross section */
    CCrossSectionData(const CCrossSectionData& other, double startWavelength, double endWavelength);

    /** Converts a CSpectrum to a CCrossSectionData */
    CCrossSectionData(const CSpectrum& spectrum);

    /** Converts an array of data into a CCrossSectionData. Only 'm_crossSection' will be filled in. */
    CCrossSectionData(const std::vector <double>& crossSectionData);

    CCrossSectionData& operator=(const CCrossSectionData& other);

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
    void Set(double* wavelength, double* crossSection, unsigned long pointNum);

    /** Sets the cross-section information to the values in the
        supplied array */
    void Set(double* crossSection, unsigned long pointNum);

    /** Sets the cross-section information to the values in the
        supplied array */
    void Set(MathFit::CVector& crossSection, unsigned long pointNum);

    /** Gets the cross section at the given pixel */
    double GetAt(unsigned int index) const;

    /** Gets the length of this cross section */
    unsigned long GetSize() const;

    /** Gets the wavelength at the given pixel */
    double GetWavelengthAt(unsigned int index) const;

    /** Locates and returns the index which corresponds to the provided wavelength.
        Returns -1 if no such index could be found. */
    double FindWavelength(double wavelength) const;

    /** Reads the cross section from a file
        @return 0 on success
        @return non-zero value on fail */
    int ReadCrossSectionFile(const std::string& fileName);

    /** Saves the contents of this cross section to file.
        @return 0 on success. */
    int SaveToFile(const std::string& filename) const;
};

// Performs a high-pass filtering of this cross section file.
// This will:
//  1) multiply by 2.5e15
//  2) exponentiate the cross section
//  3) high-pass filter (binomial, 500 passes)
//  4) log the cross section and
//  5) divide by 2.5e15 (optional)
// If is here assumed that the input cross section has a unit of cm2 / molecule.
// @param desiredOutputUnit if set to CrossSectionUnit::cm2_molecule then the reference is mutliplied by 2.5e15 
//  before the filtering and divide by the same number after. */
void HighPassFilter(CCrossSectionData& crossSection, CrossSectionUnit desiredOutputUnit);

/** Performs a high-pass filtering of a ring spectrum.
    This will:
        1) high-pass filter (binomial, 500 passes).
        2) log the cross section. */
void HighPassFilter_Ring(CCrossSectionData& crossSection);

/** Multiplies this cross section with the given constant */
int Multiply(CCrossSectionData& crossSection, double scalar);

/** Calculates the logarithm of this cross section file */
int Log(CCrossSectionData& crossSection);

/** Resamples a cross section to a given new resolution on a uniform grid.
    The result will be sampled on the same wavelength range (given by crossSection.m_wavelength)
        with a uniform step-size of 'resolution'.
    @param crossSection the cross section to resample.
    @param resolution the desired new resolution.
    @param result will be filled with the resampled cross section. */
void Resample(const CCrossSectionData& crossSection, double resolution, std::vector<double>& result);

/** Resamples a cross section to a given new resolution grid.
    The result will be sampled on the points given by resolution.
    If the newGrid stretches beyond the existing crossSection.m_wavelength then the result will be filled with zeros there
    @param crossSection the cross section to resample.
    @param newGrid the desired new resolution grid.
    @param result will be filled with the resampled cross section. */
void Resample(const CCrossSectionData& crossSection, const std::vector<double>& newGrid, std::vector<double>& result);

/** Shifts the provided data the given number of pixels (positive values corresponds to shifting towards
    higher indices). This will approximate the data using a cubic spline and then shift the spline. */
void Shift(std::vector<double>& data, double pixelCount);

/** Creates a simple gaussian function centered in the middle of the resulting CCrossSectionData
    @param gaussianSigma The sigma parameter of the Gaussian to generate
    @param deltaLambda The difference wavelength between two neighbouring points in the generated
        gaussian.m_waveLength array */
void CreateGaussian(double gaussianSigma, double deltaLambda, CCrossSectionData& gaussian);
}
