#pragma once

#include <vector>
#include <string>
#include <SpectralEvaluation/Definitions.h>

namespace novac
{
class BasicScanEvaluationResult;

/** The class CPlumeInScanProperty is used to describe
    how a plume is seen by a scan. This incorporates properties
    such as the angle at which the centre of the plume is seen
    or an estimate of how large portion of the plume is seen. */
class CPlumeInScanProperty
{
public:
    CPlumeInScanProperty() = default;
    ~CPlumeInScanProperty() = default;

    /** Assignment operator */
    CPlumeInScanProperty& operator=(const CPlumeInScanProperty& p2) = default;

    /** The centre position of the plume for the first motor (only in standard NOVAC instrument). Scan angle in degrees.*/
    double plumeCenter = NOT_A_NUMBER;

    /** The centre position of the plume, for the second motor (only available in Heidelberg type instruments). Scan angle in degrees. */
    double plumeCenter2 = NOT_A_NUMBER;

    /** The estimated error in the estimation of the plume centre position for the first motor. */
    double plumeCenterError = NOT_A_NUMBER;

    /** The estimated error in the estimation of the plume centre position for the second motor (if any). */
    double plumeCenterError2 = NOT_A_NUMBER;

    /** The offset of the scan, describes the column
        of gas in the sky-spectrum */
    double offset = 0.0;

    /** The completeness of the scan. This is a value
        ranging from 0.5 to 1.0 trying to estimate
        how large portion of the plume is visible
        in the scan. */
    double completeness = 0.0;

    /** The edges of the plume. These are the scan angles (first motor only)
        where the colmn of the plume has dropped to 1/e from it's highest value. */
    double plumeEdgeLow = NOT_A_NUMBER;
    double plumeEdgeHigh = NOT_A_NUMBER;

    /** The locations where the column values of the plume has dropped to half from it's highest value */
    double plumeHalfHigh = NOT_A_NUMBER;
    double plumeHalfLow = NOT_A_NUMBER;

};

// --------------------------------------------------------------------
// -------------- CALCULATING IF WE SEE THE PLUME ---------------------
// --------------------------------------------------------------------

/** Finds the plume in the supplied scan. Return value is true if there is a plume, otherwise false
    @param scanAngles - the scanAngles for the measurements.
    @param phi - the second scan angle for the measurements, only for Heidelberg-instruments with two motors.
    @param columns - the slant columns for the measurements. Must be from a normal scan
    @param columnErrors - the corresponding slant column errors for the measurement.
    @param badEvaluation - the result of the quality judgement of the measured columns,
                        badEvaluation[i] = true means a bad value
    @param numPoints - the number of points in the scan. Must also be the length
                        of the vectors 'columns', 'columnErrors', and 'badEvaluation'
    @param plumeOffset - the offset of the plume. Must have been calculated by a prior call to CalculatePlumeOffset.
    @param plumeProperties - Will on successful return be filled with the calculated properties of the plume (not completeness though).
    @param message - Will be filled with the reason the plume wasn't found, if it wasn't.. */
bool FindPlume(const std::vector<double>& scanAngles, const std::vector<double>& phi, const std::vector<double>& columns, const std::vector<double>& columnErrors, const std::vector<bool>& badEvaluation, long numPoints, double plumeOffset, CPlumeInScanProperty& plumeProperties, std::string* message = nullptr);

/** Finds the plume in the supplied scan. Return value is true if there is a plume, otherwise false
    @param evaluatedScan - the evaluated scan.
    @param specieIdx - the evaluated specie to get the columns for (normally zero).
    @param plumeOffset - the offset of the plume. Must have been calculated by a prior call to CalculatePlumeOffset.
    @param plumeProperties - Will on successful return be filled with the calculated properties of the plume (not completeness though).
    @param message - Will be filled with the reason the plume wasn't found, if it wasn't.. */
bool FindPlume(const BasicScanEvaluationResult& evaluatedScan, int specieIdx, double plumeOffset, CPlumeInScanProperty& plumeProperties, std::string* message = nullptr);

/** Tries to calculate the completeness of the given scan.
    This will internally call FindPlume to verify the location of the plume and sets the associated fields in plume
    The completeness is 1.0 if the entire plume can be seen and 0.0 if the plume cannot be seen at all.
    @param scanAngles - the scanAngles for the measurements.
    @param columns - the slant columns for the measurements. Must be from a normal scan
    @param columnErrors - the corresponding slant column errors for the measurement.
    @param badEvaluation - the result of the quality judgement of the measured columns,
                            badEvaluation[i] = true means a bad value
    @param offset - the offset of the scan, as calculate by CalculatePlumeOffset.
    @param numPoints - the number of points in the scan. Must also be the length
                            of the vectors 'columns', 'columnErrors', and 'badEvaluation'
    @param completeness - Will on successful return be filled with the completeness of the plume.
    @param message - Will be filled with the reason the plume wasn't found, if it wasn't..
    @return true if there is a plume, otherwise false. */
bool CalculatePlumeCompleteness(
    const std::vector<double>& scanAngles,
    const std::vector<double>& phi,
    const std::vector<double>& columns,
    const std::vector<double>& columnErrors,
    const std::vector<bool>& badEvaluation,
    double offset,
    long numPoints,
    CPlumeInScanProperty& plumeProperties,
    std::string* message = nullptr);

/** Tries to calculate the completeness of the given scan.
    This will internally call FindPlume to verify the location of the plume.
    This needs the field 'offset' in plumeProperties to be set.
    The completeness is 1.0 if the entire plume can be seen and 0.0 if the plume
    cannot be seen at all.
    Return value is true if there is a plume, otherwise false
    @param evaluatedScan - the evaluated scan.
    @param specieIdx - the index of the main specie in 'evaluatedScan', typically 0.
    @param plumeProperties - will on successful return have its 'completness' field filled in.
    @param message - Will be filled with the reason the plume wasn't found, if it wasn't..
    @return true if there is a plume in the scan, otherwise false. */
bool CalculatePlumeCompleteness(const BasicScanEvaluationResult& evaluatedScan, int specieIdx, CPlumeInScanProperty& plumeProperties, std::string* message = nullptr);

/** Calculates the 'offset' of the scan, i.e. the column amount in the sky spectrum, by judging from
    the lowest columns in the scan. */
double CalculatePlumeOffset(const std::vector<double>& columns, const std::vector<bool>& badEvaluation, long numPoints);

/** Calculates the 'offset' of the scan, i.e. the column amount in the sky spectrum, by judging from
    the lowest columns in the scan.
    This value is filled into the provided CPlumeInScanProperty and returned. */
double CalculatePlumeOffset(const BasicScanEvaluationResult& evaluatedScan, int specieIdx, CPlumeInScanProperty& plumeProperties);

}
