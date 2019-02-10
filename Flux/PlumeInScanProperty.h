#pragma once

#include "../Defintions.h"
#include <vector>

/** The class <b>CPlumeInScanProperty</b> is used to describe
    how a plume is seen by a scan. This incorporates properties
    such as the angle at which the centre of the plume is seen
    or an estimate of how large portion of the plume is seen. */
class CPlumeInScanProperty
{
public:
    CPlumeInScanProperty() = default;
    ~CPlumeInScanProperty() = default;

    /** Assignment operator */
    CPlumeInScanProperty &operator=(const CPlumeInScanProperty &p2) = default;

    /** The centre position of the plume,
        one parameter for each motor */
    double m_plumeCentre[2] = { NOT_A_NUMBER, NOT_A_NUMBER };

    /** The estimated error in the estimation of the
        plume centre position. */
    double m_plumeCentreError[2] = { NOT_A_NUMBER , NOT_A_NUMBER };

    /** The offset of the scan, describes the column
        of gas in the sky-spectrum */
    double m_offset = 0.0;

    /** The completeness of the scan. This is a value
        ranging from 0.5 to 1.0 trying to estimate
        how large portion of the plume is visible
        in the scan. */
    double m_completeness = 0.0;

    /** The edges of the plume */
    double m_plumeEdge_low = NOT_A_NUMBER;
    double m_plumeEdge_high = NOT_A_NUMBER;
};

// --------------------------------------------------------------------
// -------------- CALCULATING IF WE SEE THE PLUME ---------------------
// --------------------------------------------------------------------

/** Finds the plume in the supplied scan. Return value is true if there is a plume, otherwise false
    @param scanAngles - the scanAngles for the measurements.
    @param columns - the slant columns for the measurements. Must be from a normal scan
    @param columnErrors - the corresponding slant column errors for the measurement.
    @param badEvaluation - the result of the quality judgement of the measured columns,
                        badEvaluation[i] = true means a bad value
    @param numPoints - the number of points in the scan. Must also be the length
                        of the vectors 'columns', 'columnErrors', and 'badEvaluation'
    @param plumeCentre - Will on successful return be filled with the scan angle
                        which hits the centre of the plume
    @param plumeWidth - will on successful return be filled with the
                        estimated width of the plume (same unit as the scanAngles)
    @param plumeEdge_low - will on successful return be filled with the
                            lower edge  of the plume (same unit as the scanAngles)
    @param plumeEdge_high - will on successful return be filled with the
                            higher edge  of the plume (same unit as the scanAngles)	*/
bool FindPlume(const std::vector<double>& scanAngles, const std::vector<double>& phi, const std::vector<double>& columns, const std::vector<double>& columnErrors, const std::vector<bool>& badEvaluation, long numPoints, double &plumeCentre_alpha, double &plumeCentre_phi, double &plumeEdge_low, double &plumeEdge_high);
bool FindPlume(const double *scanAngles, const double *phi, const double *columns, const double *columnErrors, const bool *badEvaluation, long numPoints, CPlumeInScanProperty &plumeProperties);

/** Tries to calculate the completeness of the given scan.
    The completeness is 1.0 if the entire plume can be seen and 0.0 if the plume
    cannot be seen at all.
    Return value is true if there is a plume, otherwise false
    @param scanAngles - the scanAngles for the measurements.
    @param columns - the slant columns for the measurements. Must be from a normal scan
    @param columnErrors - the corresponding slant column errors for the measurement.
    @param badEvaluation - the result of the quality judgement of the measured columns,
                            badEvaluation[i] = true means a bad value
    @param numPoints - the number of points in the scan. Must also be the length
                            of the vectors 'columns', 'columnErrors', and 'badEvaluation'
    @param completeness - Will on successful return be filled with the completeness of the plume */
bool CalculatePlumeCompleteness(const std::vector<double>& scanAngles, const std::vector<double>& phi, const std::vector<double>& columns, const std::vector<double>& columnErrors, const std::vector<bool>& badEvaluation, double offset, long numPoints, double &completeness);
bool CalculatePlumeCompleteness(const double *scanAngles, const double *phi, const double *columns, const double *columnErrors, const bool *badEvaluation, double offset, long numPoints, CPlumeInScanProperty &plumeProperties);

/** Calculates the 'offset' of the scan, i.e. the column amount in the sky spectrum, by judging from 
    the lowest columns in the scan. */
double CalculatePlumeOffset(const double *columns, const bool *badEvaluation, long numPoints);
