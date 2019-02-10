#include "PlumeInScanProperty.h"
#include <algorithm>

template <class T> double Average(T array[], long nElements) {
    if (nElements <= 0)
        return 0.0;

    double sum = 0;
    for (int k = 0; k < nElements; ++k) {
        sum += array[k];
    }
    return (sum / nElements);
}

template <class T> T Min(T *pBuffer, long bufLen) {
    T minValue = pBuffer[0];
    for (long i = 1; i < bufLen; i++) {
        if (pBuffer[i] < minValue)
            minValue = pBuffer[i];
    }
    return minValue;
}

template <class T> T Max(T *pBuffer, long bufLen) {
    T maxValue = pBuffer[0];
    for (long i = 1; i < bufLen; ++i) {
        if (pBuffer[i] > maxValue)
            maxValue = pBuffer[i];
    }
    return maxValue;
}

/** This function finds the 'N' lowest values in the supplied array.
On successfull return the array 'output' will be filled with the N lowest
values in the input array, sorted in ascending order.
@param array - the array to look into
@param nElements - the number of elements in the supplied array
@param output - the output array, must be at least 'N'- elements long
@param N - the number of values to take out. 	*/
template <class T> bool FindNLowest(const T array[], long nElements, T output[], int N, int *indices = NULL) {
    for (int i = 0; i < N; ++i)
        output[i] = 1e16; // to get some initial value

                          // loop through all elements in the array
    for (int i = 0; i < nElements; ++i) {

        // compare this element with all elements in the output array.
        for (int j = 0; j < N; ++j) {
            if (array[i] < output[j]) {
                // If we found a higher value, shift all other values down one step...
                for (int k = N - 1; k > j; --k) {
                    output[k] = output[k - 1];
                    if (indices)							indices[k] = indices[k - 1];
                }
                output[j] = array[i];
                if (indices)								indices[j] = i;
                break;
            }
        }
    }
    return true;
}

/// --------------------------- FINDING THE PLUME IN ONE SCAN ---------------------------

// VERSION 1: FROM NOVACPROGRAM
bool FindPlume(const std::vector<double>& scanAngles, const std::vector<double>& phi, const std::vector<double>& columns, const std::vector<double>& columnErrors, const std::vector<bool>& badEvaluation, long numPoints, double &plumeCentre_alpha, double &plumeCentre_phi, double &plumeEdge_low, double &plumeEdge_high) {

    // There is a plume iff there is a region, where the column-values are considerably
    //	much higher than in the rest of the scan

    // Make a local copy of the values, picking out only the good ones
    std::vector<double> col(numPoints);
    std::vector<double> colE(numPoints);
    std::vector<double> angle(numPoints);
    std::vector<double> p(numPoints);

    int		nCol = 0; // <-- the number of ok column values
    for (int k = 0; k < numPoints; ++k) {
        if (!badEvaluation[k]) {
            col[nCol] = columns[k];
            colE[nCol] = columnErrors[k];
            angle[nCol] = scanAngles[k];
            p[nCol] = phi[k];
            ++nCol;
        }
    }
    if (nCol <= 5) { // <-- if too few ok points, then there's no plume
        return false;
    }

    // Try different divisions of the scan to see if there is a region of at least
    //	'minWidth' values where the column-value is considerably higher than the rest
    double highestDifference = -1e16;
    long minWidth = 5;
    int regionLow = 0, regionHigh = 0;
    for (int low = 0; low < nCol; ++low) {
        for (int high = low + minWidth; high < nCol; ++high) {
            // The width of the region has to be at least 'minWidth' values, otherwise there's no idea to search
            if (high - low < minWidth)
                continue;

            // the average column value in the region we're testing
            double avgInRegion = Average(col.data() + low, high - low);

            // the average column value outside of the region we're testing
            double avgOutRegion = (Average(col.data(), low) + Average(col.data() + high, nCol - high)) * 0.5;

            if (avgInRegion - avgOutRegion > highestDifference) {
                highestDifference = avgInRegion - avgOutRegion;
                regionLow = low;
                regionHigh = high;
            }
        }
    }

    // Calculate the average column error, for the good measurement points
    double avgColError = Average(colE.data(), nCol);

    if (highestDifference > 5 * avgColError) {
        // the plume centre is the average of the scan-angles in the 'plume-region'
        //	weighted with the column values
        double sumAngle_alpha = 0, sumAngle_phi = 0, sumWeight = 0;
        for (int k = regionLow; k < regionHigh; ++k) {
            sumAngle_alpha += angle[k] * col[k];
            sumAngle_phi += p[k] * col[k];
            sumWeight += col[k];
        }
        plumeCentre_alpha = sumAngle_alpha / sumWeight;
        plumeCentre_phi = sumAngle_phi / sumWeight; // if phi == NULL then this will be non-sense

                                                    // The edges of the plume
        plumeEdge_low = angle[0];
        plumeEdge_high = angle[nCol - 1];
        double minCol = Min(col.data(), nCol);
        double maxCol_div_e = (Max(col.data(), nCol) - minCol) * 0.3679;
        for (int k = 0; k < nCol; ++k) {
            if (angle[k] > plumeCentre_alpha) {
                break;
            }
            else if ((col[k] - minCol) < maxCol_div_e) {
                plumeEdge_low = angle[k];
            }
        }
        for (int k = nCol - 1; k >= 0; --k) {
            if (angle[k] <= plumeCentre_alpha) {
                break;
            }
            else if ((col[k] - minCol) < maxCol_div_e) {
                plumeEdge_high = angle[k];
            }
        }

        return true;
    }
    else {
        return false;
    }
}

// VERSION 2: FROM NOVACPPP
bool FindPlume(const double *scanAngles, const double *phi, const double *columns, const double *columnErrors, const bool *badEvaluation, long numPoints, CPlumeInScanProperty &plumeProperties) {

    // There is a plume iff there is a region, where the column-values are considerably
    //	much higher than in the rest of the scan

    // Make a local copy of the values, picking out only the good ones
    double *col = new double[numPoints];
    double *colE = new double[numPoints];
    double *angle = new double[numPoints];
    double *p = new double[numPoints];
    double minColumn = 1e23;
    int		nCol = 0; // <-- the number of ok column values
    for (int k = 0; k < numPoints; ++k) {
        if (!badEvaluation[k]) {
            col[nCol] = columns[k];
            colE[nCol] = columnErrors[k];
            angle[nCol] = scanAngles[k];
            p[nCol] = phi[k];
            minColumn = std::min(minColumn, columns[k]);
            ++nCol;
        }
    }
    if (nCol <= 5) { // <-- if too few ok points, then there's no plume
        delete[] col;		delete[] colE;		delete[] angle;	delete[] p;
        return false;
    }

    // Add the smallest column-value to all other columns, this to make an approximate offset
    for (int k = 0; k < nCol; ++k) {
        col[k] -= minColumn;
    }

    // Try different divisions of the scan to see if there is a region of at least
    //	'minWidth' values where the column-value is considerably higher than the rest
    double highestDifference = -1e16;
    long minWidth = 5;
    int regionLow = 0, regionHigh = 0;
    for (int low = 0; low < nCol; ++low) {
        for (int high = low + minWidth; high < nCol; ++high) {
            // The width of the region has to be at least 'minWidth' values, otherwise there's no idea to search
            //	There must also be at least 'minWidth' values outside of the region...
            if ((high - low < minWidth) || (nCol - high + low < minWidth))
                continue;

            // the average column value in the region we're testing
            double avgInRegion = Average(col + low, high - low);

            // the average column value outside of the region we're testing
            double avgOutRegion = (Average(col, low)*low + Average(col + high, nCol - high)*(nCol - high)) / (low + nCol - high);

            if (avgInRegion - avgOutRegion > highestDifference) {
                highestDifference = avgInRegion - avgOutRegion;
                regionLow = low;
                regionHigh = high;
            }
        }
    }

    // Calculate the average column error, for the good measurement points
    double avgColError = Average(colE, nCol);

    if (highestDifference > 5 * avgColError) {
        // the plume centre is the average of the scan-angles in the 'plume-region'
        //	weighted with the column values
        double sumAngle_alpha = 0, sumAngle_phi = 0, sumWeight = 0;
        for (int k = regionLow; k < regionHigh; ++k) {
            sumAngle_alpha += angle[k] * col[k];
            sumAngle_phi += p[k] * col[k];
            sumWeight += col[k];
        }
        plumeProperties.m_plumeCentre[0] = sumAngle_alpha / sumWeight;
        plumeProperties.m_plumeCentre[1] = sumAngle_phi / sumWeight; // if phi == NULL then this will be non-sense

                                                                     // The edges of the plume and the error in the estimated plume centres
        plumeProperties.m_plumeEdge_low = angle[0];
        plumeProperties.m_plumeEdge_high = angle[nCol - 1];
        double peakLow = angle[0];
        double peakHigh = angle[nCol - 1];
        double maxCol_div_e = Max(col, nCol) * 0.3679;
        double maxCol_90 = Max(col, nCol) * 0.90;

        // Search for the lower edge of the plume ...
        for (int k = 0; k < nCol - 2; ++k) {
            if (angle[k] > plumeProperties.m_plumeCentre[0]) {
                break;
            }
            else if (Average(col + k, 2) < maxCol_div_e) {
                plumeProperties.m_plumeEdge_low = angle[k];
            }
            if ((col[k] < maxCol_90) && (col[k + 1] >= maxCol_90)) {
                peakLow = angle[k];
            }
        }

        // .. and then the upper edge
        for (int k = nCol - 1; k >= 1; --k) {
            if (angle[k] <= plumeProperties.m_plumeCentre[0]) {
                break;
            }
            else if (Average(col + k - 1, 2) < maxCol_div_e) {
                plumeProperties.m_plumeEdge_high = angle[k];
            }
            if ((col[k] < maxCol_90) && (col[k - 1] >= maxCol_90)) {
                peakHigh = angle[k];
            }
        }
        plumeProperties.m_plumeCentreError[0] = (peakHigh - peakLow) / 2;

        delete[] col;	delete[] colE;
        delete[] angle; delete[] p;
        return true;
    }
    else {
        delete[] col;	delete[] colE;
        delete[] angle; delete[] p;
        return false;
    }
}


/// --------------------------- PLUME COMPLETENESS ---------------------------

// VERSION 1: FROM NOVACPROGRAM
bool CalculatePlumeCompleteness(const std::vector<double>& scanAngles, const std::vector<double>& phi, const std::vector<double>& columns, const std::vector<double>& columnErrors, const std::vector<bool>& badEvaluation, double offset, long numPoints, double &completeness) {
    double plumeCentre_alpha, plumeCentre_phi;
    double plumeEdge_low, plumeEdge_high;

    int nDataPointsToAverage = 5;

    // Check if there is a plume at all...
    bool inPlume = FindPlume(scanAngles, phi, columns, columnErrors, badEvaluation, numPoints, plumeCentre_alpha, plumeCentre_phi, plumeEdge_low, plumeEdge_high);
    if (!inPlume) {
        completeness = 0.0; // <-- no plume at all
        return false;
    }

    // Calculate the average of the 'nDataPointsToAverage' left-most values
    double avgLeft = 0.0;
    int nAverage = 0;
    for (int k = 0; k < numPoints; ++k) {
        if (!badEvaluation[k]) {
            avgLeft += columns[k] - offset;
            ++nAverage;
            if (nAverage == nDataPointsToAverage)
                break;
        }
    }
    if (nAverage < nDataPointsToAverage) {
        // not enough data-points to make an ok average, return fail
        completeness = 0.0; // <-- no plume at all
        return false;
    }
    avgLeft /= nDataPointsToAverage;

    // Calculate the average of the 'nDataPointsToAverage' right-most values
    double avgRight = 0.0;
    nAverage = 0;
    for (int k = numPoints - 1; k > 0; --k) {
        if (!badEvaluation[k]) {
            avgRight += columns[k] - offset;
            ++nAverage;
            if (nAverage == nDataPointsToAverage)
                break;
        }
    }
    if (nAverage < nDataPointsToAverage) {
        // not enough data-points to make an ok average, return fail
        completeness = 0.0; // <-- no plume at all
        return false;
    }
    avgRight /= nDataPointsToAverage;

    // Find the maximum column value
    double maxColumn = 0.0;
    for (int k = 0; k < numPoints; ++k) {
        if (!badEvaluation[k]) {
            maxColumn = std::max(maxColumn, columns[k] - offset);
        }
    }

    // The completeness
    completeness = 1.0 - 0.5 * std::max(avgLeft, avgRight) / maxColumn;
    if (completeness > 1.0)
        completeness = 1.0;

    return true;
}

// VERSION 2: FROM NOVACPPP
bool CalculatePlumeCompleteness(const double *scanAngles, const double *phi, const double *columns, const double *columnErrors, const bool *badEvaluation, double offset, long numPoints, CPlumeInScanProperty &plumeProperties) {
    int nDataPointsToAverage = 5;

    // Check if there is a plume at all...
    bool inPlume = FindPlume(scanAngles, phi, columns, columnErrors, badEvaluation, numPoints, plumeProperties);
    if (!inPlume) {
        plumeProperties.m_completeness = 0.0; // <-- no plume at all
        return false;
    }

    // Calculate the average of the 'nDataPointsToAverage' left-most values
    double avgLeft = 0.0;
    int nAverage = 0;
    for (int k = 0; k < numPoints; ++k) {
        if (!badEvaluation[k]) {
            avgLeft += columns[k] - offset;
            ++nAverage;
            if (nAverage == nDataPointsToAverage)
                break;
        }
    }
    if (nAverage < nDataPointsToAverage) {
        // not enough data-points to make an ok average, return fail
        plumeProperties.m_completeness = 0.0; // <-- no plume at all
        return false;
    }
    avgLeft /= nDataPointsToAverage;

    // Calculate the average of the 'nDataPointsToAverage' right-most values
    double avgRight = 0.0;
    nAverage = 0;
    for (int k = numPoints - 1; k > 0; --k) {
        if (!badEvaluation[k]) {
            avgRight += columns[k] - offset;
            ++nAverage;
            if (nAverage == nDataPointsToAverage)
                break;
        }
    }
    if (nAverage < nDataPointsToAverage) {
        // not enough data-points to make an ok average, return fail
        plumeProperties.m_completeness = 0.0; // <-- no plume at all
        return false;
    }
    avgRight /= nDataPointsToAverage;

    // Find the maximum column value
    double maxColumn = 0.0;
    for (int k = 0; k < numPoints; ++k) {
        if (!badEvaluation[k]) {
            maxColumn = std::max(maxColumn, columns[k] - offset);
        }
    }

    // The completeness
    plumeProperties.m_completeness = 1.0 - 0.5 * std::max(avgLeft, avgRight) / maxColumn;
    if (plumeProperties.m_completeness > 1.0)
        plumeProperties.m_completeness = 1.0;

    return true;
}


double CalculatePlumeOffset(const double *columns, const bool *badEvaluation, long numPoints) {
    double avg;
    long i;

    // calculate the offset as the average of the three lowest column values 
    //    that are not considered as 'bad' values

    double *testColumns = new double[numPoints];

    int numColumns = 0;
    for (i = 0; i < numPoints; ++i) {
        if (badEvaluation[i])
            continue;

        testColumns[numColumns++] = columns[i];
    }

    if (numColumns <= 5) {
        delete[] testColumns;
        return 0.0;
    }

    // Find the N lowest column values
    int N = (int)(0.2 * numColumns);
    double *m = new double[N];
    memset(m, (int)1e6, N * sizeof(double));
    if (FindNLowest(testColumns, numColumns, m, N)) {
        avg = Average(m, N);
        delete[] testColumns;
        delete[] m;
        return avg;
    }

    delete[] testColumns;
    delete[] m;

    // could not calculate a good offset.
    return 0;
}