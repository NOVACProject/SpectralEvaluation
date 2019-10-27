#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>
#include <SpectralEvaluation/VectorUtils.h>
#include <algorithm>
#include <cstring>
#include <sstream>

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

/// --------------------------- FINDING THE PLUME IN ONE SCAN ---------------------------

// VERSION 1: FROM NOVACPROGRAM
bool FindPlume(const std::vector<double>& scanAngles, const std::vector<double>& phi, const std::vector<double>& columns, const std::vector<double>& columnErrors, const std::vector<bool>& badEvaluation, long numPoints, CPlumeInScanProperty& plumeProperties, std::string* message)
{

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
        if (nullptr != message)
        {
            *message = "Plume not found, less than five spectra are labelled as good evaluations.";
        }
        return false;
    }

    // Try different divisions of the scan to see if there is a region of at least
    //	'minWidth' values where the column-value is considerably higher than the rest
    double highestDifference = -1e16; // maximum difference between in-plume and out-of-plume regions
    long minWidth = 5;
    int foundRegionLowIdx = 0;
    int foundRegionHighIdx = 0;
    for (int testedLowIdx = 0; testedLowIdx < nCol; ++testedLowIdx)
    {
        for (int testedHighIdx = testedLowIdx + minWidth; testedHighIdx < nCol; ++testedHighIdx)
        {
            const int testedRegionSize = testedHighIdx - testedLowIdx;

            // The width of the region has to be at least 'minWidth' values, otherwise there's no idea to search
            //  There must also be at least 'minWidth' values outside of the region...
            if ((testedRegionSize < minWidth) || (nCol - testedRegionSize < minWidth))
            {
                continue;
            }

            // the average column value in the region we're testing
            const double avgInRegion = Average(col.data() + testedLowIdx, testedRegionSize);

            // the average column value outside of the region we're testing
            // TODO: These are actually 'sum'
            const double avgOutRegion = (Average(col.data(), testedLowIdx)*testedLowIdx + Average(col.data() + testedHighIdx, nCol - testedHighIdx)*(nCol - testedHighIdx)) / (testedLowIdx + nCol - testedHighIdx);

            if (avgInRegion - avgOutRegion > highestDifference)
            {
                highestDifference = avgInRegion - avgOutRegion;
                foundRegionLowIdx = testedLowIdx;
                foundRegionHighIdx = testedHighIdx;
            }
        }
    }

    // Calculate the average column error, for the good measurement points
    double avgColError = Average(colE.data(), nCol);

    if (highestDifference > 5 * avgColError)
    {
        // the plume centre is the average of the scan-angles in the 'plume-region'
        //	weighted with the column values
        double sumAngle_alpha = 0, sumAngle_phi = 0, sumWeight = 0;
        for (int k = foundRegionLowIdx; k < foundRegionHighIdx; ++k)
        {
            sumAngle_alpha += angle[k] * col[k];
            sumAngle_phi += p[k] * col[k];
            sumWeight += col[k];
        }
        plumeProperties.plumeCenter = sumAngle_alpha / sumWeight;
        plumeProperties.plumeCenter2 = sumAngle_phi / sumWeight; // if phi == NULL then this will be non-sense

        // The edges of the plume
        // TODO: this could really be better, with interpolation between the angles...
        plumeProperties.plumeEdgeLow = angle[0];
        plumeProperties.plumeEdgeHigh = angle[nCol - 1];
        double peakLow = angle[0];
        double peakHigh = angle[nCol - 1];
        double minCol = Min(col.data(), nCol);
        double maxCol_div_e = (Max(col.data(), nCol) - minCol) * 0.3679;
        double maxCol_90 = (Max(col.data(), nCol) - minCol) * 0.90;
        double maxCol_half = (Max(col.data(), nCol) - minCol) * 0.5;
        for (int k = 0; k < nCol; ++k)
        {
            if (angle[k] > plumeProperties.plumeCenter)
            {
                break;
            }
            if ((col[k] - minCol) < maxCol_div_e)
            {
                plumeProperties.plumeEdgeLow = angle[k];
            }
            if ((col[k] - minCol) < maxCol_half)
            {
                plumeProperties.plumeHalfLow = angle[k];
            }
            if (((col[k] - minCol) < maxCol_90) && ((col[k + 1] - minCol) >= maxCol_90))
            {
                peakLow = angle[k];
            }
        }
        for (int k = nCol - 1; k >= 0; --k)
        {
            if (angle[k] <= plumeProperties.plumeCenter)
            {
                break;
            }
            if ((col[k] - minCol) < maxCol_div_e)
            {
                plumeProperties.plumeEdgeHigh = angle[k];
            }
            if ((col[k] - minCol) < maxCol_half)
            {
                plumeProperties.plumeHalfHigh = angle[k];
            }
            if (((col[k] - minCol) < maxCol_90) && ((col[k - 1] - minCol) >= maxCol_90))
            {
                peakHigh = angle[k];
            }
        }

        plumeProperties.plumeCenterError = (peakHigh - peakLow) / 2;

        return true;
    }
    else
    {
        if (nullptr != message)
        {
            std::stringstream msg;
            msg << "Plume not found, strongest plume-to-background difference: " << highestDifference << ", average column error: " << avgColError << ", plume-to-background ratio: " << highestDifference / avgColError;
            *message = msg.str();
        }
        return false;
    }
}

/// --------------------------- PLUME COMPLETENESS ---------------------------

// VERSION 1: FROM NOVACPROGRAM
bool CalculatePlumeCompleteness(const std::vector<double>& scanAngles, const std::vector<double>& phi, const std::vector<double>& columns, const std::vector<double>& columnErrors, const std::vector<bool>& badEvaluation, double offset, long numPoints, CPlumeInScanProperty &plumeProperties, std::string* message)
{

    int nDataPointsToAverage = 5;

    // Check if there is a plume at all...
    bool inPlume = FindPlume(scanAngles, phi, columns, columnErrors, badEvaluation, numPoints, plumeProperties, message);
    if (!inPlume) {
        plumeProperties.completeness = 0.0; // <-- no plume at all
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
        plumeProperties.completeness = 0.0; // <-- no plume at all
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
        plumeProperties.completeness = 0.0; // <-- no plume at all
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
    plumeProperties.completeness = 1.0 - 0.5 * std::max(avgLeft, avgRight) / maxColumn;
    if (plumeProperties.completeness > 1.0)
        plumeProperties.completeness = 1.0;

    return true;
}

double CalculatePlumeOffset(const std::vector<double>& columns, const std::vector<bool>& badEvaluation, long numPoints)
{
    // calculate the offset as the average of the three lowest column values 
    //    that are not considered as 'bad' values
    std::vector<double> testColumns;
    testColumns.reserve(columns.size());

    for (int i = 0; i < numPoints; ++i)
    {
        if (!badEvaluation[i])
        {
            testColumns.push_back(columns[i]);
        }
    }

    if (testColumns.size() <= 5)
    {
        return 0.0;
    }

    // Find the N lowest column values
    const int N = (int)(0.2 * testColumns.size());
    std::vector<double> m(N, 1e6);
    FindNLowest(testColumns, N, m);
    const double avg = Average(m.data(), N);
    return avg;
}
