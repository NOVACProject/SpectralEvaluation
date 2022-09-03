#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>
#include <SpectralEvaluation/VectorUtils.h>
#include <SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h>
#include <SpectralEvaluation/VectorUtils.h>
#include <algorithm>
#include <cstring>
#include <sstream>

namespace novac
{

    template <class T> double Average(T array[], long nElements) {
        if (nElements <= 0)
            return 0.0;

        double sum = 0;
        for (int k = 0; k < nElements; ++k) {
            sum += array[k];
        }
        return (sum / nElements);
    }

    template <class T> T Min(T* pBuffer, long bufLen) {
        T minValue = pBuffer[0];
        for (long i = 1; i < bufLen; i++) {
            if (pBuffer[i] < minValue)
                minValue = pBuffer[i];
        }
        return minValue;
    }

    template <class T> T Max(T* pBuffer, long bufLen) {
        T maxValue = pBuffer[0];
        for (long i = 1; i < bufLen; ++i) {
            if (pBuffer[i] > maxValue)
                maxValue = pBuffer[i];
        }
        return maxValue;
    }

    /// --------------------------- FINDING THE PLUME IN ONE SCAN ---------------------------

    bool FindPlume(const BasicScanEvaluationResult& evaluatedScan, int specieIdx, CPlumeInScanProperty& plumeProperties, std::string* message)
    {
        std::vector<double> scanAngles;
        std::vector<double> phi;
        std::vector<double> columns;
        std::vector<double> columnErrors;
        std::vector<bool> badEvaluation;
        const long numPoints = static_cast<long>(evaluatedScan.m_spec.size());

        for (long idx = 0; idx < numPoints; ++idx)
        {
            scanAngles.push_back(evaluatedScan.m_specInfo[idx].m_scanAngle);
            phi.push_back(evaluatedScan.m_specInfo[idx].m_scanAngle2);

            columns.push_back(evaluatedScan.m_spec[idx].m_referenceResult[specieIdx].m_column);
            columnErrors.push_back(evaluatedScan.m_spec[idx].m_referenceResult[specieIdx].m_columnError);
            badEvaluation.push_back(evaluatedScan.m_spec[idx].IsBad());
        }

        return FindPlume(scanAngles, phi, columns, columnErrors, badEvaluation, numPoints, plumeProperties, message);
    }


    // VERSION 1: FROM NOVACPROGRAM
    bool FindPlume(
        const std::vector<double>& scanAngles,
        const std::vector<double>& phi,
        const std::vector<double>& columns,
        const std::vector<double>& columnErrors,
        const std::vector<bool>& badEvaluation,
        long numPoints,
        CPlumeInScanProperty& plumeProperties,
        std::string* message)
    {

        // There is a plume iff there is a region, where the column-values are considerably
        //	much higher than in the rest of the scan

        // Make a local copy of the values, picking out only the data from the good spectra.
        std::vector<double> columnOfGoodSpectra;
        std::vector<double> columnErrorOfGoodSpectra;
        std::vector<double> angleOfGoodSpectra;
        std::vector<double> angle2OfGoodSpectra;

        columnOfGoodSpectra.reserve(numPoints);
        columnErrorOfGoodSpectra.reserve(numPoints);
        angleOfGoodSpectra.reserve(numPoints);
        angle2OfGoodSpectra.reserve(numPoints);

        int numberOfGoodspectra = 0; // <-- the number of ok column values
        for (int k = 0; k < numPoints; ++k) {
            if (!badEvaluation[k]) {
                columnOfGoodSpectra.push_back(columns[k]);
                columnErrorOfGoodSpectra.push_back(columnErrors[k]);
                angleOfGoodSpectra.push_back(scanAngles[k]);
                angle2OfGoodSpectra.push_back(phi[k]);
                ++numberOfGoodspectra;
            }
        }
        if (numberOfGoodspectra <= 5) { // <-- if too few ok points, then there's no plume
            if (nullptr != message)
            {
                *message = "Plume not found, less than five spectra are labelled as good evaluations.";
            }
            return false;
        }

        // Try different divisions of the scan to see if there is a region of at least
        //  'minWidth' values where the column-value is considerably higher than the rest
        double highestDifference = -1e16; // maximum difference between in-plume and out-of-plume regions
        long minWidth = 5;
        int foundRegionLowIdx = 0;
        int foundRegionHighIdx = 0;
        for (int testedLowIdx = 0; testedLowIdx < numberOfGoodspectra; ++testedLowIdx)
        {
            for (int testedHighIdx = testedLowIdx + minWidth; testedHighIdx < numberOfGoodspectra; ++testedHighIdx)
            {
                const int testedRegionSize = testedHighIdx - testedLowIdx;

                // The width of the region has to be at least 'minWidth' values, otherwise there's no idea to search
                //  There must also be at least 'minWidth' values outside of the region...
                if ((testedRegionSize < minWidth) || (numberOfGoodspectra - testedRegionSize < minWidth))
                {
                    continue;
                }

                // the average column value in the region we're testing
                const double avgInRegion = Average(columnOfGoodSpectra.data() + testedLowIdx, testedRegionSize);

                // the average column value outside of the region we're testing
                const double avgOutRegion = (Average(columnOfGoodSpectra.data(), testedLowIdx) * testedLowIdx + Average(columnOfGoodSpectra.data() + testedHighIdx, numberOfGoodspectra - testedHighIdx) * (numberOfGoodspectra - testedHighIdx)) / (testedLowIdx + numberOfGoodspectra - testedHighIdx);

                if (avgInRegion - avgOutRegion > highestDifference)
                {
                    highestDifference = avgInRegion - avgOutRegion;
                    foundRegionLowIdx = testedLowIdx;
                    foundRegionHighIdx = testedHighIdx;
                }
            }
        }

        // Calculate the average column error, for the good measurement points
        const double avgColError = Average(columnErrorOfGoodSpectra.data(), numberOfGoodspectra);

        if (!(highestDifference > 5 * avgColError))
        {
            if (nullptr != message)
            {
                std::stringstream msg;
                msg << "Plume not found, strongest plume-to-background difference: " << highestDifference << ", average column error: " << avgColError << ", plume-to-background ratio: " << highestDifference / avgColError;
                *message = msg.str();
            }
            return false;
        }

        // the plume centre is the average of the scan-angles in the 'plume-region'
        //	weighted with the column values
        double sumAngle_alpha = 0, sumAngle_phi = 0, sumWeight = 0;
        for (int k = foundRegionLowIdx; k < foundRegionHighIdx; ++k)
        {
            sumAngle_alpha += angleOfGoodSpectra[k] * columnOfGoodSpectra[k];
            sumAngle_phi += angle2OfGoodSpectra[k] * columnOfGoodSpectra[k];
            sumWeight += columnOfGoodSpectra[k];
        }
        plumeProperties.plumeCenter = sumAngle_alpha / sumWeight;
        plumeProperties.plumeCenter2 = sumAngle_phi / sumWeight;

        // Add a reasonability check here. The plume center must never fall outside of the range of angles but may do so
        // if there are 'bad' column values included in 'columnOfGoodSpectra'.
        plumeProperties.plumeCenter = std::max(angleOfGoodSpectra[foundRegionLowIdx], std::min(angleOfGoodSpectra[foundRegionHighIdx], plumeProperties.plumeCenter));
        plumeProperties.plumeCenter2 = std::max(angle2OfGoodSpectra[foundRegionLowIdx], std::min(angle2OfGoodSpectra[foundRegionHighIdx], plumeProperties.plumeCenter));

        // The edges of the plume
        // TODO: this could really be better, with interpolation between the angles...
        plumeProperties.plumeEdgeLow = angleOfGoodSpectra.front();
        plumeProperties.plumeEdgeHigh = angleOfGoodSpectra.back();
        double peakLow = angleOfGoodSpectra.front();
        double peakHigh = angleOfGoodSpectra.back();
        const double minCol = Min(columnOfGoodSpectra.data(), numberOfGoodspectra);
        const double maxCol = Max(columnOfGoodSpectra.data(), numberOfGoodspectra) - minCol;
        const double maxCol_div_e = maxCol * 0.3679;
        const double maxCol_90 = maxCol * 0.90;
        const double maxCol_half = maxCol * 0.5;

        // The left size of the plume
        for (int idx = 0; idx < numberOfGoodspectra - 1; ++idx)
        {
            if (angleOfGoodSpectra[idx] > plumeProperties.plumeCenter)
            {
                break;
            }
            if ((columnOfGoodSpectra[idx] - minCol) < maxCol_div_e)
            {
                plumeProperties.plumeEdgeLow = angleOfGoodSpectra[idx];
            }
            if ((columnOfGoodSpectra[idx] - minCol) < maxCol_half)
            {
                plumeProperties.plumeHalfLow = angleOfGoodSpectra[idx];
            }
            if (((columnOfGoodSpectra[idx] - minCol) < maxCol_90) && ((columnOfGoodSpectra[idx + 1] - minCol) >= maxCol_90))
            {
                peakLow = angleOfGoodSpectra[idx];
            }
        }

        // The right side of the plume
        for (int idx = numberOfGoodspectra - 1; idx > 0; --idx)
        {
            if (angleOfGoodSpectra[idx] <= plumeProperties.plumeCenter)
            {
                break;
            }
            if ((columnOfGoodSpectra[idx] - minCol) < maxCol_div_e)
            {
                plumeProperties.plumeEdgeHigh = angleOfGoodSpectra[idx];
            }
            if ((columnOfGoodSpectra[idx] - minCol) < maxCol_half)
            {
                plumeProperties.plumeHalfHigh = angleOfGoodSpectra[idx];
            }
            if (((columnOfGoodSpectra[idx] - minCol) < maxCol_90) && ((columnOfGoodSpectra[idx - 1] - minCol) >= maxCol_90))
            {
                peakHigh = angleOfGoodSpectra[idx];
            }
        }

        plumeProperties.plumeCenterError = (peakHigh - peakLow) / 2;

        return true;
    }

    /// --------------------------- PLUME COMPLETENESS ---------------------------

    bool CalculatePlumeCompleteness(const BasicScanEvaluationResult& evaluatedScan, int specieIdx, CPlumeInScanProperty& plumeProperties, std::string* message)
    {
        std::vector<double> scanAngles;
        std::vector<double> phi;
        std::vector<double> columns;
        std::vector<double> columnErrors;
        std::vector<bool> badEvaluation;
        const long numPoints = static_cast<long>(evaluatedScan.m_spec.size());

        for (long idx = 0; idx < numPoints; ++idx)
        {
            scanAngles.push_back(evaluatedScan.m_specInfo[idx].m_scanAngle);
            phi.push_back(evaluatedScan.m_specInfo[idx].m_scanAngle2);

            columns.push_back(evaluatedScan.m_spec[idx].m_referenceResult[specieIdx].m_column);
            columnErrors.push_back(evaluatedScan.m_spec[idx].m_referenceResult[specieIdx].m_columnError);
            badEvaluation.push_back(evaluatedScan.m_spec[idx].IsBad());
        }

        return CalculatePlumeCompleteness(scanAngles, phi, columns, columnErrors, badEvaluation, plumeProperties.offset, numPoints, plumeProperties, message);
    }

    // VERSION 1: FROM NOVACPROGRAM
    bool CalculatePlumeCompleteness(
        const std::vector<double>& scanAngles,
        const std::vector<double>& phi,
        const std::vector<double>& columns,
        const std::vector<double>& columnErrors,
        const std::vector<bool>& badEvaluation,
        double offset,
        long numPoints,
        CPlumeInScanProperty& plumeProperties,
        std::string* message)
    {

        int nDataPointsToAverage = 5;

        // Check if there is a plume at all...
        bool inPlume = ::novac::FindPlume(scanAngles, phi, columns, columnErrors, badEvaluation, numPoints, plumeProperties, message);
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

    double CalculatePlumeOffset(const BasicScanEvaluationResult& evaluatedScan, int specieIdx, CPlumeInScanProperty& plumeProperties)
    {
        std::vector<double> columns;
        std::vector<bool> badEvaluation;
        const long numPoints = static_cast<long>(evaluatedScan.m_spec.size());

        for (long idx = 0; idx < numPoints; ++idx)
        {
            columns.push_back(evaluatedScan.m_spec[idx].m_referenceResult[specieIdx].m_column);
            badEvaluation.push_back(evaluatedScan.m_spec[idx].IsBad());
        }

        plumeProperties.offset = CalculatePlumeOffset(columns, badEvaluation, numPoints);

        return plumeProperties.offset;
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

}
