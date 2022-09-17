#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>
#include <SpectralEvaluation/VectorUtils.h>
#include <SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h>
#include <SpectralEvaluation/VectorUtils.h>
#include <algorithm>
#include <cstring>
#include <sstream>

namespace novac
{
    // Helper structure, collecting together data which we need for selecting spectra in the scan.
    struct ScanEvaluationData
    {
        int indexInScan = 0;
        double scanAngle = 0.0; // the angle of the first motor, in degrees.
        double scanAngle2 = 0.0; // the angle of the second motor, only for Mark-II (heidelberg) instruments
        double offsetCorrectedColumn = 0.0;
        double columnError = 0.0;
        double peakSaturation = 0.0;
    };

    double AverageColumn(const std::vector< ScanEvaluationData>& data, int startIdx, int numberOfElements)
    {
        if (numberOfElements <= 0)
        {
            return 0.0;
        }

        double sum = 0;
        for (int k = startIdx; k < startIdx + numberOfElements; ++k)
        {
            sum += data[k].offsetCorrectedColumn;
        }

        return (sum / numberOfElements);
    }

    double AverageColumnError(const std::vector< ScanEvaluationData>& data, int startIdx, int numberOfElements)
    {
        if (numberOfElements <= 0)
        {
            return 0.0;
        }

        double sum = 0;
        for (int k = startIdx; k < startIdx + numberOfElements; ++k)
        {
            sum += data[k].columnError;
        }

        return (sum / numberOfElements);
    }

    double MinColumn(const std::vector< ScanEvaluationData>& data)
    {
        double minValue = data[0].offsetCorrectedColumn;
        for (const auto& eval : data)
        {
            minValue = std::min(minValue, eval.offsetCorrectedColumn);
        }
        return minValue;
    }

    double MaxColumn(const std::vector< ScanEvaluationData>& data)
    {
        double minValue = data[0].offsetCorrectedColumn;
        for (const auto& eval : data)
        {
            minValue = std::max(minValue, eval.offsetCorrectedColumn);
        }
        return minValue;
    }

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

    /// --------------------------- CALCULATING THE PLUME OFFSET ---------------------------

    // Calculates the offset of the plume given all the evaluations of the 'good' spectra in the scan.
    double CalculatePlumeOffsetFromGoodColumnValues(const std::vector<double>& goodColumns)
    {
        // Find the N lowest goodColumns values
        const int N = (int)(0.2 * goodColumns.size());
        std::vector<double> m(N, 1e6);
        FindNLowest(goodColumns, N, m);
        return Average(m.data(), N);
    }

    double CalculatePlumeOffset(const std::vector< ScanEvaluationData>& evaluation)
    {
        // calculate the offset as the average of the three lowest offsetCorrectedColumn values 
        //    that are not considered as 'bad' values
        std::vector<double> columns;
        columns.reserve(evaluation.size());

        for (const auto& eval : evaluation)
        {
            columns.push_back(eval.offsetCorrectedColumn);
        }

        return CalculatePlumeOffsetFromGoodColumnValues(columns);
    }

    double CalculatePlumeOffset(const BasicScanEvaluationResult& evaluatedScan, int specieIdx, CPlumeInScanProperty& plumeProperties)
    {
        std::vector<double> columns;
        const long numPoints = static_cast<long>(evaluatedScan.m_spec.size());

        for (long idx = 0; idx < numPoints; ++idx)
        {
            if (!evaluatedScan.m_spec[idx].IsBad())
            {
                columns.push_back(evaluatedScan.m_spec[idx].m_referenceResult[specieIdx].m_column);
            }
        }

        plumeProperties.offset = CalculatePlumeOffsetFromGoodColumnValues(columns);

        return plumeProperties.offset;
    }

    double CalculatePlumeOffset(const std::vector<double>& columns, const std::vector<bool>& badEvaluation, long numPoints)
    {
        // calculate the offset as the average of the three lowest offsetCorrectedColumn values 
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

        // Find the N lowest offsetCorrectedColumn values
        const int N = (int)(0.2 * testColumns.size());
        std::vector<double> m(N, 1e6);
        FindNLowest(testColumns, N, m);
        const double avg = Average(m.data(), N);
        return avg;
    }

    /// --------------------------- FINDING THE PLUME IN ONE SCAN ---------------------------

    bool FindPlume(const std::vector< ScanEvaluationData>& evaluation, CPlumeInScanProperty& plumeProperties, std::string* message)
    {
        const int numberOfGoodspectra = static_cast<int>(evaluation.size());

        if (numberOfGoodspectra <= 5) { // <-- if too few ok points, then there's no plume
            if (nullptr != message)
            {
                *message = "Plume not found, less than five spectra are labelled as good evaluations.";
            }
            return false;
        }

        // If the plume offset hasn't already been calculated then do so now.
        if (std::abs(plumeProperties.offset) < 1.0)
        {
            plumeProperties.offset = CalculatePlumeOffset(evaluation);
        }

        // Try different divisions of the scan to see if there is a region of at least
        //  'minWidth' values where the offsetCorrectedColumn-value is considerably higher than the rest
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

                // the average offsetCorrectedColumn value in the region we're testing
                const double avgInRegion = AverageColumn(evaluation, testedLowIdx, testedRegionSize);

                // the average offsetCorrectedColumn value outside of the region we're testing
                const double avgOutRegion = (AverageColumn(evaluation, 0, testedLowIdx) * testedLowIdx + AverageColumn(evaluation, testedHighIdx, numberOfGoodspectra - testedHighIdx) * (numberOfGoodspectra - testedHighIdx)) / (testedLowIdx + numberOfGoodspectra - testedHighIdx);

                if (avgInRegion - avgOutRegion > highestDifference)
                {
                    highestDifference = avgInRegion - avgOutRegion;
                    foundRegionLowIdx = testedLowIdx;
                    foundRegionHighIdx = testedHighIdx;
                }
            }
        }

        // Calculate the average offsetCorrectedColumn error, for the good measurement points
        const double avgColError = AverageColumnError(evaluation, 0, numberOfGoodspectra);

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
        //	weighted with the offsetCorrectedColumn values
        double sumAngle_alpha = 0, sumAngle_phi = 0, sumWeight = 0;
        for (int k = foundRegionLowIdx; k < foundRegionHighIdx; ++k)
        {
            sumAngle_alpha += evaluation[k].scanAngle * evaluation[k].offsetCorrectedColumn;
            sumAngle_phi += evaluation[k].scanAngle2 * evaluation[k].offsetCorrectedColumn;
            sumWeight += evaluation[k].offsetCorrectedColumn;
        }
        plumeProperties.plumeCenter = sumAngle_alpha / sumWeight;
        plumeProperties.plumeCenter2 = sumAngle_phi / sumWeight;

        // Add a reasonability check here. The plume center must never fall outside of the range of angles but may do so
        // if there are 'bad' offsetCorrectedColumn values included in 'columnOfGoodSpectra'.
        plumeProperties.plumeCenter = std::max(evaluation[foundRegionLowIdx].scanAngle, std::min(evaluation[foundRegionHighIdx].scanAngle, plumeProperties.plumeCenter));
        plumeProperties.plumeCenter2 = std::max(evaluation[foundRegionLowIdx].scanAngle2, std::min(evaluation[foundRegionHighIdx].scanAngle2, plumeProperties.plumeCenter));

        // The edges of the plume
        // TODO: this could really be better, with interpolation between the angles...
        plumeProperties.plumeEdgeLow = evaluation.front().scanAngle;
        plumeProperties.plumeEdgeHigh = evaluation.back().scanAngle;
        double peakLow = evaluation.front().scanAngle;
        double peakHigh = evaluation.back().scanAngle;
        const double maxCol = MaxColumn(evaluation);
        const double maxCol_div_e = maxCol * 0.3679;
        const double maxCol_90 = maxCol * 0.90;
        const double maxCol_half = maxCol * 0.5;

        // The left size of the plume
        for (int idx = 0; idx < numberOfGoodspectra - 1; ++idx)
        {
            if (evaluation[idx].scanAngle > plumeProperties.plumeCenter)
            {
                break;
            }
            if (evaluation[idx].offsetCorrectedColumn < maxCol_div_e)
            {
                plumeProperties.plumeEdgeLow = evaluation[idx].scanAngle;
            }
            if (evaluation[idx].offsetCorrectedColumn < maxCol_half)
            {
                plumeProperties.plumeHalfLow = evaluation[idx].scanAngle;
            }
            if ((evaluation[idx].offsetCorrectedColumn < maxCol_90) && (evaluation[idx + 1].offsetCorrectedColumn >= maxCol_90))
            {
                peakLow = evaluation[idx].scanAngle;
            }
        }

        // The right side of the plume
        for (int idx = numberOfGoodspectra - 1; idx > 0; --idx)
        {
            if (evaluation[idx].scanAngle <= plumeProperties.plumeCenter)
            {
                break;
            }
            if (evaluation[idx].offsetCorrectedColumn < maxCol_div_e)
            {
                plumeProperties.plumeEdgeHigh = evaluation[idx].scanAngle;
            }
            if (evaluation[idx].offsetCorrectedColumn < maxCol_half)
            {
                plumeProperties.plumeHalfHigh = evaluation[idx].scanAngle;
            }
            if ((evaluation[idx].offsetCorrectedColumn < maxCol_90) && (evaluation[idx - 1].offsetCorrectedColumn >= maxCol_90))
            {
                peakHigh = evaluation[idx].scanAngle;
            }
        }

        plumeProperties.plumeCenterError = (peakHigh - peakLow) / 2;

        return true;

    }

    bool FindPlume(const BasicScanEvaluationResult& evaluatedScan, int specieIdx, double plumeOffset, CPlumeInScanProperty& plumeProperties, std::string* message)
    {
        std::vector< ScanEvaluationData> evaluation;
        const long numPoints = static_cast<long>(evaluatedScan.m_spec.size());

        for (long idx = 0; idx < numPoints; ++idx)
        {
            if (evaluatedScan.m_spec[idx].IsBad())
            {
                continue;
            }

            ScanEvaluationData data;
            data.scanAngle = evaluatedScan.m_specInfo[idx].m_scanAngle;
            data.scanAngle2 = evaluatedScan.m_specInfo[idx].m_scanAngle2;
            data.offsetCorrectedColumn = evaluatedScan.m_spec[idx].m_referenceResult[specieIdx].m_column - plumeOffset;
            data.columnError = evaluatedScan.m_spec[idx].m_referenceResult[specieIdx].m_columnError;

            evaluation.push_back(data);
        }

        return FindPlume(evaluation, plumeProperties, message);
    }


    // VERSION 1: FROM NOVACPROGRAM
    bool FindPlume(
        const std::vector<double>& scanAngles,
        const std::vector<double>& phi,
        const std::vector<double>& columns,
        const std::vector<double>& columnErrors,
        const std::vector<bool>& badEvaluation,
        long numPoints,
        double plumeOffset,
        CPlumeInScanProperty& plumeProperties,
        std::string* message)
    {
        std::vector< ScanEvaluationData> evaluation;

        for (long idx = 0; idx < numPoints; ++idx)
        {
            if (badEvaluation[idx])
            {
                continue;
            }

            ScanEvaluationData data;
            data.scanAngle = scanAngles[idx];
            data.scanAngle2 = phi[idx];
            data.offsetCorrectedColumn = columns[idx] - plumeOffset;
            data.columnError = columnErrors[idx];

            evaluation.push_back(data);
        }

        return FindPlume(evaluation, plumeProperties, message);
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
        bool inPlume = ::novac::FindPlume(scanAngles, phi, columns, columnErrors, badEvaluation, numPoints, offset, plumeProperties, message);
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

        // Find the maximum offsetCorrectedColumn value
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
}
