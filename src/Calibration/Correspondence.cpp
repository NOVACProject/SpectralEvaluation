#include <SpectralEvaluation/Calibration/Correspondence.h>
#include <SpectralEvaluation/VectorUtils.h>

#include <stdexcept>

namespace novac
{

// TODO: Unit-test these mehtods somewhere
bool AllMeasuredPointsAreUnique(const std::vector<Correspondence>& correspondences)
{
    if (correspondences.size() == 0)
    {
        return true;
    }

    std::vector<size_t> selecedMeasuredPoints;
    selecedMeasuredPoints.push_back(correspondences[0].measuredIdx);
    for (size_t correspondenceIdx = 1; correspondenceIdx < correspondences.size(); ++correspondenceIdx)
    {
        if (Contains(selecedMeasuredPoints, correspondences[correspondenceIdx].measuredIdx))
        {
            return false;
        }
        selecedMeasuredPoints.push_back(correspondences[correspondenceIdx].measuredIdx);
    }

    return true;
}

bool AllTheoreticalPointsAreUnique(const std::vector<Correspondence>& correspondences)
{
    if (correspondences.size() == 0)
    {
        return true;
    }

    std::vector<size_t> selecedTheoreticalPoints;
    selecedTheoreticalPoints.push_back(correspondences[0].theoreticalIdx);
    for (size_t correspondenceIdx = 1; correspondenceIdx < correspondences.size(); ++correspondenceIdx)
    {
        if (Contains(selecedTheoreticalPoints, correspondences[correspondenceIdx].theoreticalIdx))
        {
            return false;
        }
        selecedTheoreticalPoints.push_back(correspondences[correspondenceIdx].theoreticalIdx);
    }

    return true;
}

void ListAllPossibleKeypointCombinationsStartingWith(
    size_t number,
    const std::vector<Correspondence>& allCorrespondences,
    std::vector<size_t>& selectedIndices,
    std::vector<std::vector<Correspondence>>& result)
{
    if (selectedIndices.size() == number)
    {
        // We're done searching for indices. Check that this is a unique and valid combination of correspondences.
        std::vector<Correspondence> selectedCorrespondences;
        for (size_t idx : selectedIndices)
        {
            selectedCorrespondences.push_back(allCorrespondences[idx]);
        }
        if (AllMeasuredPointsAreUnique(selectedCorrespondences) && AllTheoreticalPointsAreUnique(selectedCorrespondences))
        {
            result.push_back(selectedCorrespondences);
        }
        return;
    }

    for (size_t firstMeasuredIdx = selectedIndices.back() + 1; firstMeasuredIdx < allCorrespondences.size(); ++firstMeasuredIdx)
    {
        selectedIndices.push_back(firstMeasuredIdx);
        ListAllPossibleKeypointCombinationsStartingWith(number, allCorrespondences, selectedIndices, result);
        selectedIndices.pop_back();
    }

}

std::vector<std::vector<Correspondence>> ListAllPossibleKeypointCombinations(
    size_t number,
    const std::vector<Correspondence>& allCorrespondences)
{
    std::vector<std::vector<Correspondence>> result;

    for (size_t firstMeasuredIdx = 0; firstMeasuredIdx < allCorrespondences.size(); ++firstMeasuredIdx)
    {
        std::vector<size_t> selectedIndices;
        selectedIndices.push_back(firstMeasuredIdx);

        ListAllPossibleKeypointCombinationsStartingWith(number, allCorrespondences, selectedIndices, result);
    }

    return result;
}

void SelectMaybeInliers(size_t number, const std::vector<Correspondence>& allCorrespondences, std::mt19937& randomGenerator, std::vector<Correspondence>& result)
{
    if (number == allCorrespondences.size())
    {
        result = allCorrespondences;
        return;
    }
    else if (number > allCorrespondences.size())
    {
        throw std::invalid_argument("Cannot select a larger number of inliers than the entire data set.");
    }

    result.resize(number);
    std::uniform_int_distribution<size_t> uniform_dist(0, allCorrespondences.size() - 1);

    // "Randomly" select correspondences with the restriction that:
    //  a given measured or theoretical point must not be present twice in the selected result
    std::vector<size_t> selectedIndices;
    selectedIndices.reserve(number);
    std::vector<size_t> selecedMeasuredPoints;
    selecedMeasuredPoints.reserve(number);
    std::vector<size_t> selecedTheoreticalPoints;
    selecedTheoreticalPoints.reserve(number);

    while (selectedIndices.size() < number)
    {
        const size_t suggestedIdx = uniform_dist(randomGenerator);

        if (!Contains(selecedMeasuredPoints, allCorrespondences[suggestedIdx].measuredIdx) &&
            !Contains(selecedTheoreticalPoints, allCorrespondences[suggestedIdx].theoreticalIdx))
        {
            selecedMeasuredPoints.push_back(allCorrespondences[suggestedIdx].measuredIdx);
            selecedTheoreticalPoints.push_back(allCorrespondences[suggestedIdx].theoreticalIdx);
            selectedIndices.push_back(suggestedIdx);
        }
    }

    size_t resultIdx = 0;
    for (size_t idx : selectedIndices)
    {
        result[resultIdx] = allCorrespondences[idx];
        ++resultIdx;
    }

    return;
}

int FindCorrespondenceWithTheoreticalIdx(const std::vector<Correspondence>& allEntries, size_t theoreticalIdxOfEntryToFind)
{
    for (int ii = 0; ii < static_cast<int>(allEntries.size()); ++ii)
    {
        if (allEntries[ii].theoreticalIdx == theoreticalIdxOfEntryToFind)
        {
            return ii;
        }
    }

    return -1; // not found
}

double MeasurePixelSpan(const std::vector<Correspondence>& correspondences)
{
    if (correspondences.size() == 0)
    {
        return 0.0;
    }
    double minPixel = correspondences[0].measuredValue;
    double maxPixel = correspondences[0].measuredValue;

    for (size_t ii = 1; ii < correspondences.size(); ++ii)
    {
        minPixel = std::min(minPixel, correspondences[ii].measuredValue);
        maxPixel = std::max(maxPixel, correspondences[ii].measuredValue);
    }

    return (maxPixel - minPixel);
}

std::vector<std::vector<Correspondence>> ArrangeByMeasuredKeypoint(const std::vector<Correspondence>& allCorrespondences)
{
    std::vector<std::vector<Correspondence>> result;
    result.reserve(100); // initial guess

    for (const Correspondence& c : allCorrespondences)
    {
        if (c.measuredIdx >= result.size())
        {
            result.resize(1 + c.measuredIdx);
        }
        result[c.measuredIdx].push_back(c);
    }

    return result;
}

std::vector<bool> ListInliers(const std::vector<Correspondence>& selectedValues, const std::vector<Correspondence>& allValues)
{
    std::vector<bool> result(allValues.size(), false);

    // TODO: This is not efficient!
    for (const auto& selectedValue : selectedValues)
    {
        for (size_t ii = 0; ii < allValues.size(); ++ii)
        {
            if (allValues[ii].measuredIdx == selectedValue.measuredIdx && allValues[ii].theoreticalIdx == selectedValue.theoreticalIdx)
            {
                result[ii] = true;
                break;
            }
        }
    }

    return result;
}


}
