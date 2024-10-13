#include <SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h>
#include <SpectralEvaluation/StringUtils.h>

namespace novac
{

int BasicScanEvaluationResult::AppendResult(const CEvaluationResult& evalRes, const CSpectrumInfo& specInfo)
{
    // Append the evaluationresult to the end of the 'm_spec'-vector
    m_spec.push_back(CEvaluationResult(evalRes));

    // Append the spectral information to the end of the 'm_specInfo'-vector
    m_specInfo.push_back(CSpectrumInfo(specInfo));

    // Increase the numbers of spectra in this result-set.
    ++m_specNum;
    return 0;
}

int BasicScanEvaluationResult::RemoveResult(unsigned int specIndex)
{
    if (specIndex >= m_specNum)
    {
        return 1; // not a valid index
    }

    // Remove the desired value
    m_specInfo.erase(begin(m_specInfo) + specIndex);

    // Decrease the number of values in the list
    m_specNum -= 1;

    return 0;
}

void BasicScanEvaluationResult::InitializeArrays(long size)
{
    if (size < 0 || size > 1024)
    {
        return;
    }

    m_spec.reserve(size);
    m_specInfo.reserve(size);
}

int BasicScanEvaluationResult::GetSpecieIndex(const std::string& specieName) const
{
    // If there are no spectra, there can be no species
    if (m_spec.size() <= 0)
    {
        return -1;
    }

    // if there's only one specie, assume that this is the correct one
    if (m_spec[0].m_referenceResult.size() == 1)
    {
        return 0;
    }

    for (size_t i = 0; i < m_spec[0].m_referenceResult.size(); ++i)
    {
        if (EqualsIgnoringCase(m_spec[0].m_referenceResult[i].m_specieName, specieName))
        {
            return (int)i;
        }
    }

    return -1;
}

std::vector<double> GetColumns(const BasicScanEvaluationResult& result, int specieIndex)
{
    std::vector<double> columns;
    if (specieIndex < 0 || result.m_spec.size() == 0 || result.m_spec[0].m_referenceResult.size() <= (size_t)specieIndex)
    {
        return columns;
    }

    columns.reserve(result.m_spec.size());

    for (size_t ii = 0; ii < result.m_spec.size(); ++ii)
    {
        columns.push_back(result.m_spec[ii].m_referenceResult[specieIndex].m_column);
    }

    return columns;
}
}
