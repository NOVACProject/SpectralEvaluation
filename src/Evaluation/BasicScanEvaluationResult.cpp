#include <SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h>
#include <SpectralEvaluation/StringUtils.h>

using namespace Evaluation;

int BasicScanEvaluationResult::GetSpecieIndex(const char* specieName) const
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
    if (specieIndex < 0 || result.m_spec.size() == 0 || result.m_spec[0].m_referenceResult.size() <= specieIndex)
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