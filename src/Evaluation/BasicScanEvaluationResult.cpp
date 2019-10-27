#include <SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h>
#include <SpectralEvaluation/Utils.h>

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
