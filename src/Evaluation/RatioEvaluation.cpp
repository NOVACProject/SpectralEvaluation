#include <SpectralEvaluation/Evaluation/RatioEvaluation.h>
#include <SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h>
#include <SpectralEvaluation/Evaluation/EvaluationResult.h>
#include <SpectralEvaluation/Evaluation/EvaluationBase.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>

namespace Evaluation
{
    RatioEvaluation::RatioEvaluation(const RatioEvaluationSettings& settings)
        : m_settings(settings)
    {
    }

    std::vector<Ratio> RatioEvaluation::Run(FileHandler::CScanFileHandler& scan)
    {
        std::vector<Ratio> result;
        
        if (!IsSuitableScanForRatioEvaluation(m_settings, m_masterResult, m_masterResultProperties))
        {
            return result;
        }

        std::vector<int> referenceSpectra;
        std::vector<int> plumeSpectra;
        SelectSpectraForRatioEvaluation(m_settings, m_masterResult, m_masterResultProperties, referenceSpectra, plumeSpectra);
        if (referenceSpectra.size() < (size_t)m_settings.minNumberOfReferenceSpectra ||
            plumeSpectra.size() < (size_t)m_settings.minNumberOfSpectraInPlume)
        {
            return result;
        }

        CSpectrum inPlumeSpectrum;
        AverageSpectra(scan, plumeSpectra, inPlumeSpectrum);

        CSpectrum outOfPlumeSpectrum;
        AverageSpectra(scan, referenceSpectra, outOfPlumeSpectrum);

        // TODO: Correct dark

        // Calculate the SO2 column in the spectrum
        double masterColumn = 0.0;
        {
            CEvaluationBase eval{ m_masterFitWindow };
            eval.SetSkySpectrum(outOfPlumeSpectrum);
            if (eval.Evaluate(inPlumeSpectrum))
            {
                // TODO: Handle errors here...
            }
            else
            {
                masterColumn = eval.m_result.m_referenceResult[0].m_column;
            }
        }

        for (const CFitWindow& window : m_referenceFit)
        {
            CEvaluationBase eval{window};
            eval.SetSkySpectrum(outOfPlumeSpectrum);
            if (eval.Evaluate(inPlumeSpectrum))
            {
                // Handle errors here
            }
            else
            {
                Ratio r;
                r.minorResult   = eval.m_result.m_referenceResult[0].m_column;
                r.majorResult   = masterColumn;
                r.ratio         = eval.m_result.m_referenceResult[0].m_column / masterColumn;
                // TODO: Error estimation
            }
        }

        return result;
    }

    bool IsSuitableScanForRatioEvaluation(const RatioEvaluationSettings& settings, const BasicScanEvaluationResult& scanResult, const CPlumeInScanProperty& properties)
    {
        if (scanResult.m_spec.size() < settings.minNumberOfSpectraInPlume + settings.minNumberOfReferenceSpectra)
        {
            return false; // not enough spectra
        }
        if (properties.completeness < 0.9)
        {
            return false;
        }
        return true; // TODO: Add more checks...
    }

    void SelectSpectraForRatioEvaluation(const RatioEvaluationSettings& settings, const BasicScanEvaluationResult& scanResult, const CPlumeInScanProperty& properties, std::vector<int>& referenceSpectra, std::vector<int>& inPlumeSpectra)
    {
        referenceSpectra.clear();
        inPlumeSpectra.clear();

        const int so2SpecieIndex = 0; // todo: figure this out

        for (size_t idx = 0; idx < scanResult.m_spec.size(); ++idx)
        {
            if (scanResult.m_specInfo[idx].m_scanAngle < properties.plumeEdgeLow ||
                scanResult.m_specInfo[idx].m_scanAngle > properties.plumeEdgeHigh)
            {
                referenceSpectra.push_back((int)idx);
            }
            else if(scanResult.m_spec[idx].m_referenceResult[so2SpecieIndex].m_column >= settings.minInPlumeColumn)
            {
                inPlumeSpectra.push_back((int)idx);
            }
        }

        return;
    }

    int AverageSpectra(FileHandler::CScanFileHandler& scan, const std::vector<int>& indices, CSpectrum& result)
    {
        if (indices.size() == 0)
        {
            return 0;
        }

        scan.GetSpectrum(result, indices[0]);
        int nofAveragedSpectra = 1;

        for (size_t ii = 1; ii < indices.size(); ++ii)
        {
            CSpectrum tmpSpec;
            if (scan.GetSpectrum(tmpSpec, indices[ii]))
            {
                result.Add(tmpSpec);
                ++nofAveragedSpectra;
            }
        }

        result.Div((double)nofAveragedSpectra);

        return nofAveragedSpectra;
    }
}