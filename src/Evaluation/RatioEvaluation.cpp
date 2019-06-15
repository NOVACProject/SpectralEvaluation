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
        for (int specIdx : plumeSpectra)
        {
            CSpectrum tmpSpec;
            if (scan.GetSpectrum(tmpSpec, specIdx))
            {
                inPlumeSpectrum.Add(tmpSpec);
            }
        }

        CSpectrum outOfPlumeSpectrum;
        for (int specIdx : referenceSpectra)
        {
            CSpectrum tmpSpec;
            if (scan.GetSpectrum(tmpSpec, specIdx))
            {
                outOfPlumeSpectrum.Add(tmpSpec);
            }
        }


        for (const CFitWindow& window : m_referenceFit)
        {
            CEvaluationBase eval{window};
            eval.SetSkySpectrum(outOfPlumeSpectrum);
            eval.Evaluate(inPlumeSpectrum);
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

}