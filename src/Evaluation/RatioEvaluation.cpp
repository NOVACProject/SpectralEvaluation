#include <SpectralEvaluation/Evaluation/RatioEvaluation.h>
#include <SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h>
#include <SpectralEvaluation/Evaluation/DarkSpectrum.h>
#include <SpectralEvaluation/Configuration/DarkSettings.h>
#include <SpectralEvaluation/Evaluation/EvaluationResult.h>
#include <SpectralEvaluation/Evaluation/EvaluationBase.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>

#include <numeric>

namespace Evaluation
{
    RatioEvaluation::RatioEvaluation(const RatioEvaluationSettings& settings, const Configuration::CDarkSettings& darkSettings)
        : m_darkSettings(darkSettings), m_settings(settings)
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

        CSpectrum inPlumeDark;
        if (!Evaluation::GetDark(scan, inPlumeSpectrum, m_darkSettings, inPlumeDark, m_lastErrorMessage))
        {
            // TODO: Handle errors
        }
        if (inPlumeDark.NumSpectra() > 0 && !m_averagedSpectra)
        {
            inPlumeDark.Div(inPlumeDark.NumSpectra());
        }
        inPlumeSpectrum.Sub(inPlumeDark);



        CSpectrum outOfPlumeSpectrum;
        AverageSpectra(scan, referenceSpectra, outOfPlumeSpectrum);

        CSpectrum outOfPlumeDark;
        if (!Evaluation::GetDark(scan, outOfPlumeSpectrum, m_darkSettings, outOfPlumeDark, m_lastErrorMessage))
        {
            // TODO: Handle errors
        }
        if (outOfPlumeDark.NumSpectra() > 0 && !m_averagedSpectra)
        {
            outOfPlumeDark.Div(outOfPlumeDark.NumSpectra());
        }
        outOfPlumeSpectrum.Sub(outOfPlumeDark);



        // Calculate the SO2 column in the spectrum
        double masterColumn = 0.0;
        double masterColumnError = 0.0;
        std::string masterSpecieName;
        {
            CEvaluationBase eval{ m_masterFitWindow };
            eval.SetSkySpectrum(outOfPlumeSpectrum);
            if (eval.Evaluate(inPlumeSpectrum))
            {
                // TODO: Handle errors here...
            }
            else
            {
                masterColumn        = eval.m_result.m_referenceResult[0].m_column;
                masterColumnError   = eval.m_result.m_referenceResult[0].m_columnError;
                masterSpecieName    = m_masterFitWindow.ref[0].m_specieName;
            }
        }


        // Now finally calculate the ratios
        for (const CFitWindow& window : m_referenceFit)
        {
            CFitWindow windowCopy = window;
            ReadReferences(windowCopy);
            CEvaluationBase eval{ windowCopy };
            eval.SetSkySpectrum(outOfPlumeSpectrum);
            if (eval.Evaluate(inPlumeSpectrum))
            {
                // Handle errors here
            }
            else
            {
                Ratio r;
                r.minorResult       = eval.m_result.m_referenceResult[0].m_column;
                r.minorError        = eval.m_result.m_referenceResult[0].m_columnError;
                r.minorSpecieName   = window.ref[0].m_specieName;

                r.majorResult       = masterColumn;
                r.majorError        = masterColumnError;
                r.majorSpecieName   = masterSpecieName;

                r.ratio             = eval.m_result.m_referenceResult[0].m_column / masterColumn;
                r.error             = std::sqrt( std::pow(r.minorError / r.minorResult, 2.0) + std::pow(r.majorError / r.majorResult, 2.0) );

                result.push_back(r);
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
        if (properties.completeness < 0.7)
        {
            return false;
        }
        return true; // TODO: Add more checks...
    }

    // Calculates the average column value of the given specie in the index [startIdx, endIdx]
    double AverageColumnValue(const BasicScanEvaluationResult& scanResult, int specieIndex, size_t startIdx, size_t endIdx)
    {
        double sum = 0.0;
        for (size_t ii = startIdx; ii < endIdx; ++ii)
        {
            sum += scanResult.m_spec[ii].m_referenceResult[specieIndex].m_column;
        }
        return sum / (double)(endIdx - startIdx);
    }

    void SelectSpectraForRatioEvaluation(const RatioEvaluationSettings& settings, const BasicScanEvaluationResult& scanResult, const CPlumeInScanProperty& properties, std::vector<int>& referenceSpectra, std::vector<int>& inPlumeSpectra)
    {
        referenceSpectra.clear();
        inPlumeSpectra.clear();

        const size_t referenceRegionWidth = 10U;

        if (scanResult.m_spec.size() <= referenceRegionWidth)
        {
            return; // Cannot retrieve a region, too few spectra...
        }

        const int so2SpecieIndex = 0; // TODO: figure this out

        // Step 1, find the in-plume region.
        for (size_t idx = 0; idx < scanResult.m_spec.size(); ++idx)
        {
            if (scanResult.m_specInfo[idx].m_scanAngle >= properties.plumeEdgeLow &&
                scanResult.m_specInfo[idx].m_scanAngle <= properties.plumeEdgeHigh &&
                scanResult.m_spec[idx].m_referenceResult[so2SpecieIndex].m_column >= settings.minInPlumeColumn)
            {
                inPlumeSpectra.push_back((int)idx);
            }
        }

        // Step 2, find the reference region (10 adjacent spectra with lowest avg column value)
        size_t startIdxOfReferenceRegion = 0U;
        double lowestMeanColumnValue = std::numeric_limits<double>::max();
        for (size_t startIdxCandidate = 0U; startIdxCandidate < scanResult.m_spec.size() - referenceRegionWidth; ++startIdxCandidate)
        {
            const double meanColumnValue = AverageColumnValue(scanResult, so2SpecieIndex, startIdxCandidate, startIdxCandidate + referenceRegionWidth);
            if (meanColumnValue < lowestMeanColumnValue)
            {
                lowestMeanColumnValue = meanColumnValue;
                startIdxOfReferenceRegion = startIdxCandidate;
            }
        }

        referenceSpectra = std::vector<int>(referenceRegionWidth);
        std::iota(begin(referenceSpectra), end(referenceSpectra), (int)startIdxOfReferenceRegion);

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