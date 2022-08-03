#include <SpectralEvaluation/Evaluation/RatioEvaluation.h>
#include <SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h>
#include <SpectralEvaluation/Evaluation/DarkSpectrum.h>
#include <SpectralEvaluation/Configuration/DarkSettings.h>
#include <SpectralEvaluation/Evaluation/EvaluationResult.h>
#include <SpectralEvaluation/Evaluation/BasicMath.h>
#include <SpectralEvaluation/Evaluation/DoasFit.h>
#include <SpectralEvaluation/Evaluation/DoasFitPreparation.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/StringUtils.h>
#include <SpectralEvaluation/Math/SpectrumMath.h>

#include <numeric>

namespace novac
{

    bool IsSuitableScanForRatioEvaluation(const RatioEvaluationSettings& settings, const BasicScanEvaluationResult& scanResult, const CPlumeInScanProperty& properties)
    {
        if (static_cast<int>(scanResult.m_spec.size()) < settings.minNumberOfSpectraInPlume + settings.minNumberOfReferenceSpectra)
        {
            return false; // not enough spectra
        }
        if (properties.completeness < 0.7)
        {
            return false;
        }
        return true; // TODO: Add more checks...
    }

    int AddAsSky(const std::vector<double>& referenceData, CFitWindow& window, SHIFT_TYPE shiftOption = SHIFT_TYPE::SHIFT_FIX)
    {
        int indexOfSkySpectrum = window.nRef;

        window.ref[window.nRef].m_data = std::make_unique<novac::CCrossSectionData>(referenceData);
        window.ref[window.nRef].m_specieName = "sky";
        window.ref[window.nRef].m_columnOption = novac::SHIFT_TYPE::SHIFT_FIX;
        window.ref[window.nRef].m_columnValue = -1.0;
        window.ref[window.nRef].m_squeezeOption = novac::SHIFT_TYPE::SHIFT_FIX;
        window.ref[window.nRef].m_squeezeValue = 1.0;
        window.ref[window.nRef].m_shiftOption = shiftOption;
        window.ref[window.nRef].m_shiftValue = 0.0;
        window.nRef += 1;

        return indexOfSkySpectrum;
    }

    void AddAsReference(const std::vector<double>& referenceData, CFitWindow& window, const std::string& name, int linkShiftToIdx = -1)
    {
        window.ref[window.nRef].m_data = std::make_unique<novac::CCrossSectionData>(referenceData);
        window.ref[window.nRef].m_specieName = name;
        window.ref[window.nRef].m_columnOption = novac::SHIFT_TYPE::SHIFT_FREE;
        window.ref[window.nRef].m_columnValue = 1.0;
        window.ref[window.nRef].m_squeezeOption = novac::SHIFT_TYPE::SHIFT_FIX;
        window.ref[window.nRef].m_squeezeValue = 1.0;

        if (linkShiftToIdx >= 0)
        {
            window.ref[window.nRef].m_shiftOption = SHIFT_TYPE::SHIFT_LINK;
            window.ref[window.nRef].m_shiftValue = linkShiftToIdx;
        }
        else
        {
            window.ref[window.nRef].m_shiftOption = SHIFT_TYPE::SHIFT_FIX;
            window.ref[window.nRef].m_shiftValue = 0.0;
        }

        window.nRef += 1;
    }

    int FindReferenceIndex(const CFitWindow& window, const std::string& nameToFind)
    {
        for (int ii = 0; ii < window.nRef; ++ii)
        {
            if (EqualsIgnoringCase(window.ref[ii].m_specieName, nameToFind))
            {
                return ii;
            }
        }

        return -1;
    }

    std::vector<double> PrepareRingLambda4Spectrum(const std::vector<double>& ringSpectrum, const std::vector<double>& wavelength)
    {
        if (ringSpectrum.size() != wavelength.size())
        {
            throw std::invalid_argument("Cannot calculate a ringxlambda4 spectrum if the length of the ring spectrum differs from the wavelength calibration length.");
        }

        std::vector<double> result;
        result.resize(ringSpectrum.size());

        for (size_t ii = 0; ii < ringSpectrum.size(); ++ii)
        {
            const double lambda = wavelength[ii];
            result[ii] = ringSpectrum[ii] * std::pow(lambda, 4.0);
        }

        return result;
    }

    RatioEvaluation::RatioEvaluation(const RatioEvaluationSettings& settings, const Configuration::CDarkSettings& darkSettings)
        : m_darkSettings(darkSettings), m_settings(settings)
    {
    }

    std::vector<Ratio> RatioEvaluation::Run(IScanSpectrumSource& scan, std::string* errorMessage)
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
            if (errorMessage != nullptr)
            {
                *errorMessage = "Too few suitable spectra for in-plume or out-of-plume";
            }

            return result;
        }

        // Create the in-plume spectrum
        CSpectrum inPlumeSpectrum;
        AverageSpectra(scan, plumeSpectra, inPlumeSpectrum);

        CSpectrum inPlumeDark;
        if (!::novac::GetDark(scan, inPlumeSpectrum, m_darkSettings, inPlumeDark, m_lastErrorMessage))
        {
            if (errorMessage != nullptr)
            {
                *errorMessage = m_lastErrorMessage;
            }
            return result;
        }
        if (inPlumeDark.NumSpectra() > 0 && !m_averagedSpectra)
        {
            inPlumeDark.Div(inPlumeDark.NumSpectra());
        }
        inPlumeSpectrum.Sub(inPlumeDark);


        // Create the out-of-plume spectrum
        CSpectrum outOfPlumeSpectrum;
        AverageSpectra(scan, referenceSpectra, outOfPlumeSpectrum);

        CSpectrum outOfPlumeDark;
        if (!::novac::GetDark(scan, outOfPlumeSpectrum, m_darkSettings, outOfPlumeDark, m_lastErrorMessage))
        {
            if (errorMessage != nullptr)
            {
                *errorMessage = m_lastErrorMessage;
            }
            return result;
        }
        if (outOfPlumeDark.NumSpectra() > 0 && !m_averagedSpectra)
        {
            outOfPlumeDark.Div(outOfPlumeDark.NumSpectra());
        }
        outOfPlumeSpectrum.Sub(outOfPlumeDark);

        // Get the pixel-to-wavelength calibration of the out-of-plume-spectrum.
        // TODO: Decide on a better input for this.
        if (m_masterFitWindow.fraunhoferRef.m_data != nullptr && m_masterFitWindow.fraunhoferRef.m_data->m_waveLength.size() == outOfPlumeSpectrum.m_length)
        {
            outOfPlumeSpectrum.m_wavelength = m_masterFitWindow.fraunhoferRef.m_data->m_waveLength;
        }

        /* Notice that there are several improvements which can be done to the fit here, some TODO:s
        * TODO: If the references are calibrated towards the sky spectrum (could be an input option here) then we can link the shift of all the references to sky.
        * TODO: The shift of the sky can be linked to the calculated Ring and RingxLambda4
        */

        // Create the DOAS setup for the master (SO2) evaluation
        CFitWindow localSO2FitWindow = m_masterFitWindow;

        const auto filteredOutOfPlumeSpectrum = DoasFitPreparation::PrepareSkySpectrum(outOfPlumeSpectrum, FIT_TYPE::FIT_POLY);
        const int skySpectrumIdx = AddAsSky(filteredOutOfPlumeSpectrum, localSO2FitWindow, SHIFT_TYPE::SHIFT_FREE);

        const auto ringSpectrum = DoasFitPreparation::PrepareRingSpectrum(outOfPlumeSpectrum, FIT_TYPE::FIT_POLY);
        AddAsReference(ringSpectrum, localSO2FitWindow, "ring", skySpectrumIdx);

        const auto ringLambda4Spectrum = PrepareRingLambda4Spectrum(ringSpectrum, outOfPlumeSpectrum.m_wavelength);
        AddAsReference(ringLambda4Spectrum, localSO2FitWindow, "ringLambda4", skySpectrumIdx);

        const auto offsetPoly = DoasFitPreparation::PrepareIntensitySpacePolynomial(outOfPlumeSpectrum);
        AddAsReference(offsetPoly, localSO2FitWindow, "offset", skySpectrumIdx);

        DoasFit doas;
        doas.Setup(localSO2FitWindow);

        const auto filteredMeasuredSpectrum = DoasFitPreparation::PrepareMeasuredSpectrum(inPlumeSpectrum, outOfPlumeSpectrum, FIT_TYPE::FIT_POLY);

        // Calculate the SO2 column in the spectrum
        DoasResult so2DoasResult;
        try
        {
            doas.Run(filteredMeasuredSpectrum.data(), filteredMeasuredSpectrum.size(), so2DoasResult);
        }
        catch (DoasFitException&)
        {
            // TODO: Handle this error properly.
            return result;
        }


        // Now finally calculate the ratios
        for (const CFitWindow& window : m_referenceFit)
        {
            // Create the DOAS setup for the master evaluation
            CFitWindow broFitWindow = window;

            // Fix SO2 to the column retrieved in the first result
            const int so2RefIdx = FindReferenceIndex(broFitWindow, so2DoasResult.referenceResult[0].name);
            if (so2RefIdx < 0)
            {
                if (errorMessage != nullptr)
                {
                    *errorMessage = "Failed to find the SO2 reference in BrO fit window";
                }
                break;
            }
            broFitWindow.ref[so2RefIdx].m_columnOption = SHIFT_TYPE::SHIFT_FIX;
            broFitWindow.ref[so2RefIdx].m_columnValue = so2DoasResult.referenceResult[0].column;

            // Add the out-of-plume spectrum as a reference
            const int indexOfSkySpectrumInBrO = AddAsSky(filteredOutOfPlumeSpectrum, broFitWindow, SHIFT_TYPE::SHIFT_FREE);

            AddAsReference(ringSpectrum, broFitWindow, "ring", indexOfSkySpectrumInBrO);
            AddAsReference(ringLambda4Spectrum, localSO2FitWindow, "ringLambda4", skySpectrumIdx);
            AddAsReference(offsetPoly, broFitWindow, "offset", skySpectrumIdx);

            DoasFit broDoasSetup;
            broDoasSetup.Setup(broFitWindow);

            DoasResult broDoasResult;
            try
            {
                broDoasSetup.Run(filteredMeasuredSpectrum.data(), filteredMeasuredSpectrum.size(), broDoasResult);
            }
            catch (DoasFitException&)
            {
                // TODO: Handle this error properly.
                continue;
            }

            // we can now calculate a ratio
            Ratio r;
            r.minorResult = broDoasResult.referenceResult[0].column;
            r.minorError = broDoasResult.referenceResult[0].columnError;
            r.minorSpecieName = broDoasResult.referenceResult[0].name;

            r.majorResult = so2DoasResult.referenceResult[0].column;
            r.majorError = so2DoasResult.referenceResult[0].columnError;
            r.majorSpecieName = so2DoasResult.referenceResult[0].name;

            r.ratio = r.minorResult / r.majorResult;
            r.error = std::abs(r.ratio) * std::sqrt(std::pow(r.minorError / r.minorResult, 2.0) + std::pow(r.majorError / r.majorResult, 2.0));

            result.push_back(r);
        }
        return result;
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
}