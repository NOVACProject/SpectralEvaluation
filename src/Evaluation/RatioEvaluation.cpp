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
#include <sstream>

namespace novac
{
    struct PreparedInputSpectraForDoasEvaluation
    {
        std::vector<double> filteredInPlumeSpectrum;

        std::vector<double> filteredOutOfPlumespectrum;

        std::vector<double> ringSpectrum;

        std::vector<double> ringLambda4Spectrum;

        std::vector<double> intensityOffsetSpectrum;
    };

    bool IsSuitableScanForRatioEvaluation(const RatioEvaluationSettings& settings, const BasicScanEvaluationResult& scanResult, const CPlumeInScanProperty& properties, std::string& errorMessage)
    {
        if (static_cast<int>(scanResult.m_spec.size()) < settings.minNumberOfSpectraInPlume + settings.minNumberOfReferenceSpectra)
        {
            errorMessage = "Too few spectra in scan for creating adding an in-plume-spectrum and an out-of-plume-spectrum.";
            return false; // not enough spectra
        }
        if (properties.completeness < 0.7)
        {
            errorMessage = "Plume completeness below threshold of 0.7";
            return false;
        }
        return true; // TODO: Add more checks...
    }

    void RatioEvaluation::DarkCorrectSpectrum(IScanSpectrumSource& scan, CSpectrum& spectrum) const
    {
        std::string errorMessage;
        auto correspondingDarkSpectrum = std::make_unique<CSpectrum>();
        if (!::novac::GetDark(scan, spectrum, m_darkSettings, *correspondingDarkSpectrum, errorMessage))
        {
            return throw std::invalid_argument(errorMessage);
        }
        if (correspondingDarkSpectrum->NumSpectra() > 0 && !m_averagedSpectra)
        {
            correspondingDarkSpectrum->Div(correspondingDarkSpectrum->NumSpectra());
        }
        spectrum.Sub(*correspondingDarkSpectrum);
    }

    void AddRingSpectraAsReferences(CFitWindow& localSO2FitWindow, const std::vector<double>& ringSpectrum, const std::vector<double>& ringLambda4Spectrum, int skySpectrumIdx)
    {
        if (localSO2FitWindow.ringCalculation != RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING)
        {
            AddAsReference(localSO2FitWindow, ringSpectrum, "ring", skySpectrumIdx);

            if (localSO2FitWindow.ringCalculation == RING_CALCULATION_OPTION::CALCULATE_RING_X2)
            {
                AddAsReference(localSO2FitWindow, ringLambda4Spectrum, "ringLambda4", skySpectrumIdx);
            }
        }
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

    Ratio RatioEvaluation::CalculateRatio(const DoasResult& majorWindow, const DoasResult& minorWindow)
    {
        Ratio ratio;
        ratio.minorResult = minorWindow.referenceResult[0].column;
        ratio.minorError = minorWindow.referenceResult[0].columnError;
        ratio.minorSpecieName = minorWindow.referenceResult[0].name;

        ratio.majorResult = majorWindow.referenceResult[0].column;
        ratio.majorError = majorWindow.referenceResult[0].columnError;
        ratio.majorSpecieName = majorWindow.referenceResult[0].name;

        ratio.ratio = ratio.minorResult / ratio.majorResult;
        ratio.error = std::abs(ratio.ratio) * std::sqrt(std::pow(ratio.minorError / ratio.minorResult, 2.0) + std::pow(ratio.majorError / ratio.majorResult, 2.0));

        return ratio;
    }

    DoasResult RunEvaluation(const CFitWindow& window, PreparedInputSpectraForDoasEvaluation& spectra)
    {
        CFitWindow localCopyOfWindow = window;

        int skySpectrumIdx = 0;
        if (window.fitType != FIT_TYPE::FIT_HP_DIV)
        {
            skySpectrumIdx = AddAsSky(localCopyOfWindow, spectra.filteredOutOfPlumespectrum, SHIFT_TYPE::SHIFT_FIX);
        }

        AddRingSpectraAsReferences(localCopyOfWindow, spectra.ringSpectrum, spectra.ringLambda4Spectrum, skySpectrumIdx);

        if (localCopyOfWindow.includeIntensitySpacePolyominal)
        {
            AddAsReference(localCopyOfWindow, spectra.intensityOffsetSpectrum, "offset", skySpectrumIdx);
        }

        DoasFit doas;
        doas.Setup(localCopyOfWindow);

        DoasResult doasResult;
        doas.Run(spectra.filteredInPlumeSpectrum.data(), spectra.filteredInPlumeSpectrum.size(), doasResult);

        return doasResult;
    }

    RatioEvaluation::RatioEvaluation(const RatioEvaluationSettings& settings, const Configuration::CDarkSettings& darkSettings)
        : m_darkSettings(darkSettings), m_settings(settings)
    {
    }

    std::vector<Ratio> RatioEvaluation::Run(IScanSpectrumSource& scan, RatioEvaluationDebugInformation* debugInfo)
    {
        if (debugInfo != nullptr)
        {
            return this->Run(scan, *debugInfo);
        }
        else
        {
            RatioEvaluationDebugInformation localDebugInfo;
            return this->Run(scan, localDebugInfo);
        }
    }

    std::vector<Ratio> RatioEvaluation::Run(IScanSpectrumSource& scan, RatioEvaluationDebugInformation& debugInfo)
    {
        std::vector<Ratio> result;
        debugInfo.doasResults.clear();

        try
        {
            ValidateSetup();

            if (!IsSuitableScanForRatioEvaluation(m_settings, m_masterResult, m_masterResultProperties, debugInfo.errorMessage))
            {
                return result;
            }

            SelectSpectraForRatioEvaluation(m_settings, m_masterResult, m_masterResultProperties, debugInfo.outOfPlumeSpectra, debugInfo.plumeSpectra);

            if (debugInfo.outOfPlumeSpectra.size() < (size_t)m_settings.minNumberOfReferenceSpectra ||
                debugInfo.plumeSpectra.size() < (size_t)m_settings.minNumberOfSpectraInPlume)
            {
                std::stringstream message;
                message << "Too few suitable spectra for in-plume or out-of-plume. ";
                message << "In plume: " << debugInfo.plumeSpectra.size() << " selected and " << m_settings.minNumberOfSpectraInPlume << " required. ";
                message << "Out of plume: " << debugInfo.outOfPlumeSpectra.size() << " selected and " << m_settings.minNumberOfReferenceSpectra << " required. ";
                debugInfo.errorMessage = message.str();

                return result;
            }

            // Create the in-plume spectrum
            AverageSpectra(scan, debugInfo.plumeSpectra, debugInfo.inPlumeSpectrum);
            DarkCorrectSpectrum(scan, debugInfo.inPlumeSpectrum);

            // Create the out-of-plume spectrum
            AverageSpectra(scan, debugInfo.outOfPlumeSpectra, debugInfo.outOfPlumeSpectrum);
            DarkCorrectSpectrum(scan, debugInfo.outOfPlumeSpectrum);

            // Get the pixel-to-wavelength calibration of the out-of-plume-spectrum.
            // TODO: Decide on a better input for this.
            if (m_masterFitWindow.fraunhoferRef.m_data != nullptr && m_masterFitWindow.fraunhoferRef.m_data->m_waveLength.size() == debugInfo.outOfPlumeSpectrum.m_length)
            {
                debugInfo.outOfPlumeSpectrum.m_wavelength = m_masterFitWindow.fraunhoferRef.m_data->m_waveLength;
            }
            else if (m_masterFitWindow.ref[0].m_data->m_waveLength.size() == m_masterFitWindow.ref[0].m_data->m_crossSection.size())
            {
                debugInfo.outOfPlumeSpectrum.m_wavelength = m_masterFitWindow.ref[0].m_data->m_waveLength;
            }
            else
            {
                debugInfo.errorMessage = "Failed to retrieve a pixel-to-wavelength calibration for the setup";
                return result;
            }

            /* Notice that there are several improvements which can be done to the fit here, some TODO:s
            * TODO: If the references are calibrated towards the sky spectrum (could be an input option here) then we can link the shift of all the references to sky.
            * TODO: The shift of the sky can be linked to the calculated Ring and RingxLambda4
            */
            PreparedInputSpectraForDoasEvaluation filteredSpectra;

            if (m_masterFitWindow.fitType != FIT_TYPE::FIT_HP_DIV)
            {
                filteredSpectra.filteredOutOfPlumespectrum = DoasFitPreparation::PrepareSkySpectrum(debugInfo.outOfPlumeSpectrum, m_masterFitWindow.fitType);
            }
            else
            {
                DoasFitPreparation::RemoveOffset(debugInfo.outOfPlumeSpectrum);
            }
            filteredSpectra.filteredInPlumeSpectrum = DoasFitPreparation::PrepareMeasuredSpectrum(debugInfo.inPlumeSpectrum, debugInfo.outOfPlumeSpectrum, m_masterFitWindow.fitType);
            filteredSpectra.intensityOffsetSpectrum = DoasFitPreparation::PrepareIntensitySpacePolynomial(debugInfo.outOfPlumeSpectrum);

            if (AnyFitWindowRequiresRingSpectrum())
            {
                filteredSpectra.ringSpectrum = DoasFitPreparation::PrepareRingSpectrum(debugInfo.outOfPlumeSpectrum, m_masterFitWindow.fitType);
                filteredSpectra.ringLambda4Spectrum = PrepareRingLambda4Spectrum(filteredSpectra.ringSpectrum, debugInfo.outOfPlumeSpectrum.m_wavelength);
            }

            // Calculate the SO2 column in the spectrum
            DoasResult so2DoasResult = RunEvaluation(m_masterFitWindow, filteredSpectra);
            debugInfo.doasResults.push_back(so2DoasResult);

            // Evaluate for the minor species and calculate the ratios
            for (const CFitWindow& window : m_referenceFit)
            {
                CFitWindow broFitWindow = window;

                // Fix SO2 to the column retrieved in the first result (if it is included in this window)
                const int so2RefIdx = FindReferenceIndex(broFitWindow, so2DoasResult.referenceResult[0].name);
                if (so2RefIdx > 0)
                {
                    broFitWindow.ref[so2RefIdx].m_columnOption = SHIFT_TYPE::SHIFT_FIX;
                    broFitWindow.ref[so2RefIdx].m_columnValue = so2DoasResult.referenceResult[0].column;
                }

                DoasResult broDoasResult = RunEvaluation(broFitWindow, filteredSpectra);
                debugInfo.doasResults.push_back(broDoasResult);

                Ratio r = CalculateRatio(so2DoasResult, broDoasResult);

                result.push_back(r);
            }
            return result;
        }
        catch (DoasFitException& ex)
        {
            debugInfo.errorMessage = "Error when performing DOAS fit ";
            if (ex.m_fitWindowName.size() > 0)
            {
                debugInfo.errorMessage += "(" + ex.m_fitWindowName + ")";
            }
            debugInfo.errorMessage += "Message: " + std::string(ex.what());
            return result;
        }
        catch (std::exception& ex)
        {
            debugInfo.errorMessage = ex.what();
            return result;
        }
    }

    bool RatioEvaluation::AnyFitWindowRequiresRingSpectrum() const
    {
        if (m_masterFitWindow.ringCalculation != RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING)
        {
            return true;
        }
        for (const auto& window : m_referenceFit)
        {
            if (window.ringCalculation != RING_CALCULATION_OPTION::DO_NOT_CALCULATE_RING)
            {
                return true;
            }
        }

        // no window requiring a Ring spectrum was found.
        return false;
    }

    void RatioEvaluation::ValidateSetup() const
    {
        for (const auto& window : m_referenceFit)
        {
            if (window.fitType != m_masterFitWindow.fitType)
            {
                throw std::invalid_argument("Mixed types of fit windows. All Doas fits must have the same fit type.");
            }
        }
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