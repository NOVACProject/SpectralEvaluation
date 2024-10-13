#include <algorithm>
#include <fstream>
#include <sstream>
#include <limits>

#include <SpectralEvaluation/Evaluation/PlumeSpectrumSelector.h>
#include <SpectralEvaluation/File/SpectrumIO.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h>
#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>
#include <SpectralEvaluation/Spectra/SpectrometerModel.h>
#include <SpectralEvaluation/Math/SpectrumMath.h>
#include <SpectralEvaluation/Configuration/RatioEvaluationSettings.h>

namespace novac
{

std::string FormatDate(const CSpectrumInfo& spectrumInfo)
{
    char buffer[64];
    snprintf(buffer, 64, "%04d%02d%02d", spectrumInfo.m_startTime.year, spectrumInfo.m_startTime.month, spectrumInfo.m_startTime.day);
    return std::string(buffer);
}

std::string FormatTimestamp(const CDateTime& time)
{
    char buffer[64];
    snprintf(buffer, 64, "%02d%02d%02d", (int)time.hour, (int)time.minute, (int)time.second);
    return std::string(buffer);
}

std::unique_ptr<PlumeSpectra> PlumeSpectrumSelector::CreatePlumeSpectra(
    IScanSpectrumSource& originalScanFile,
    const BasicScanEvaluationResult& originalScanResult,
    const CPlumeInScanProperty& plumeProperties,
    const Configuration::RatioEvaluationSettings& settings,
    const SpectrometerModel& spectrometerModel,
    int mainSpecieIndex,
    std::string* errorMessage)
{
    m_mainSpecieIndex = mainSpecieIndex;

    auto skySpectrum = std::make_unique<CSpectrum>();
    if (originalScanFile.GetSky(*skySpectrum))
    {
        if (errorMessage != nullptr)
        {
            *errorMessage = "Failed to read sky spectrum of scan.";
        }
        return nullptr;
    }
    auto darkSpectrum = std::make_unique<CSpectrum>();
    if (originalScanFile.GetDark(*darkSpectrum))
    {
        if (errorMessage != nullptr)
        {
            *errorMessage = "Failed to read dark spectrum of scan.";
        }
        return nullptr;
    }

    if (!IsSuitableScanForRatioEvaluation(originalScanResult, plumeProperties, settings, errorMessage))
    {
        return nullptr;
    }

    auto result = std::make_unique<PlumeSpectra>();

    // Concatenate the important information on the scan, keeping only the good spectra.
    std::vector< InitialEvaluationData> evaluationData;
    ExtractEvaluationData(originalScanFile, *darkSpectrum, originalScanResult, settings, spectrometerModel, evaluationData, result->rejectedSpectrumIndices);
    if (static_cast<int>(evaluationData.size()) < settings.minNumberOfSpectraInPlume + settings.numberOfSpectraOutsideOfPlume)
    {
        if (errorMessage != nullptr)
        {
            std::stringstream msg;
            msg << "Too few spectra with good intensity in scan (" << evaluationData.size() << ") at least ";
            msg << settings.minNumberOfSpectraInPlume + settings.numberOfSpectraOutsideOfPlume << " required.";
            *errorMessage = msg.str();
        }
        return nullptr;
    }

    const bool spectraSelected = SelectSpectra(
        evaluationData,
        plumeProperties,
        settings,
        result->referenceSpectrumIndices,
        result->inPlumeSpectrumIndices,
        errorMessage);

    if (!spectraSelected)
    {
        return result;
    }

    result->referenceSpectrum = std::make_unique<CSpectrum>();
    result->inPlumeSpectrum = std::make_unique<CSpectrum>();
    result->darkSpectrum = std::make_unique<CSpectrum>();

    // Create the spectra
    AverageSpectra(originalScanFile, result->referenceSpectrumIndices, *result->referenceSpectrum, false);
    result->referenceSpectrum->m_info.m_name = "sky";
    result->referenceSpectrum->m_info.m_scanIndex = 0;

    result->darkSpectrum = std::move(darkSpectrum);
    result->darkSpectrum->m_info.m_scanIndex = 1;

    AverageSpectra(originalScanFile, result->inPlumeSpectrumIndices, *result->inPlumeSpectrum, false);
    result->inPlumeSpectrum->m_info.m_name = "plume";
    result->inPlumeSpectrum->m_info.m_scanIndex = 2;

    result->skySpectrumInfo = skySpectrum->m_info;

    return result;
}

void PlumeSpectrumSelector::CreatePlumeSpectrumFile(
    IScanSpectrumSource& originalScanFile,
    const BasicScanEvaluationResult& originalScanResult,
    const CPlumeInScanProperty& plumeProperties,
    const Configuration::RatioEvaluationSettings& settings,
    const SpectrometerModel& spectrometerModel,
    int mainSpecieIndex,
    const std::string& outputDirectory,
    std::string* errorMessage)
{

    novac::LogContext context; // TODO: Get from input

    const auto spectra = this->CreatePlumeSpectra(originalScanFile, originalScanResult, plumeProperties, settings, spectrometerModel, mainSpecieIndex, errorMessage);
    if (spectra == nullptr)
    {
        return;
    }

    std::stringstream spectrumOutputFileName;
    spectrumOutputFileName << outputDirectory << "/PlumeSpectra_" << spectra->skySpectrumInfo.m_device;
    spectrumOutputFileName << "_" << FormatDate(spectra->skySpectrumInfo);
    spectrumOutputFileName << "_" << FormatTimestamp(spectra->skySpectrumInfo.m_startTime);
    spectrumOutputFileName << "_0.pak";

    // Save the spectra as 1) sky, 2) dark and 3) inplume. This simulates the usual order of the spectra in novac data.
    CSpectrumIO spectrumWriter;
    spectrumWriter.AddSpectrumToFile(spectrumOutputFileName.str(), *spectra->referenceSpectrum, nullptr, 0, true);
    spectrumWriter.AddSpectrumToFile(spectrumOutputFileName.str(), *spectra->darkSpectrum);
    spectrumWriter.AddSpectrumToFile(spectrumOutputFileName.str(), *spectra->inPlumeSpectrum);

    std::stringstream textOutputFileName;
    textOutputFileName << outputDirectory << "/PlumeSpectra_" << spectra->skySpectrumInfo.m_device;
    textOutputFileName << "_" << FormatDate(spectra->skySpectrumInfo);
    textOutputFileName << "_" << FormatTimestamp(spectra->skySpectrumInfo.m_startTime);
    textOutputFileName << "_0.txt";

    std::ofstream textOutput(textOutputFileName.str());
    textOutput << "InPlume: " << std::endl;
    for (int idx : spectra->inPlumeSpectrumIndices)
    {
        CSpectrum spectrum;
        originalScanFile.GetSpectrum(context, (int)idx, spectrum);
        textOutput << idx << "\t" << spectrum.m_info.m_scanAngle << "\t" << FormatTimestamp(spectrum.m_info.m_startTime) << std::endl;
    }
    textOutput << std::endl;

    textOutput << "Reference: " << std::endl;
    for (int idx : spectra->referenceSpectrumIndices)
    {
        CSpectrum spectrum;
        originalScanFile.GetSpectrum(context, (int)idx, spectrum);
        textOutput << idx << "\t" << spectrum.m_info.m_scanAngle << "\t" << FormatTimestamp(spectrum.m_info.m_startTime) << std::endl;
    }
    textOutput << std::endl;

    textOutput << "Plume Properties: " << std::endl;
    textOutput << "  Completeness: " << plumeProperties.completeness << std::endl;
    textOutput << "  Center: " << plumeProperties.plumeCenter << std::endl;
    textOutput << "  Offset: " << plumeProperties.offset << std::endl;
    textOutput << "  Low Edge: " << plumeProperties.plumeEdgeLow << std::endl;
    textOutput << "  High Edge: " << plumeProperties.plumeEdgeHigh << std::endl;
    textOutput << "  HWHM Low: " << plumeProperties.plumeHalfLow << std::endl;
    textOutput << "  HWHM High: " << plumeProperties.plumeHalfHigh << std::endl;
}

void PlumeSpectrumSelector::ExtractEvaluationData(
    IScanSpectrumSource& originalScanFile,
    const CSpectrum& darkSpectrum,
    const BasicScanEvaluationResult& originalScanResult,
    const Configuration::RatioEvaluationSettings& settings,
    const SpectrometerModel& spectrometermodel,
    std::vector< InitialEvaluationData>& evaluationData,
    std::vector<std::pair<int, std::string>>& rejectedIndices) const
{
    novac::LogContext context; // TODO: Get from input

    evaluationData.clear();
    rejectedIndices.clear();

    for (size_t ii = 0; ii < originalScanResult.m_spec.size(); ++ii)
    {
        if (originalScanResult.m_spec[ii].IsBad())
        {
            rejectedIndices.push_back(std::make_pair((int)ii, "bad evaluation"));
            continue;
        }

        InitialEvaluationData data;
        data.indexInScan = static_cast<int>(ii);
        data.scanAngle = originalScanResult.m_specInfo[ii].m_scanAngle;
        data.offsetCorrectedColumn = originalScanResult.m_spec[ii].m_referenceResult[m_mainSpecieIndex].m_column;

        CSpectrum spectrum;
        if (0 == originalScanFile.GetSpectrum(context, (int)ii, spectrum))
        {
            GetIntensityOfSpectrum(spectrum, darkSpectrum, spectrometermodel, data.peakSaturation, data.peakSaturationAfterDarkCorrection);
        }

        // Already now, remove spectra which do not fulfill the intensity criteria
        if (data.peakSaturation > settings.maxSaturationRatio)
        {
            rejectedIndices.push_back(std::make_pair((int)ii, "saturated spectrum"));
        }
        else if (data.peakSaturationAfterDarkCorrection < settings.minSaturationRatio)
        {
            rejectedIndices.push_back(std::make_pair((int)ii, "dark spectrum"));
        }
        else
        {
            evaluationData.push_back(data);
        }
    }
}

bool PlumeSpectrumSelector::SelectSpectra(
    const std::vector< InitialEvaluationData>& evaluationData,
    const CPlumeInScanProperty& properties,
    const Configuration::RatioEvaluationSettings& settings,
    std::vector<int>& referenceSpectra,
    std::vector<int>& inPlumeSpectra,
    std::string* errorMessage)
{
    referenceSpectra.clear();
    inPlumeSpectra.clear();

    // Find a proprosal for the in-plume region.
    inPlumeSpectra = FindSpectraInPlume(evaluationData, properties, settings);
    std::sort(begin(inPlumeSpectra), end(inPlumeSpectra));

    if (static_cast<int>(inPlumeSpectra.size()) < settings.minNumberOfSpectraInPlume)
    {
        if (errorMessage != nullptr)
        {
            *errorMessage = "Too few good spectra in plume.";
        }
        return false;
    }

    // Get the minimum column for the selected in-plume spectra
    double minInPlumeColumn = std::numeric_limits<double>::max();
    for (int scanIdx : inPlumeSpectra)
    {
        minInPlumeColumn = std::min(minInPlumeColumn, ColumnAtScanIndex(evaluationData, scanIdx));
    }

    // Find the reference spectra as the spectra with lowest column value (ignore the sky and the dark spectra here..) 
    referenceSpectra = FindSpectraOutOfPlume(evaluationData, settings, inPlumeSpectra, minInPlumeColumn);
    std::sort(begin(referenceSpectra), end(referenceSpectra));

    if (static_cast<int>(referenceSpectra.size()) < settings.numberOfSpectraOutsideOfPlume)
    {
        if (errorMessage != nullptr)
        {
            *errorMessage = "Too few good spectra out of plume.";
        }
        return false;
    }

    return true;
}

bool PlumeSpectrumSelector::IsSuitableScanForRatioEvaluation(
    const BasicScanEvaluationResult& originalScanResult,
    const CPlumeInScanProperty& properties,
    const Configuration::RatioEvaluationSettings& settings,
    std::string* errorMessage)
{
    if (static_cast<int>(originalScanResult.m_spec.size()) < settings.minNumberOfSpectraInPlume + settings.numberOfSpectraOutsideOfPlume)
    {
        if (errorMessage != nullptr)
        {
            std::stringstream msg;
            msg << "Too few spectra in scan (" << originalScanResult.m_spec.size() << ") at least ";
            msg << settings.minNumberOfSpectraInPlume + settings.numberOfSpectraOutsideOfPlume << " required.";
            *errorMessage = msg.str();
        }
        return false;
    }
    else if (properties.completeness + 0.01 < settings.minimumPlumeCompleteness)
    {
        if (errorMessage != nullptr)
        {
            std::stringstream msg;
            msg << "Plume completeness (" << properties.completeness << ") below limit of (" << settings.minimumPlumeCompleteness << ")";
            *errorMessage = msg.str();
        }
        return false;
    }
    else if (settings.requireVisiblePlumeEdges && std::abs(properties.plumeHalfLow - NOT_A_NUMBER) < 1.0)
    {
        if (errorMessage != nullptr)
        {
            std::stringstream msg;
            msg << "Column needs to drop at least 50% on lower flank";
            *errorMessage = msg.str();
        }
        return false;
    }
    else if (settings.requireVisiblePlumeEdges && std::abs(properties.plumeHalfHigh - NOT_A_NUMBER) < 1.0)
    {
        if (errorMessage != nullptr)
        {
            std::stringstream msg;
            msg << "Column needs to drop at least 50% on upper flank";
            *errorMessage = msg.str();
        }
        return false;
    }
    else if (std::abs(properties.plumeHalfLow - NOT_A_NUMBER) < 1.0 && std::abs(properties.plumeHalfHigh - NOT_A_NUMBER) < 1.0)
    {
        if (errorMessage != nullptr)
        {
            std::stringstream msg;
            msg << "Column needs to drop at least 50% on either lower or upper flank";
            *errorMessage = msg.str();
        }
        return false;
    }

    return true;
}

double PlumeSpectrumSelector::ColumnAtScanIndex(const std::vector< InitialEvaluationData> evaluationData, int scanIndex) const
{
    for (const auto& data : evaluationData)
    {
        if (data.indexInScan == scanIndex)
        {
            return data.offsetCorrectedColumn;
        }
    }
    return 0.0;
}

std::vector<int> PlumeSpectrumSelector::FindSpectraInPlume(
    const std::vector< InitialEvaluationData>& evaluationData,
    const CPlumeInScanProperty& properties,
    const Configuration::RatioEvaluationSettings& settings)
{
    // add all spectra in the scan-angle range [plumeHalfLow, plumeHalfHigh]
    std::vector<std::pair<int, double>> spectrumColumnVsIndex;
    const double minimumScanAngle = std::max(settings.minimumScanAngle, properties.plumeHalfLow);
    const double maximumScanAngle = std::min(settings.maximumScanAngle, (properties.plumeHalfHigh > -200.0) ? properties.plumeHalfHigh : +90.0);
    for (int idx = 0; idx < (int)evaluationData.size(); ++idx)
    {
        if (evaluationData[idx].scanAngle > minimumScanAngle + 0.1 &&
            evaluationData[idx].scanAngle < maximumScanAngle - 0.1)
        {
            spectrumColumnVsIndex.push_back(std::pair<int, double>(evaluationData[idx].indexInScan, evaluationData[idx].offsetCorrectedColumn));
        }
    }
    // sort the values in order of descending column
    std::sort(begin(spectrumColumnVsIndex), end(spectrumColumnVsIndex), [&](const std::pair<int, double>& p1, const std::pair<int, double>& p2)
        {
            return p1.second > p2.second;
        });

    std::vector<int> scanIndices;
    for (const auto& item : spectrumColumnVsIndex)
    {
        scanIndices.push_back(item.first);

        if (static_cast<int>(scanIndices.size()) == settings.maxNumberOfSpectraInPlume)
        {
            return scanIndices;
        }
    }

    return scanIndices;
}

std::vector<int> PlumeSpectrumSelector::FindSpectraOutOfPlume(
    const std::vector< InitialEvaluationData>& evaluationData,
    const Configuration::RatioEvaluationSettings& settings,
    const std::vector<int>& inPlumeProposal,
    double minimumColumnInPlume)
{
    std::vector<int> indices;

    // First list all good spectra, which are not already selected to be in the plume, and then list these in order of increasing column.
    std::vector<std::pair<int, double>> spectrumColumnVsIndex;
    for (const auto& data : evaluationData)
    {
        if (std::find(begin(inPlumeProposal), end(inPlumeProposal), data.indexInScan) != end(inPlumeProposal))
        {
            continue;
        }

        // Only include spectra with an evaluation column smaller than the smallest in-plume column (in order to have contrast between in-plume and out-of-plume)
        if (data.offsetCorrectedColumn + settings.minimumInPlumeColumnDifference < minimumColumnInPlume)
        {
            spectrumColumnVsIndex.push_back(std::pair<int, double>(data.indexInScan, data.offsetCorrectedColumn));
        }
    }
    // sort the values in order of ascending column
    std::sort(begin(spectrumColumnVsIndex), end(spectrumColumnVsIndex), [&](const std::pair<int, double>& p1, const std::pair<int, double>& p2)
        {
            return p1.second < p2.second;
        });

    std::vector<int> referenceSpectraProposal;
    for (const auto& item : spectrumColumnVsIndex)
    {
        referenceSpectraProposal.push_back(item.first);

        if (static_cast<int>(referenceSpectraProposal.size()) == settings.numberOfSpectraOutsideOfPlume)
        {
            return referenceSpectraProposal;
        }
    }

    return referenceSpectraProposal;
}

void PlumeSpectrumSelector::GetIntensityOfSpectrum(
    const CSpectrum& spectrum,
    const CSpectrum& darkSpectrum,
    const SpectrometerModel& spectrometerModel,
    double& peakSaturation,
    double& peakSaturationAfterDarkCorrection) const
{
    peakSaturation = GetMaximumSaturationRatioOfSpectrum(spectrum, spectrometerModel, 1, spectrometerModel.numberOfPixels - 1);

    // dark correct and check for low intensity
    CSpectrum spectrumClone(spectrum);
    spectrumClone.Sub(darkSpectrum);

    peakSaturationAfterDarkCorrection = GetMaximumSaturationRatioOfSpectrum(spectrumClone, spectrometerModel, 1, spectrometerModel.numberOfPixels - 1);
}

}