#include <algorithm>
#include <fstream>
#include <sstream>
#include <cstdio>

#include <SpectralEvaluation/Evaluation/PlumeSpectrumSelector.h>
#include <SpectralEvaluation/File/SpectrumIO.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h>
#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>
#include <SpectralEvaluation/Spectra/SpectrometerModel.h>
#include <SpectralEvaluation/Math/SpectrumMath.h>
#include <SpectralEvaluation/Configuration/RatioEvaluationSettings.h>

using namespace novac;

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
    const BasicScanEvaluationResult& scanResult,
    const CPlumeInScanProperty& plumeProperties,
    const Configuration::RatioEvaluationSettings settings,
    int mainSpecieIndex,
    std::string* errorMessage)
{
    std::vector<int> referenceSpectrumIndices;
    std::vector<int> inPlumeSpectrumIndices;

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

    // Get some parameters regarding the scan and the spectrometer
    const auto model = CSpectrometerDatabase::GetInstance().GetModel(skySpectrum->m_info.m_specModelName);

    if (!IsSuitableScanForRatioEvaluation(*skySpectrum, *darkSpectrum, scanResult, plumeProperties, model, settings, errorMessage))
    {
        return nullptr;
    }

    SelectSpectra(
        originalScanFile,
        *darkSpectrum,
        scanResult,
        plumeProperties,
        model,
        settings,
        referenceSpectrumIndices,
        inPlumeSpectrumIndices,
        errorMessage);

    if (referenceSpectrumIndices.size() == 0 || inPlumeSpectrumIndices.size() == 0)
    {
        return nullptr;
    }

    auto result = std::make_unique<PlumeSpectra>();
    result->referenceSpectrum = std::make_unique<CSpectrum>();
    result->inPlumeSpectrum = std::make_unique<CSpectrum>();
    result->darkSpectrum = std::make_unique<CSpectrum>();

    // Create the spectra
    AverageSpectra(originalScanFile, referenceSpectrumIndices, *result->referenceSpectrum, false);
    result->referenceSpectrum->m_info.m_name = "sky";
    result->referenceSpectrum->m_info.m_scanIndex = 0;

    result->darkSpectrum = std::move(darkSpectrum);
    result->darkSpectrum->m_info.m_scanIndex = 1;

    AverageSpectra(originalScanFile, inPlumeSpectrumIndices, *result->inPlumeSpectrum, false);
    result->inPlumeSpectrum->m_info.m_name = "plume";
    result->inPlumeSpectrum->m_info.m_scanIndex = 2;

    result->referenceSpectrumIndices = std::move(referenceSpectrumIndices);
    result->inPlumeSpectrumIndices = std::move(inPlumeSpectrumIndices);

    result->skySpectrumInfo = skySpectrum->m_info;

    return result;
}

void PlumeSpectrumSelector::CreatePlumeSpectrumFile(
    IScanSpectrumSource& originalScanFile,
    const BasicScanEvaluationResult& scanResult,
    const CPlumeInScanProperty& plumeProperties,
    const Configuration::RatioEvaluationSettings settings,
    int mainSpecieIndex,
    const std::string& outputDirectory,
    std::string* errorMessage)
{
    const auto spectra = this->CreatePlumeSpectra(originalScanFile, scanResult, plumeProperties, settings, mainSpecieIndex, errorMessage);
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
        originalScanFile.GetSpectrum((int)idx, spectrum);
        textOutput << idx << "\t" << spectrum.m_info.m_scanAngle << "\t" << FormatTimestamp(spectrum.m_info.m_startTime) << std::endl;
    }
    textOutput << std::endl;

    textOutput << "Reference: " << std::endl;
    for (int idx : spectra->referenceSpectrumIndices)
    {
        CSpectrum spectrum;
        originalScanFile.GetSpectrum((int)idx, spectrum);
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

void PlumeSpectrumSelector::SelectSpectra(
    IScanSpectrumSource& scanFile,
    const CSpectrum& darkSpectrum,
    const BasicScanEvaluationResult& scanResult,
    const CPlumeInScanProperty& properties,
    const SpectrometerModel& spectrometerModel,
    const Configuration::RatioEvaluationSettings settings,
    std::vector<int>& referenceSpectra,
    std::vector<int>& inPlumeSpectra,
    std::string* errorMessage)
{
    referenceSpectra.clear();
    inPlumeSpectra.clear();

    if (static_cast<int>(scanResult.m_spec.size()) <= settings.numberOfSpectraOutsideOfPlume + settings.minNumberOfSpectraInPlume)
    {
        return; // Cannot retrieve a region, too few spectra...
    }

    // Find a proprosal for the in-plume region.
    auto inPlumeProposal = FindSpectraInPlume(scanResult, properties, settings);
    inPlumeProposal = FilterSpectraUsingIntensity(inPlumeProposal, scanFile, darkSpectrum, spectrometerModel, settings);
    if (static_cast<int>(inPlumeProposal.size()) < settings.minNumberOfSpectraInPlume)
    {
        if (errorMessage != nullptr)
        {
            *errorMessage = "Too few good spectra in plume.";
        }
        return;
    }

    // Find the reference spectra as the spectra with lowest column value (ignore the sky and the dark spectra here..) 
    auto referenceProposal = FindSpectraOutOfPlume(scanFile, darkSpectrum, scanResult, spectrometerModel, settings, inPlumeProposal);
    if (static_cast<int>(referenceProposal.size()) < settings.numberOfSpectraOutsideOfPlume)
    {
        if (errorMessage != nullptr)
        {
            *errorMessage = "Too few good spectra out of plume.";
        }
        return;
    }

    std::sort(begin(inPlumeProposal), end(inPlumeProposal));
    std::sort(begin(referenceProposal), end(referenceProposal));

    inPlumeSpectra = std::move(inPlumeProposal);
    referenceSpectra = std::move(referenceProposal);

    return;
}

bool PlumeSpectrumSelector::IsSuitableScanForRatioEvaluation(
    const CSpectrum& skySpectrum,
    const CSpectrum& darkSpectrum,
    const BasicScanEvaluationResult& scanResult,
    const CPlumeInScanProperty& properties,
    const SpectrometerModel& spectrometerModel,
    const Configuration::RatioEvaluationSettings settings,
    std::string* errorMessage)
{
    if (static_cast<int>(scanResult.m_spec.size()) < settings.minNumberOfSpectraInPlume + settings.numberOfSpectraOutsideOfPlume)
    {
        if (errorMessage != nullptr)
        {
            std::stringstream msg;
            msg << "Too few spectra in scan (" << scanResult.m_spec.size() << ") at least ";
            msg << settings.minNumberOfSpectraInPlume + settings.numberOfSpectraOutsideOfPlume << " required.";
            *errorMessage = msg.str();
        }
        return false;
    }
    else if (properties.completeness < settings.minimumPlumeCompleteness)
    {
        if (errorMessage != nullptr)
        {
            std::stringstream msg;
            msg << "Plume completeness (" << properties.completeness << ") below limit of (" << settings.minimumPlumeCompleteness << ")";
            *errorMessage = msg.str();
        }
        return false;
    }
    else if (std::abs(properties.plumeHalfLow - NOT_A_NUMBER) < 1.0)
    {
        if (errorMessage != nullptr)
        {
            std::stringstream msg;
            msg << "Column needs to drop at least 50% on lower flank";
            *errorMessage = msg.str();
        }
        return false;
    }
    else if (std::abs(properties.plumeHalfHigh - NOT_A_NUMBER) < 1.0)
    {
        if (errorMessage != nullptr)
        {
            std::stringstream msg;
            msg << "Column needs to drop at least 50% on upper flank";
            *errorMessage = msg.str();
        }
        return false;
    }

    if (!SpectrumFulfillsIntensityRequirement(skySpectrum, darkSpectrum, spectrometerModel, settings))
    {
        if (errorMessage != nullptr)
        {
            std::stringstream msg;
            msg << "Sky spectrum does not have saturation ratio in required range " << settings.minSaturationRatio << " to " << settings.maxSaturationRatio;
            *errorMessage = msg.str();
        }
        return false;
    }

    return true;
}

std::vector<int> PlumeSpectrumSelector::FindSpectraInPlume(
    const BasicScanEvaluationResult& scanResult,
    const CPlumeInScanProperty& properties,
    const Configuration::RatioEvaluationSettings settings)
{
    std::vector<int> indices;

    // add all spectra in the scan-angle range [plumeHalfLow, plumeHalfHigh]
    const double minimumScanAngle = std::max(settings.minimumScanAngle, properties.plumeHalfLow);
    const double maximumScanAngle = std::min(settings.maximumScanAngle, properties.plumeHalfHigh);
    for (int idx = 0; idx < (int)scanResult.m_spec.size(); ++idx)
    {
        if (scanResult.m_specInfo[idx].m_scanAngle > minimumScanAngle + 0.1 &&
            scanResult.m_specInfo[idx].m_scanAngle < maximumScanAngle - 0.1)
        {
            indices.push_back(idx);
        }
    }

    // limit the number of spectra
    while (indices.size() > static_cast<size_t>(settings.maxNumberOfSpectraInPlume))
    {
        auto first = begin(indices);
        auto last = begin(indices) + indices.size() - 1;

        if (scanResult.m_spec[*first].m_referenceResult[m_mainSpecieIndex].m_column >
            scanResult.m_spec[*last].m_referenceResult[m_mainSpecieIndex].m_column)
        {
            indices.erase(last);
        }
        else
        {
            indices.erase(first);
        }
    }

    return indices;
}

std::vector<int> PlumeSpectrumSelector::FindSpectraOutOfPlume(
    IScanSpectrumSource& scanFile,
    const CSpectrum& darkSpectrum,
    const BasicScanEvaluationResult& scanResult,
    const SpectrometerModel& spectrometerModel,
    const Configuration::RatioEvaluationSettings settings,
    const std::vector<int>& inPlumeProposal)
{
    std::vector<int> indices;

    // First list all good spectra, which are not already selected to be in the plume, and then list these in order of increasing column.
    std::vector<std::pair<int, double>> spectrumColumnVsIndex;
    for (int idx = 2U; idx < scanResult.m_spec.size(); ++idx)
    {
        if (scanResult.m_spec[idx].IsBad())
        {
            continue;
        }
        if (std::find(begin(inPlumeProposal), end(inPlumeProposal), idx) != end(inPlumeProposal))
        {
            continue;
        }
        spectrumColumnVsIndex.push_back(std::pair<int, double>(idx, scanResult.m_spec[idx].m_referenceResult[m_mainSpecieIndex].m_column));
    }
    std::sort(begin(spectrumColumnVsIndex), end(spectrumColumnVsIndex), [&](const std::pair<int, double>& p1, const std::pair<int, double>& p2)
    {
        return p1.second < p2.second;
    });

    std::vector<int> referenceSpectraProposal;
    CSpectrum spectrum;
    for (const auto& item : spectrumColumnVsIndex)
    {
        int idx = item.first;

        if (0 == scanFile.GetSpectrum((int)idx, spectrum) && SpectrumFulfillsIntensityRequirement(spectrum, darkSpectrum, spectrometerModel, settings))
        {
            referenceSpectraProposal.push_back(idx);

            if (referenceSpectraProposal.size() == static_cast<int>(settings.numberOfSpectraOutsideOfPlume))
            {
                return referenceSpectraProposal;
            }
        }
    }

    return referenceSpectraProposal;
}

bool PlumeSpectrumSelector::SpectrumFulfillsIntensityRequirement(
    const CSpectrum& spectrum,
    const CSpectrum& darkSpectrum,
    const SpectrometerModel& spectrometerModel,
    const Configuration::RatioEvaluationSettings settings)
{
    double maxSaturationRatio = GetMaximumSaturationRatioOfSpectrum(spectrum, spectrometerModel);
    if (maxSaturationRatio > settings.maxSaturationRatio)
    {
        return false;
    }

    // dark correct and check for low intensity
    CSpectrum spectrumClone(spectrum);
    spectrumClone.Sub(darkSpectrum);

    maxSaturationRatio = GetMaximumSaturationRatioOfSpectrum(spectrum, spectrometerModel);
    if (maxSaturationRatio < settings.minSaturationRatio)
    {
        return false;
    }

    return true;
}

std::vector<int> PlumeSpectrumSelector::FilterSpectraUsingIntensity(
    const std::vector<int>& proposedIndices,
    IScanSpectrumSource& scanFile,
    const CSpectrum& darkSpectrum,
    const SpectrometerModel& spectrometerModel,
    const Configuration::RatioEvaluationSettings settings)
{
    std::vector<int> result;
    result.reserve(proposedIndices.size());

    CSpectrum spectrum;
    for (int idx : proposedIndices)
    {
        if (0 == scanFile.GetSpectrum((int)idx, spectrum) && SpectrumFulfillsIntensityRequirement(spectrum, darkSpectrum, spectrometerModel, settings))
        {
            result.push_back(idx);
        }
    }

    return result;
}
