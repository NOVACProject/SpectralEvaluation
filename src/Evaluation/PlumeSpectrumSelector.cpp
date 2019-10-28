#include <iomanip>
#include <numeric>
#include <fstream>
#include <sstream>

#include <SpectralEvaluation/Evaluation/PlumeSpectrumSelector.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h>
#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>

using namespace Evaluation;

std::string FormatDate(const CSpectrumInfo& spectrumInfo)
{
    char buffer[64];
    sprintf_s(buffer, 64, "%04d%02d%02d", spectrumInfo.m_startTime.year, spectrumInfo.m_startTime.month, spectrumInfo.m_startTime.day);
    return std::string(buffer);
}

std::string FormatTimestamp(const CDateTime& time)
{
    char buffer[64];
    sprintf_s(buffer, 64, "%02d%02d%02d", (int)time.hour, (int)time.minute, (int)time.second);
    return std::string(buffer);
}

void PlumeSpectrumSelector::CreatePlumeSpectrumFile(
    FileHandler::CScanFileHandler& originalScanFile,
    const BasicScanEvaluationResult& scanResult,
    const CPlumeInScanProperty& properties,
    int mainSpecieIndex,
    const std::string& outputDirectory)
{
    std::vector<size_t> referenceSpectrumIndices;
    std::vector<size_t> inPlumeSpectrumIndices;

    PlumeSpectrumSelectionSettings defaultSettings;

    if (!IsSuitableScanForRatioEvaluation(originalScanFile, defaultSettings, scanResult, properties))
    {
        return;
    }

    SelectSpectra(
        defaultSettings,
        originalScanFile,
        scanResult,
        properties,
        mainSpecieIndex,
        referenceSpectrumIndices,
        inPlumeSpectrumIndices);

    if (referenceSpectrumIndices.size() > 0 && inPlumeSpectrumIndices.size() > 0)
    {
        // Create and save the spectra
        CSpectrum referenceSpectrum;
        originalScanFile.AddSpectra(referenceSpectrumIndices, referenceSpectrum);
        referenceSpectrum.m_info.m_name = "sky"; // TODO: find a better name for this

        CSpectrum inPlumeSpectrum;
        originalScanFile.AddSpectra(inPlumeSpectrumIndices, inPlumeSpectrum);
        inPlumeSpectrum.m_info.m_name = "plume"; // TODO: find a better name for this

        CSpectrum darkSpectrum;
        originalScanFile.GetDark(darkSpectrum);

        CSpectrum skySpectrum;
        originalScanFile.GetSky(skySpectrum);

        std::stringstream spectrumOutputFileName;
        spectrumOutputFileName << std::setw(2) << std::setfill('0');
        spectrumOutputFileName << outputDirectory << "/PlumeSpectra_" << darkSpectrum.m_info.m_device;
        spectrumOutputFileName << "_" << FormatDate(skySpectrum.m_info);
        spectrumOutputFileName << "_" << FormatTimestamp(skySpectrum.m_info.m_startTime);
        spectrumOutputFileName << "_0.pak";

        SpectrumIO::CSpectrumIO spectrumWriter;
        spectrumWriter.AddSpectrumToFile(spectrumOutputFileName.str(), inPlumeSpectrum);
        spectrumWriter.AddSpectrumToFile(spectrumOutputFileName.str(), referenceSpectrum);
        spectrumWriter.AddSpectrumToFile(spectrumOutputFileName.str(), darkSpectrum);

        std::stringstream textOutputFileName;
        textOutputFileName << std::setw(2) << std::setfill('0');
        textOutputFileName << outputDirectory << "/PlumeSpectra_" << darkSpectrum.m_info.m_device;
        textOutputFileName << "_" << FormatDate(skySpectrum.m_info);
        textOutputFileName << "_" << FormatTimestamp(skySpectrum.m_info.m_startTime);
        textOutputFileName << "_0.txt";

        std::ofstream textOutput(textOutputFileName.str());
        textOutput << "InPlume: " << std::endl;
        for (size_t idx : inPlumeSpectrumIndices)
        {
            CSpectrum spectrum;
            originalScanFile.GetSpectrum(spectrum, (long)idx);
            textOutput << idx << "\t" << spectrum.m_info.m_scanAngle << "\t" << FormatTimestamp(spectrum.m_info.m_startTime) << std::endl;
        }

        textOutput << "Reference: " << std::endl;
        for (size_t idx : referenceSpectrumIndices)
        {
            CSpectrum spectrum;
            originalScanFile.GetSpectrum(spectrum, (long)idx);
            textOutput << idx << "\t" << spectrum.m_info.m_scanAngle << "\t" << FormatTimestamp(spectrum.m_info.m_startTime) << std::endl;
        }
        textOutput << std::endl;

        textOutput << "Plume Properties: " << std::endl;
        textOutput << "  Completeness: " << properties.completeness << std::endl;
        textOutput << "  Center: " << properties.plumeCenter << std::endl;
        textOutput << "  Offset: " << properties.offset << std::endl;
        textOutput << "  Low Edge: " << properties.plumeEdgeLow << std::endl;
        textOutput << "  High Edge: " << properties.plumeEdgeHigh << std::endl;
        textOutput << "  HWHM Low: " << properties.plumeHalfLow << std::endl;
        textOutput << "  HWHM High: " << properties.plumeHalfHigh << std::endl;

    }
}

void PlumeSpectrumSelector::SelectSpectra(
    const PlumeSpectrumSelectionSettings& settings,
    FileHandler::CScanFileHandler& scanFile,
    const BasicScanEvaluationResult& scanResult,
    const CPlumeInScanProperty& properties,
    int mainSpecieIndex,
    std::vector<size_t>& referenceSpectra,
    std::vector<size_t>& inPlumeSpectra)
{
    referenceSpectra.clear();
    inPlumeSpectra.clear();

    if (scanResult.m_spec.size() <= settings.minNumberOfSpectraOutsideOfPlume)
    {
        return; // Cannot retrieve a region, too few spectra...
    }

    // Find a proprosal for the in-plume region.
    auto inPlumeProposal = FindSpectraInPlume(scanResult, properties, mainSpecieIndex, settings);
    inPlumeProposal = FilterSpectraUsingIntensity(inPlumeProposal, scanFile, settings);

    // Step 2, find the reference region (10 adjacent spectra with lowest avg column value)
    const size_t referenceRegionWidth = 10U;
    size_t startIdxOfReferenceRegion = 0U;
    double lowestMeanColumnValue = std::numeric_limits<double>::max();
    for (size_t startIdxCandidate = 0U; startIdxCandidate < scanResult.m_spec.size() - referenceRegionWidth; ++startIdxCandidate)
    {
        const double meanColumnValue = AverageColumnValue(scanResult, mainSpecieIndex, startIdxCandidate, startIdxCandidate + referenceRegionWidth);
        if (meanColumnValue < lowestMeanColumnValue)
        {
            lowestMeanColumnValue = meanColumnValue;
            startIdxOfReferenceRegion = startIdxCandidate;
        }
    }

    referenceSpectra = std::vector<size_t>(referenceRegionWidth);
    std::iota(begin(referenceSpectra), end(referenceSpectra), (int)startIdxOfReferenceRegion);

    inPlumeSpectra = std::move(inPlumeProposal);

    return;
}

double PlumeSpectrumSelector::AverageColumnValue(
    const BasicScanEvaluationResult& scanResult,
    int specieIndex,
    size_t startIdx,
    size_t endIdx)
{
    double sum = 0.0;
    for (size_t ii = startIdx; ii < endIdx; ++ii)
    {
        sum += scanResult.m_spec[ii].m_referenceResult[specieIndex].m_column;
    }
    return sum / (double)(endIdx - startIdx);
}

bool PlumeSpectrumSelector::IsSuitableScanForRatioEvaluation(
    FileHandler::CScanFileHandler& scanFile,
    const PlumeSpectrumSelectionSettings& settings,
    const BasicScanEvaluationResult& scanResult,
    const CPlumeInScanProperty& properties)
{
    if (scanResult.m_spec.size() < settings.minNumberOfSpectraInPlume + settings.minNumberOfSpectraOutsideOfPlume)
    {
        return false; // not enough spectra
    }

    if (properties.completeness < 0.7)
    {
        return false; // need to see the plume
    }
    else if (std::abs(properties.plumeHalfLow - NOT_A_NUMBER) < 1.0 || std::abs(properties.plumeHalfHigh - NOT_A_NUMBER) < 1.0)
    {
        return false; // need to fall down to 50% of the peak-value on both sides
    }

    CSpectrum skySpectrum;
    if (scanFile.GetSky(skySpectrum))
    {
        return false;
    }

    auto model = CSpectrometerDatabase::GetInstance().GetModel(skySpectrum.m_info.m_specModelName);

    // Check the saturation ratio of the sky spectrum, making sure that it is not saturated.
    double maxSaturationRatio = GetMaximumSaturationRatioOfSpectrum(skySpectrum, model);

    if (maxSaturationRatio > settings.maxSaturationRatio)
    {
        return false;
    }

    // dark correct and check for low intensity
    CSpectrum darkSpectrum;
    if (scanFile.GetDark(darkSpectrum))
    {
        return false;
    }
    skySpectrum.Sub(darkSpectrum); // TODO: check if the number of co-adds are identical...

    maxSaturationRatio = GetMaximumSaturationRatioOfSpectrum(skySpectrum, model);
    if (maxSaturationRatio < settings.minSaturationRatio)
    {
        return false;
    }

    return true;
}

std::vector<size_t> PlumeSpectrumSelector::FindSpectraInPlume(
    const BasicScanEvaluationResult& scanResult,
    const CPlumeInScanProperty& properties,
    int mainSpecieIndex,
    const PlumeSpectrumSelectionSettings& settings)
{
    std::vector<size_t> indices;

    for (size_t idx = 0; idx < scanResult.m_spec.size(); ++idx)
    {
        if (scanResult.m_specInfo[idx].m_scanAngle > properties.plumeHalfLow + 0.1 &&
            scanResult.m_specInfo[idx].m_scanAngle < properties.plumeHalfHigh - 0.1 &&
            scanResult.m_spec[idx].m_referenceResult[mainSpecieIndex].m_column - properties.offset >= settings.minInPlumeColumn)
        {
            indices.push_back(idx);
        }
    }

    // limit the number of spectra to 10
    while (indices.size() > 10)
    {
        auto first = begin(indices);
        auto last = begin(indices) + indices.size() - 1;

        if (scanResult.m_spec[*first].m_referenceResult[mainSpecieIndex].m_column - properties.offset >
            scanResult.m_spec[*last].m_referenceResult[mainSpecieIndex].m_column - properties.offset)
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

std::vector<size_t> PlumeSpectrumSelector::FilterSpectraUsingIntensity(
    const std::vector<size_t>& proposedIndices,
    FileHandler::CScanFileHandler& scanFile,
    const PlumeSpectrumSelectionSettings& settings)
{
    std::vector<size_t> result;
    result.reserve(proposedIndices.size());

    // Get the dark-spectrum from the file. This is used to dark-correct the spectra
    //  _and_ to identify the spectrometer model for judging the maximum intensity.
    CSpectrum darkSpectrum;
    if (scanFile.GetDark(darkSpectrum))
    {
        return result; // failure.
    }
    auto model = CSpectrometerDatabase::GetInstance().GetModel(darkSpectrum.m_info.m_specModelName);

    CSpectrum spectrum;
    for (size_t idx : proposedIndices)
    {
        if (scanFile.GetSpectrum(spectrum, (long)idx))
        {
            double maxSaturationRatio = GetMaximumSaturationRatioOfSpectrum(spectrum, model);

            if (maxSaturationRatio > settings.maxSaturationRatio)
            {
                continue; // ignore this spectrum
            }

            // dark correct and check for low intensity
            spectrum.Sub(darkSpectrum); // TODO: check if the number of co-adds are identical...

            maxSaturationRatio = GetMaximumSaturationRatioOfSpectrum(spectrum, model);
            if (maxSaturationRatio < settings.minSaturationRatio)
            {
                continue; // ignore this spectrum.
            }

            result.push_back(idx);
        }
    }

    return result;
}
