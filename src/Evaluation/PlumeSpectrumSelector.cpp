#include <algorithm>
#include <fstream>
#include <sstream>
#include <cstdio>

#include <SpectralEvaluation/Evaluation/PlumeSpectrumSelector.h>
#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h>
#include <SpectralEvaluation/Flux/PlumeInScanProperty.h>
#include <SpectralEvaluation/Spectra/SpectrometerModel.h>

using namespace Evaluation;

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

void PlumeSpectrumSelector::CreatePlumeSpectrumFile(
    FileHandler::CScanFileHandler& originalScanFile,
    const BasicScanEvaluationResult& scanResult,
    const CPlumeInScanProperty& properties,
    int mainSpecieIndex,
    const std::string& outputDirectory)
{
    std::vector<size_t> referenceSpectrumIndices;
    std::vector<size_t> inPlumeSpectrumIndices;

    m_settings = PlumeSpectrumSelectionSettings();
    m_mainSpecieIndex = mainSpecieIndex;

    CSpectrum skySpectrum;
    if (originalScanFile.GetSky(skySpectrum))
    {
        return; // cannot read spectra from the file, ignore it.
    }
    CSpectrum darkSpectrum;
    if (originalScanFile.GetDark(darkSpectrum))
    {
        return; // cannot read spectra from the file, ignore it.
    }

    // Get some parameters regarding the scan and the spectrometer
    auto model = CSpectrometerDatabase::GetInstance().GetModel(skySpectrum.m_info.m_specModelName);
    this->m_maximumSpectrometerIntensity = model.maximumIntensity;

    if (!IsSuitableScanForRatioEvaluation(skySpectrum, darkSpectrum, scanResult, properties))
    {
        return;
    }

    SelectSpectra(
        originalScanFile,
        darkSpectrum,
        scanResult,
        properties,
        referenceSpectrumIndices,
        inPlumeSpectrumIndices);

    if (referenceSpectrumIndices.size() > 0 && inPlumeSpectrumIndices.size() > 0)
    {
        // Create and save the spectra
        CSpectrum referenceSpectrum;
        originalScanFile.AddSpectra(referenceSpectrumIndices, referenceSpectrum);
        referenceSpectrum.m_info.m_name = "sky";
        referenceSpectrum.m_info.m_scanIndex = 0;

        darkSpectrum.m_info.m_scanIndex = 1;

        CSpectrum inPlumeSpectrum;
        originalScanFile.AddSpectra(inPlumeSpectrumIndices, inPlumeSpectrum);
        inPlumeSpectrum.m_info.m_name = "plume";
        referenceSpectrum.m_info.m_scanIndex = 2;

        std::stringstream spectrumOutputFileName;
        spectrumOutputFileName << outputDirectory << "/PlumeSpectra_" << darkSpectrum.m_info.m_device;
        spectrumOutputFileName << "_" << FormatDate(skySpectrum.m_info);
        spectrumOutputFileName << "_" << FormatTimestamp(skySpectrum.m_info.m_startTime);
        spectrumOutputFileName << "_0.pak";

        SpectrumIO::CSpectrumIO spectrumWriter;
        spectrumWriter.AddSpectrumToFile(spectrumOutputFileName.str(), referenceSpectrum, nullptr, 0, true);
        spectrumWriter.AddSpectrumToFile(spectrumOutputFileName.str(), darkSpectrum);
        spectrumWriter.AddSpectrumToFile(spectrumOutputFileName.str(), inPlumeSpectrum);

        std::stringstream textOutputFileName;
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
        textOutput << std::endl;

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
    FileHandler::CScanFileHandler& scanFile,
    const CSpectrum& darkSpectrum,
    const BasicScanEvaluationResult& scanResult,
    const CPlumeInScanProperty& properties,
    std::vector<size_t>& referenceSpectra,
    std::vector<size_t>& inPlumeSpectra)
{
    referenceSpectra.clear();
    inPlumeSpectra.clear();

    if (scanResult.m_spec.size() <= m_settings.numberOfSpectraOutsideOfPlume + m_settings.minNumberOfSpectraInPlume)
    {
        return; // Cannot retrieve a region, too few spectra...
    }

    // Find a proprosal for the in-plume region.
    auto inPlumeProposal = FindSpectraInPlume(scanResult, properties);
    inPlumeProposal = FilterSpectraUsingIntensity(inPlumeProposal, scanFile, darkSpectrum);
    if (inPlumeProposal.size() < m_settings.minNumberOfSpectraInPlume)
    {
        return;
    }

    // Find the reference spectra as the spectra with lowest column value (ignore the sky and the dark spectra here..) 
    auto referenceProposal = FindSpectraOutOfPlume(scanFile, darkSpectrum, scanResult, inPlumeProposal);
    if (referenceProposal.size() < m_settings.numberOfSpectraOutsideOfPlume)
    {
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
    const CPlumeInScanProperty& properties)
{
    if (scanResult.m_spec.size() < m_settings.minNumberOfSpectraInPlume + m_settings.numberOfSpectraOutsideOfPlume)
    {
        return false; // not enough spectra
    }
    else if (properties.completeness < 0.7)
    {
        return false; // need to see the plume
    }
    else if (std::abs(properties.plumeHalfLow - NOT_A_NUMBER) < 1.0 || std::abs(properties.plumeHalfHigh - NOT_A_NUMBER) < 1.0)
    {
        return false; // need to fall down to 50% of the peak-value on both sides
    }

    if (!SpectrumFulfillsIntensityRequirement(skySpectrum, darkSpectrum))
    {
        return false;
    }

    return true;
}

std::vector<size_t> PlumeSpectrumSelector::FindSpectraInPlume(
    const BasicScanEvaluationResult& scanResult,
    const CPlumeInScanProperty& properties)
{
    std::vector<size_t> indices;

    // add all spectra in the scan-angle range [plumeHalfLow, plumeHalfHigh]
    const double minimumScanAngle = std::max(m_settings.minimumScanAngle, properties.plumeHalfLow);
    const double maximumScanAngle = std::min(m_settings.maximumScanAngle, properties.plumeHalfHigh);
    for (size_t idx = 0; idx < scanResult.m_spec.size(); ++idx)
    {
        if (scanResult.m_specInfo[idx].m_scanAngle > minimumScanAngle + 0.1 &&
            scanResult.m_specInfo[idx].m_scanAngle < maximumScanAngle - 0.1)
        {
            indices.push_back(idx);
        }
    }

    // limit the number of spectra
    while (indices.size() > m_settings.maxNumberOfSpectraInPlume)
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

std::vector<size_t> PlumeSpectrumSelector::FindSpectraOutOfPlume(
    FileHandler::CScanFileHandler& scanFile,
    const CSpectrum& darkSpectrum,
    const BasicScanEvaluationResult& scanResult,
    const std::vector<size_t>& inPlumeProposal)
{
    std::vector<size_t> indices;

    std::vector<std::pair<size_t, double>> spectrumColumnVsIndex;
    for (size_t idx = 2U; idx < scanResult.m_spec.size(); ++idx)
    {
        if (std::find(begin(inPlumeProposal), end(inPlumeProposal), idx) == end(inPlumeProposal))
        {
            spectrumColumnVsIndex.push_back(std::pair<size_t, double>(idx, scanResult.m_spec[idx].m_referenceResult[m_mainSpecieIndex].m_column));
        }
    }
    std::sort(begin(spectrumColumnVsIndex), end(spectrumColumnVsIndex), [&](const std::pair<size_t, double>& p1, const std::pair<size_t, double>& p2)
    {
        return p1.second < p2.second;
    });

    std::vector<size_t> referenceSpectraProposal;
    CSpectrum spectrum;
    for (auto item : spectrumColumnVsIndex)
    {
        size_t idx = item.first;

        if (scanFile.GetSpectrum(spectrum, (long)idx) && SpectrumFulfillsIntensityRequirement(spectrum, darkSpectrum))
        {
            referenceSpectraProposal.push_back(idx);

            if (referenceSpectraProposal.size() == m_settings.numberOfSpectraOutsideOfPlume)
            {
                return referenceSpectraProposal;
            }
        }
    }

    return referenceSpectraProposal;
}

bool PlumeSpectrumSelector::SpectrumFulfillsIntensityRequirement(
    const CSpectrum& spectrum,
    const CSpectrum& darkSpectrum)
{
    double maxSaturationRatio = GetMaximumSaturationRatioOfSpectrum(spectrum, this->m_maximumSpectrometerIntensity);

    if (maxSaturationRatio > m_settings.maxSaturationRatio)
    {
        return false;
    }

    // dark correct and check for low intensity
    CSpectrum spectrumClone(spectrum);
    spectrumClone.Sub(darkSpectrum);

    maxSaturationRatio = GetMaximumSaturationRatioOfSpectrum(spectrum, this->m_maximumSpectrometerIntensity);
    if (maxSaturationRatio < m_settings.minSaturationRatio)
    {
        return false;
    }

    return true;
}

std::vector<size_t> PlumeSpectrumSelector::FilterSpectraUsingIntensity(
    const std::vector<size_t>& proposedIndices,
    FileHandler::CScanFileHandler& scanFile,
    const CSpectrum& darkSpectrum)
{
    std::vector<size_t> result;
    result.reserve(proposedIndices.size());

    CSpectrum spectrum;
    for (size_t idx : proposedIndices)
    {
        if (scanFile.GetSpectrum(spectrum, (long)idx) && SpectrumFulfillsIntensityRequirement(spectrum, darkSpectrum))
        {
            result.push_back(idx);
        }
    }

    return result;
}
