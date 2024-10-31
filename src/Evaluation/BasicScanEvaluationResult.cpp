#include <SpectralEvaluation/Evaluation/BasicScanEvaluationResult.h>
#include <SpectralEvaluation/StringUtils.h>
#include <SpectralEvaluation/Geometry.h>
#include <SpectralEvaluation/StringUtils.h>

namespace novac
{

int BasicScanEvaluationResult::AppendResult(const CEvaluationResult& evalRes, const CSpectrumInfo& specInfo)
{
    // Append the evaluationresult to the end of the 'm_spec'-vector
    m_spec.push_back(CEvaluationResult(evalRes));

    // Append the spectral information to the end of the 'm_specInfo'-vector
    m_specInfo.push_back(CSpectrumInfo(specInfo));

    // Increase the numbers of spectra in this result-set.
    ++m_specNum;
    return 0;
}

int BasicScanEvaluationResult::RemoveResult(unsigned int specIndex)
{
    if (specIndex >= m_specNum)
    {
        return 1; // not a valid index
    }

    // Remove the desired value
    m_specInfo.erase(begin(m_specInfo) + specIndex);

    // Decrease the number of values in the list
    m_specNum -= 1;

    return 0;
}

void BasicScanEvaluationResult::InitializeArrays(long size)
{
    if (size < 0 || size > 1024)
    {
        return;
    }

    m_spec.reserve(size);
    m_specInfo.reserve(size);
}

int BasicScanEvaluationResult::GetSpecieIndex(const std::string& specieName) const
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

CGPSData BasicScanEvaluationResult::GetLocation() const
{
    for (unsigned int k = 0; k < m_specNum; ++k)
    {
        const CSpectrumInfo& info = m_specInfo[k];
        if (std::abs(info.m_gps.m_longitude) > 1e-2)
        {
            return info.m_gps;
        }
    }

    return CGPSData();
}

std::vector<double> GetColumns(const BasicScanEvaluationResult& result, int specieIndex)
{
    std::vector<double> columns;
    if (specieIndex < 0 || result.m_spec.size() == 0 || result.m_spec[0].m_referenceResult.size() <= (size_t)specieIndex)
    {
        return columns;
    }

    columns.reserve(result.m_spec.size());

    for (size_t ii = 0; ii < result.m_spec.size(); ++ii)
    {
        columns.push_back(result.m_spec[ii].m_referenceResult[specieIndex].m_column);
    }

    return columns;
}

std::vector<double> GetColumnErrors(const BasicScanEvaluationResult& result, int specieIndex)
{
    std::vector<double> errors;
    if (specieIndex < 0 || result.m_spec.size() == 0 || result.m_spec[0].m_referenceResult.size() <= (size_t)specieIndex)
    {
        return errors;
    }

    errors.reserve(result.m_spec.size());

    for (size_t ii = 0; ii < result.m_spec.size(); ++ii)
    {
        errors.push_back(result.m_spec[ii].m_referenceResult[specieIndex].m_columnError);
    }

    return errors;
}

#pragma region Measurement Mode

MeasurementMode CheckMeasurementMode(const BasicScanEvaluationResult& result)
{
    if (IsStratosphereMeasurement(result))
    {
        return MeasurementMode::Stratosphere;
    }
    else if (IsWindMeasurement(result))
    {
        return MeasurementMode::Windspeed;
    }
    else if (IsDirectSunMeasurement(result))
    {
        return MeasurementMode::DirectSun;
    }
    else if (IsCompositionMeasurement(result))
    {
        return MeasurementMode::Composition;
    }
    else
    {
        return MeasurementMode::Flux;
    }
}

bool IsStratosphereMeasurement(const BasicScanEvaluationResult& result)
{
    // Check so that the measurement is long enough, but not too long
    if (result.NumberOfEvaluatedSpectra() < 3 || result.NumberOfEvaluatedSpectra() > 50)
    {
        return false;
    }

    // If the measurement started at a time when the Solar Zenith Angle 
    // was larger than 75 degrees then it is not a wind-speed measurement
    const novac::SolarPosition sun = novac::GetSunPosition(result.m_specInfo[0].m_startTime, result.GetLocation());

    // It is here assumed that the measurement is a stratospheric measurment
    // if there are more than 3 repetitions in the zenith positon
    int nRepetitions = 0; // <-- the number of repetitions in one position
    for (unsigned int k = 0; k < result.NumberOfEvaluatedSpectra(); ++k)
    {
        if (std::abs(result.m_specInfo[k].m_scanAngle) < 1e-2)
            ++nRepetitions;
        else
        {
            nRepetitions = 0;
        }

        if (nRepetitions > 50)
        {
            if (std::abs(sun.zenithAngle.value) > 75.0)
                return true;
            else
                return false;
        }
    }

    if (nRepetitions > 2)
        return true;

    return false;
}

bool IsFluxMeasurement(const BasicScanEvaluationResult& result)
{
    return (MeasurementMode::Flux == CheckMeasurementMode(result));
}

bool IsWindMeasurement(const BasicScanEvaluationResult& result)
{
    return IsWindMeasurement_Gothenburg(result) || IsWindMeasurement_Heidelberg(result);
}

bool IsWindMeasurement_Gothenburg(const BasicScanEvaluationResult& result)
{
    // Check so that the measurement is long enough
    if (result.NumberOfEvaluatedSpectra() < 52)
    {
        return false;
    }

    // If the measurement started at a time when the Solar Zenith Angle 
    // was larger than 85 degrees then it is not a wind-speed measurement
    const novac::SolarPosition sun = novac::GetSunPosition(result.m_specInfo[0].m_startTime, result.GetLocation());

    if (std::abs(sun.zenithAngle.value) >= 85.0)
    {
        return false;
    }

    // Check if this is a wind-measurement in the Gothenburg method...
    int nRepetitions = 0; // <-- the number of repetitions in one position
    float lastPos = result.m_specInfo[3].m_scanAngle;
    float lastPos2 = result.m_specInfo[3].m_scanAngle2;

    // It is here assumed that the measurement is a wind speed measurment
    // if there are more then 50 repetitions in one measurement positon
    for (unsigned int k = 4; k < result.NumberOfEvaluatedSpectra(); ++k)
    {
        float pos = result.m_specInfo[k].m_scanAngle;
        float pos2 = result.m_specInfo[k].m_scanAngle2;

        if ((std::abs(pos - lastPos) < 1e-2) && (std::abs(pos2 - lastPos2) < 1e-2))
        {
            ++nRepetitions;
        }
        else
        {
            nRepetitions = 0;
            lastPos = pos;
            lastPos2 = pos2;
        }

        if (nRepetitions > 50)
        {
            return true;
        }
    }

    return false;
}

bool IsWindMeasurement_Heidelberg(const BasicScanEvaluationResult& result)
{
    // Check so that the measurement is long enough
    if (result.NumberOfEvaluatedSpectra() < 52)
    {
        return false;
    }

    // Check if the channel-number is equal to 0
    if (result.m_specInfo[0].m_channel > 0)
    {
        return false;
    }

    // If the measurement started at a time when the Solar Zenith Angle 
    // was larger than 75 degrees then it is not a wind-speed measurement
    const novac::SolarPosition sun = novac::GetSunPosition(result.m_specInfo[0].m_startTime, result.GetLocation());

    if (std::abs(sun.zenithAngle.value) >= 75.0)
    {
        return false;
    }

    // Check if this is a wind-measurement in the Heidelberg method...
    int nRepetitions = 0; // <-- the number of repetitions in one position
    float scanAngle[2] = { result.m_specInfo[3].m_scanAngle, result.m_specInfo[4].m_scanAngle };
    float scanAngle2[2] = { result.m_specInfo[3].m_scanAngle2, result.m_specInfo[4].m_scanAngle2 };
    int  scanIndex = 0;

    // It is here assumed that the measurement is a wind speed measurement
    // if there are more then 50 repetitions in one measurement positon
    for (unsigned int k = 5; k < result.NumberOfEvaluatedSpectra(); ++k)
    {
        float pos = result.m_specInfo[k].m_scanAngle;
        float pos2 = result.m_specInfo[k].m_scanAngle2;

        if ((std::abs(pos - scanAngle[scanIndex]) < 1e-2) && (std::abs(pos2 - scanAngle2[scanIndex]) < 1e-2))
        {
            ++nRepetitions;
            scanIndex = (scanIndex + 1) % 2;
        }
        else
        {
            return false;
        }

        if (nRepetitions > 50)
        {
            return true;
        }
    }

    return false;
}

bool IsDirectSunMeasurement(const BasicScanEvaluationResult& result)
{
    // It is here assumed that the measurement is a direct-sun measurment
    // if there is at least 1 spectrum with the name 'direct_sun'
    for (unsigned int k = 5; k < result.NumberOfEvaluatedSpectra(); ++k)
    {
        if (EqualsIgnoringCase(result.m_specInfo[k].m_name, "direct_sun"))
        {
            return true;
        }
    }

    return false;
}

bool IsLunarMeasurement(const BasicScanEvaluationResult& result)
{
    int nFound = 0;

    // It is here assumed that the measurement is a lunar measurment
    // if there is at least 1 spectrum with the name 'lunar'
    for (unsigned int k = 5; k < result.NumberOfEvaluatedSpectra(); ++k)
    {
        if (EqualsIgnoringCase(result.m_specInfo[k].m_name, "lunar"))
        {
            ++nFound;
            if (nFound == 5)
            {
                return true;
            }
        }
    }

    return false;
}

bool IsCompositionMeasurement(const BasicScanEvaluationResult& result)
{
    // It is here assumed that the measurement is a composition measurment
    // if there is 
    //  * at least 1 spectrum with the name 'comp'
    //  * no more than 15 spectra in total
    if (result.NumberOfEvaluatedSpectra() > 15)
    {
        return false;
    }

    for (unsigned int k = 0; k < result.NumberOfEvaluatedSpectra(); ++k)
    {
        if (EqualsIgnoringCase(result.m_specInfo[k].m_name, "comp"))
        {
            return true;
        }
    }

    return false;
}

#pragma endregion Measurement mode

} // namespace novac
