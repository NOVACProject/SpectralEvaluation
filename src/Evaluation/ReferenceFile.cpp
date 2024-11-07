#include <sstream>
#include <cmath>
#include <limits>
#include <SpectralEvaluation/Evaluation/ReferenceFile.h>
#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Calibration/ReferenceSpectrumConvolution.h>
#include <SpectralEvaluation/File/File.h>

namespace novac
{

CReferenceFile& CReferenceFile::operator=(const CReferenceFile& other)
{
    this->m_path = other.m_path;
    this->m_crossSectionFile = other.m_crossSectionFile;
    this->m_slitFunctionFile = other.m_slitFunctionFile;
    this->m_wavelengthCalibrationFile = other.m_wavelengthCalibrationFile;
    this->m_specieName = other.m_specieName;
    this->m_columnOption = other.m_columnOption;
    this->m_columnValue = other.m_columnValue;
    this->m_columnMaxValue = other.m_columnMaxValue;
    this->m_shiftOption = other.m_shiftOption;
    this->m_shiftValue = other.m_shiftValue;
    this->m_shiftMaxValue = other.m_shiftMaxValue;
    this->m_squeezeOption = other.m_squeezeOption;
    this->m_squeezeValue = other.m_squeezeValue;
    this->m_squeezeMaxValue = other.m_squeezeMaxValue;
    this->m_isFiltered = other.m_isFiltered;
    this->m_include = other.m_include;

    if (other.m_data != nullptr)
    {
        this->m_data.reset(new CCrossSectionData(*other.m_data));
    }

    return *this;
}

CReferenceFile& CReferenceFile::operator=(CReferenceFile&& other) noexcept
{
    this->m_path = std::move(other.m_path);
    this->m_crossSectionFile = std::move(other.m_crossSectionFile);
    this->m_slitFunctionFile = std::move(other.m_slitFunctionFile);
    this->m_wavelengthCalibrationFile = std::move(other.m_wavelengthCalibrationFile);
    this->m_specieName = std::move(other.m_specieName);
    this->m_columnOption = other.m_columnOption;
    this->m_columnValue = other.m_columnValue;
    this->m_columnMaxValue = other.m_columnMaxValue;
    this->m_shiftOption = other.m_shiftOption;
    this->m_shiftValue = other.m_shiftValue;
    this->m_shiftMaxValue = other.m_shiftMaxValue;
    this->m_squeezeOption = other.m_squeezeOption;
    this->m_squeezeValue = other.m_squeezeValue;
    this->m_squeezeMaxValue = other.m_squeezeMaxValue;
    this->m_isFiltered = other.m_isFiltered;
    this->m_include = other.m_include;

    if (other.m_data != nullptr)
    {
        this->m_data = std::move(other.m_data);
    }

    return *this;
}

CReferenceFile::CReferenceFile(const CReferenceFile& other)
    : m_specieName(other.m_specieName),
    m_path(other.m_path),
    m_crossSectionFile(other.m_crossSectionFile),
    m_slitFunctionFile(other.m_slitFunctionFile),
    m_wavelengthCalibrationFile(other.m_wavelengthCalibrationFile),
    m_gasFactor(other.m_gasFactor),
    m_columnOption(other.m_columnOption),
    m_columnValue(other.m_columnValue),
    m_columnMaxValue(other.m_columnMaxValue),
    m_shiftOption(other.m_shiftOption),
    m_shiftValue(other.m_shiftValue),
    m_shiftMaxValue(other.m_shiftMaxValue),
    m_squeezeOption(other.m_squeezeOption),
    m_squeezeValue(other.m_squeezeValue),
    m_squeezeMaxValue(other.m_squeezeMaxValue),
    m_isFiltered(other.m_isFiltered),
    m_include(other.m_include)
{
    if (other.m_data != nullptr)
    {
        this->m_data.reset(new CCrossSectionData(*other.m_data));
    }
}

CReferenceFile::CReferenceFile(CReferenceFile&& other) noexcept
    : m_specieName(other.m_specieName),
    m_path(other.m_path),
    m_crossSectionFile(other.m_crossSectionFile),
    m_slitFunctionFile(other.m_slitFunctionFile),
    m_wavelengthCalibrationFile(other.m_wavelengthCalibrationFile),
    m_gasFactor(other.m_gasFactor),
    m_columnOption(other.m_columnOption),
    m_columnValue(other.m_columnValue),
    m_columnMaxValue(other.m_columnMaxValue),
    m_shiftOption(other.m_shiftOption),
    m_shiftValue(other.m_shiftValue),
    m_shiftMaxValue(other.m_shiftMaxValue),
    m_squeezeOption(other.m_squeezeOption),
    m_squeezeValue(other.m_squeezeValue),
    m_squeezeMaxValue(other.m_squeezeMaxValue),
    m_isFiltered(other.m_isFiltered),
    m_include(other.m_include)
{
    if (other.m_data != nullptr)
    {
        this->m_data = std::move(other.m_data);
    }
}

CReferenceFile::CReferenceFile(const CCrossSectionData& contents)
    : m_data(new CCrossSectionData(contents))
{
}

void CReferenceFile::SetColumn(SHIFT_TYPE option, double value, double value2)
{
    this->m_columnOption = option;
    this->m_columnValue = value;

    if (SHIFT_TYPE::SHIFT_LIMIT == option)
    {
        this->m_columnMaxValue = value2;
    }
}

void CReferenceFile::SetShift(SHIFT_TYPE option, double value, double value2)
{
    this->m_shiftOption = option;
    this->m_shiftValue = value;

    if (SHIFT_TYPE::SHIFT_LIMIT == option)
    {
        this->m_shiftMaxValue = value2;
    }
}

void CReferenceFile::SetSqueeze(SHIFT_TYPE option, double value, double value2)
{
    this->m_squeezeOption = option;
    this->m_squeezeValue = value;

    if (SHIFT_TYPE::SHIFT_LIMIT == option)
    {
        this->m_squeezeMaxValue = value2;
    }
}

void CReferenceFile::ReadCrossSectionDataFromFile()
{
    if (m_path.size() == 0)
    {
        std::stringstream msg;
        msg << "Attempted to read reference for '" << this->Name() << "' from disk but the path has not been set.";
        throw InvalidReferenceException(msg.str());
    }
    else if (!IsExistingFile(m_path))
    {
        std::stringstream msg;
        msg << "Attempted to read reference file from disk but the file does not exist. Path: '" << m_path << "'";
        throw InvalidReferenceException(msg.str());
    }

    m_data.reset(new CCrossSectionData());
    if (m_data->ReadCrossSectionFile(m_path))
    {
        m_data.reset();

        std::stringstream msg;
        msg << "Attempted to read reference for '" << this->Name() << "' from disk. Path: '" << m_path << "'";
        throw InvalidReferenceException(msg.str());
    }
}

int CReferenceFile::ConvolveReference()
{
    m_data.reset(new CCrossSectionData());

    if (!::novac::ConvolveReference(m_wavelengthCalibrationFile, m_slitFunctionFile, m_crossSectionFile, *m_data))
    {
        m_data.reset();
        return 1;
    }

    return 0;
}

static std::string NameAndPathOfReference(const CReferenceFile& ref)
{
    std::stringstream message;
    message << "'" << ref.Name() << "' (" << ref.m_path << ")";
    return message.str();
}

static bool AreEqual(double a, double b, double epsilon)
{
    return (std::abs(a - b) <= epsilon * std::max(std::abs(a), std::abs(b)));
}

void CReferenceFile::VerifyReferenceValues(int fromIndex, int toIndex) const
{
    if (fromIndex < 0 || toIndex < 0 || toIndex <= fromIndex)
    {
        std::stringstream message;
        message << "Invalid call to 'VerifyReferenceValues'. From: " << fromIndex << " and To: " << toIndex << " are out of range.";
        throw InvalidReferenceException(message.str());
    }

    if (this->m_data == nullptr)
    {
        std::stringstream message;
        message << "Invalid reference " << NameAndPathOfReference(*this) << ". Data is null.";
        throw InvalidReferenceException(message.str());
    }
    if (this->m_data->m_crossSection.size() < (size_t)toIndex)
    {
        std::stringstream message;
        message << "Invalid setup of reference " << NameAndPathOfReference(*this) << ". Data does not cover fit region.";
        throw InvalidReferenceException(message.str());
    }

    fromIndex = std::max(fromIndex, 0);
    toIndex = std::min(std::max(toIndex, fromIndex + 1), int(this->m_data->GetSize()));

    double minValue = this->m_data->m_crossSection[fromIndex];
    double maxValue = this->m_data->m_crossSection[fromIndex];

    for (int idx = fromIndex; idx < toIndex; ++idx)
    {
        if (std::isnan(this->m_data->m_crossSection[idx]))
        {
            std::stringstream message;
            message << "Invalid reference " << NameAndPathOfReference(*this) << ". Data contains NaN.";
            throw InvalidReferenceException(message.str());
        }
        minValue = std::min(minValue, this->m_data->m_crossSection[idx]);
        maxValue = std::max(maxValue, this->m_data->m_crossSection[idx]);
    }

    if (AreEqual(minValue, maxValue, std::numeric_limits<double>::epsilon()))
    {
        std::stringstream message;
        message << "Invalid reference " << NameAndPathOfReference(*this) << ". Data has constant value " << maxValue << " in fit region (" << fromIndex << " to " << toIndex << ").";
        throw InvalidReferenceException(message.str());
    }

    return;
}

std::string CReferenceFile::Name() const
{
    if (this->m_specieName.length() > 0)
    {
        return this->m_specieName;
    }
    return std::string("unnamed reference");
}

}

