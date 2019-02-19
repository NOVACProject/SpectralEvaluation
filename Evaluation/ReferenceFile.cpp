#include "ReferenceFile.h"
#include "CrossSectionData.h"
#include "../Spectra/ReferenceSpectrumConvolution.h"

namespace Evaluation
{

CReferenceFile::CReferenceFile()
{
}
  
CReferenceFile::~CReferenceFile()
{
}

/** assignment operator */
CReferenceFile &CReferenceFile::operator=(const CReferenceFile &other)
{
    this->m_path                        = std::string(other.m_path);
    this->m_crossSectionFile            = std::string(other.m_crossSectionFile);
    this->m_slitFunctionFile            = std::string(other.m_slitFunctionFile);
    this->m_wavelengthCalibrationFile   = std::string(other.m_wavelengthCalibrationFile);

    this->m_specieName          = std::string(other.m_specieName);
    this->m_columnOption        = other.m_columnOption;
    this->m_columnValue         = other.m_columnValue;
    this->m_columnMaxValue      = other.m_columnMaxValue;
    this->m_shiftOption         = other.m_shiftOption;
    this->m_shiftValue          = other.m_shiftValue;
    this->m_shiftMaxValue       = other.m_shiftMaxValue;
    this->m_squeezeOption       = other.m_squeezeOption;
    this->m_squeezeValue        = other.m_squeezeValue;
    this->m_squeezeMaxValue     = other.m_squeezeMaxValue;
    this->m_isFiltered          = other.m_isFiltered;
    
    if (other.m_data != nullptr)
    {
        this->m_data = new CCrossSectionData(*other.m_data);
    }

    return *this;
}

CReferenceFile::CReferenceFile(const CReferenceFile& other)
{
    this->m_path                        = std::string(other.m_path);
    this->m_crossSectionFile            = std::string(other.m_crossSectionFile);
    this->m_slitFunctionFile            = std::string(other.m_slitFunctionFile);
    this->m_wavelengthCalibrationFile   = std::string(other.m_wavelengthCalibrationFile);

    this->m_specieName      = std::string(other.m_specieName);
    this->m_columnOption    = other.m_columnOption;
    this->m_columnValue     = other.m_columnValue;
    this->m_columnMaxValue  = other.m_columnMaxValue;
    this->m_shiftOption     = other.m_shiftOption;
    this->m_shiftValue      = other.m_shiftValue;
    this->m_shiftMaxValue   = other.m_shiftMaxValue;
    this->m_squeezeOption   = other.m_squeezeOption;
    this->m_squeezeValue    = other.m_squeezeValue;
    this->m_squeezeMaxValue = other.m_squeezeMaxValue;
    this->m_isFiltered      = other.m_isFiltered;

    if (other.m_data != nullptr)
    {
        this->m_data = new CCrossSectionData(*other.m_data);
    }
}


/** Setting the column */
void CReferenceFile::SetColumn(SHIFT_TYPE option, double value, double value2)
{
    this->m_columnOption = option;
    this->m_columnValue  = value;

    if(SHIFT_LIMIT == option)
    {
        this->m_columnMaxValue = value2;
    }
}

/** Setting the shift */
void CReferenceFile::SetShift(SHIFT_TYPE option, double value, double value2)
{
    this->m_shiftOption = option;
    this->m_shiftValue  = value;

    if(SHIFT_LIMIT == option)
    {
        this->m_shiftMaxValue = value2;
    }
}

/** Setting the squeeze */
void CReferenceFile::SetSqueeze(SHIFT_TYPE option, double value, double value2)
{
    this->m_squeezeOption = option;
    this->m_squeezeValue  = value;

    if(SHIFT_LIMIT == option)
    {
        this->m_squeezeMaxValue = value2;
    }
}

int CReferenceFile::ReadCrossSectionDataFromFile()
{
    if (m_path.size() == 0)
    {
        return 1;
    }

    if (m_data != nullptr)
    {
        delete m_data;
    }

    m_data = new CCrossSectionData();
    if (m_data->ReadCrossSectionFile(m_path))
    {
        delete m_data;
        m_data = nullptr;
        return 1;
    }

    return 0;
}

int CReferenceFile::ConvolveReference()
{
    if (m_data != nullptr)
    {
        delete m_data;
    }
    m_data = new CCrossSectionData();

    if (! Evaluation::ConvolveReference(m_wavelengthCalibrationFile, m_slitFunctionFile, m_crossSectionFile, *m_data))
    {
        delete m_data;
        m_data = nullptr;
        return 1;
    }

    return 0;
}

}
