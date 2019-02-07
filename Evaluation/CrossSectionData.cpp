#include "CrossSectionData.h"
#include <fstream>

namespace Evaluation
{

CCrossSectionData::CCrossSectionData()
{
}

CCrossSectionData::~CCrossSectionData()
{
}

CCrossSectionData &CCrossSectionData::operator=(const CCrossSectionData &xs2)
{
    this->m_length = xs2.m_length;

    // copy the data of the arrays
    this->m_crossSection = std::vector<double>(begin(xs2.m_crossSection), end(xs2.m_crossSection));
    this->m_waveLength   = std::vector<double>(begin(xs2.m_waveLength), end(xs2.m_waveLength));

    return *this;
}


void SetAtGrow(std::vector<double>& v, int index, double value)
{
    if ((size_t)index >= v.size())
    {
        v.resize(index + 1);
    }
    v[index] = value;
}

/** Sets the reference information at the given pixel */
void CCrossSectionData::SetAt(int index, double wavel, double value)
{
    SetAtGrow(m_waveLength, index, wavel);
    SetAtGrow(m_crossSection, index, value);
}

void CCrossSectionData::Set(double *wavelength, double *crossSection, unsigned long pointNum)
{
    if(nullptr == wavelength) throw std::invalid_argument("Cannot set the cross section data using a null wavelength data pointer.");
    if (nullptr == wavelength) throw std::invalid_argument("Cannot set the cross section data using a null data pointer.");

    this->m_length = pointNum;
    m_waveLength.resize(pointNum);
    m_crossSection.resize(pointNum);

    for(unsigned int k = 0; k < pointNum; ++k)
    {
        this->m_waveLength[k]   = wavelength[k];
        this->m_crossSection[k] = crossSection[k];
    }
}

void CCrossSectionData::Set(double *crossSection, unsigned long pointNum)
{
    this->m_length = pointNum;

    m_waveLength.resize(pointNum);
    m_crossSection.resize(pointNum);

    for(unsigned int k = 0; k < pointNum; ++k)
    {
        double lambda = (double)k;

        this->m_waveLength[k]   = lambda;
        this->m_crossSection[k] = crossSection[k];
    }
}

void CCrossSectionData::Set(MathFit::CVector &crossSection, unsigned long pointNum)
{
    this->m_length = pointNum;

    m_waveLength.resize(pointNum);
    m_crossSection.resize(pointNum);

    for(unsigned int k = 0; k < pointNum; ++k)
    {
        const double value       = crossSection.GetAt(k);
        this->m_crossSection[k]  = value;
    }
}

double CCrossSectionData::GetAt(unsigned int index) const
{
    if(index > m_length) {
        return 0.0;
    } else {
        return m_crossSection.at(index);
    }
}

unsigned long CCrossSectionData::GetSize() const
{
    return this->m_length;
}

double CCrossSectionData::GetWavelengthAt(unsigned int index) const
{
    if(index > this->m_length) {
        return 0.0;
    } else{
        return m_waveLength.at(index);
    }
}

int CCrossSectionData::ReadCrossSectionFile(const std::string &fileName)
{
    this->m_waveLength.clear();
    this->m_crossSection.clear();
    this->m_length = 0U;

    std::ifstream fileRef;

    fileRef.open(fileName, std::ios_base::in);
    if(!fileRef.is_open())
    {
        std::cout << "ERROR: Cannot open reference file: %s", fileName.c_str();
        return 1;
    }

    int valuesReadNum = 0;
    int nColumns = 1;

    // read reference spectrum into the 'fValue's array
    const int maxSize = 65536;
    std::vector<char> tmpBuffer(maxSize);
    while(!fileRef.eof())
    {
        fileRef.getline(tmpBuffer.data(), maxSize);

        // this construction enables us to read files with both one or two columns
        double fValue1 = 0.0;
        double fValue2 = 0.0;
        int nColumns   = sscanf(tmpBuffer.data(), "%lf\t%lf", &fValue1, &fValue2);

        // check so that we actually could read the data
        if(nColumns < 1 || nColumns > 2)
        {
            break;
        }

        ++valuesReadNum;

        if (nColumns == 2)
        {
            m_waveLength.push_back(fValue1);
            m_crossSection.push_back(fValue2);
        }
        else
        {
            m_crossSection.push_back(fValue1);
        }
    }

    if (valuesReadNum == 0)
    {
        return 1; // failed to read any lines
    }
    fileRef.close();

    return 0;
}

}
