#include <SpectralEvaluation/Evaluation/CrossSectionData.h>
#include <SpectralEvaluation/Evaluation/BasicMath.h>
#include <SpectralEvaluation/Fit/Vector.h>
#include <SpectralEvaluation/Fit/CubicSplineFunction.h>
#include <SpectralEvaluation/Spectra/Grid.h>
#include <fstream>
#include <numeric>

namespace Evaluation
{

CCrossSectionData::CCrossSectionData()
{
}

CCrossSectionData::CCrossSectionData(const CCrossSectionData& other)
    : m_crossSection(begin(other.m_crossSection), end(other.m_crossSection)),
      m_waveLength(begin(other.m_waveLength), end(other.m_waveLength))
{
}

CCrossSectionData &CCrossSectionData::operator=(const CCrossSectionData &other)
{
    this->m_crossSection = std::vector<double>(begin(other.m_crossSection), end(other.m_crossSection));
    this->m_waveLength   = std::vector<double>(begin(other.m_waveLength), end(other.m_waveLength));

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
    if (nullptr == crossSection) throw std::invalid_argument("Cannot set the cross section data using a null data pointer.");

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
    if(index >= m_crossSection.size())
    {
        return 0.0;
    }
    else
    {
        return m_crossSection.at(index);
    }
}

unsigned long CCrossSectionData::GetSize() const
{
    return (unsigned long)this->m_crossSection.size();
}

double CCrossSectionData::GetWavelengthAt(unsigned int index) const
{
    if(index >= m_waveLength.size())
    {
        return 0.0;
    }
    else
    {
        return m_waveLength.at(index);
    }
}

int CCrossSectionData::ReadCrossSectionFile(const std::string &fileName)
{
    this->m_waveLength.clear();
    this->m_crossSection.clear();

    std::ifstream fileRef;

    fileRef.open(fileName, std::ios_base::in);
    if(!fileRef.is_open())
    {
        std::cout << "ERROR: Cannot open reference file: %s", fileName.c_str();
        return 1;
    }

    int valuesReadNum = 0;

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


int HighPassFilter(CCrossSectionData& crossSection, bool scaleToPpmm)
{
    CBasicMath mathObject;

    const int length = (int)crossSection.m_crossSection.size();

    mathObject.Mul(crossSection.m_crossSection.data(), length, -2.5e15);
    mathObject.Delog(crossSection.m_crossSection.data(), length);
    mathObject.HighPassBinomial(crossSection.m_crossSection.data(), length, 500);
    mathObject.Log(crossSection.m_crossSection.data(), length);

    if(!scaleToPpmm)
    {
        mathObject.Div(crossSection.m_crossSection.data(), length, 2.5e15);
    }

    return 0;
}

int HighPassFilter_Ring(CCrossSectionData& crossSection)
{
    CBasicMath mathObject;

    const int length = (int)crossSection.m_crossSection.size();

    mathObject.HighPassBinomial(crossSection.m_crossSection.data(), length, 500);
    mathObject.Log(crossSection.m_crossSection.data(), length);

    return 0;
}


int Multiply(CCrossSectionData& crossSection, double scalar)
{
    CBasicMath mathObject;

    mathObject.Mul(crossSection.m_crossSection.data(), (int)crossSection.m_crossSection.size(), scalar);


    return 0;
}

int Log(CCrossSectionData& crossSection)
{
    CBasicMath mathObject;

    mathObject.Log(crossSection.m_crossSection.data(), (int)crossSection.m_crossSection.size());

    return 0;
}


void Resample(const CCrossSectionData& slf, double resolution, std::vector<double>& resampledSlf)
{
    const double xMin = slf.m_waveLength.front();
    const double xMax = slf.m_waveLength.back();

    std::vector<double> xCopy(begin(slf.m_waveLength), end(slf.m_waveLength)); // a non-const local copy
    std::vector<double> yCopy(begin(slf.m_crossSection), end(slf.m_crossSection)); // a non-const local copy

    MathFit::CVector slfX(xCopy.data(), (int)xCopy.size(), 1, false);
    MathFit::CVector slfY(yCopy.data(), (int)yCopy.size(), 1, false);

    // Create a spline from the slit-function.
    MathFit::CCubicSplineFunction spline(slfX, slfY);

    // Create a new grid for the SLF with the same resolution as the 'grid' but with the same xMin and xMax values
    UniformGrid newGridForSlf;
    newGridForSlf.minValue = slf.m_waveLength.front();
    newGridForSlf.maxValue = slf.m_waveLength.back();
    newGridForSlf.length   = 1 + (size_t)(std::round((newGridForSlf.maxValue - newGridForSlf.minValue) / resolution));

    // do the resampling...
    resampledSlf.resize(newGridForSlf.length);
    for (size_t ii = 0; ii < newGridForSlf.length; ++ii)
    {
        const double x = newGridForSlf.At(ii);

        if (x >= xMin && x <= xMax)
        {
            resampledSlf[ii] = spline.GetValue(x);
        }
        else
        {
            resampledSlf[ii] = 0.0;
        }
    }
}

void Shift(std::vector<double>& data, double pixelCount)
{
    std::vector<double> xData(data.size(), 0.0);
    std::iota(begin(xData), end(xData), 0.0);

    std::vector<double> yData(begin(data), end(data));

    MathFit::CVector slfX(xData.data(), (int)xData.size(), 1, false);
    MathFit::CVector slfY(yData.data(), (int)yData.size(), 1, false);

    CBasicMath math;
    math.ShiftAndSqueeze(slfX, slfY, 0.0, pixelCount, 1.0);

    data.resize(slfY.GetSize());
    for (int ii = 0; ii < slfY.GetSize(); ++ii)
    {
        data[ii] = slfY.GetAt(ii);
    }
}

}
