#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Spectra/SpectrometerModel.h>

#include <algorithm>

#include <cstdarg>
#include <cassert>
#include <cstring>
#include <cmath>

#undef min
#undef max

CSpectrum::CSpectrum()
    : m_length(0)
{
    memset(m_data, 0, sizeof(double) * MAX_SPECTRUM_LENGTH);
}

CSpectrum::CSpectrum(const std::vector<double>& spectralData)
    : m_length((long)spectralData.size())
{
    memcpy(this->m_data, spectralData.data(), sizeof(double) * std::min(spectralData.size(), (size_t)MAX_SPECTRUM_LENGTH));
}

CSpectrum::CSpectrum(const std::vector<double>& wavelength, const std::vector<double>& spectralData)
    : m_length((long)spectralData.size()),
    m_wavelength{ begin(wavelength), end(wavelength) }
{
    memcpy(this->m_data, spectralData.data(), sizeof(double) * std::min(spectralData.size(), (size_t)MAX_SPECTRUM_LENGTH));
}

CSpectrum::CSpectrum(const CSpectrum& other)
    : m_length(other.m_length), 
      m_info(other.m_info)
{
    memcpy(this->m_data, &other.m_data, sizeof(double) * MAX_SPECTRUM_LENGTH);

    if (other.m_wavelength.size() > 0)
    {
        this->m_wavelength = std::vector<double>(begin(other.m_wavelength), end(other.m_wavelength));
    }
}

CSpectrum::CSpectrum(CSpectrum&& other)
    : m_length(other.m_length),
    m_info(std::move(other.m_info))
{
    // this still has to be copied
    memcpy(this->m_data, &other.m_data, sizeof(double) * MAX_SPECTRUM_LENGTH);

    if (other.m_wavelength.size() > 0)
    {
        this->m_wavelength = std::move(other.m_wavelength);
    }
}

CSpectrum& CSpectrum::operator=(const CSpectrum& other)
{
    this->m_info = other.m_info;
    this->m_length = other.m_length;
    memcpy(this->m_data, &other.m_data, sizeof(double) * MAX_SPECTRUM_LENGTH);

    if (other.m_wavelength.size() > 0)
    {
        this->m_wavelength = std::vector<double>(begin(other.m_wavelength), end(other.m_wavelength));
    }

    return *this;
}

CSpectrum& CSpectrum::operator=(CSpectrum&& other)
{
    this->m_info = std::move(other.m_info);
    this->m_length = other.m_length;

    // This still has to be copied..
    memcpy(this->m_data, &other.m_data, sizeof(double) * MAX_SPECTRUM_LENGTH);

    if (other.m_wavelength.size() > 0)
    {
        this->m_wavelength = std::move(other.m_wavelength);
    }

    return *this;
}

int CSpectrum::AssertRange(long &fromPixel, long &toPixel) const
{
    /* Check the input */
    assert(fromPixel >= 0 && toPixel >= 0);

    toPixel = std::min(toPixel, m_length - 1);
    fromPixel = std::min(fromPixel, toPixel - 1);

    assert(fromPixel <= toPixel);

    return 0;
}

double CSpectrum::MaxValue(long fromPixel, long toPixel) const
{
    fromPixel /= this->m_info.m_interlaceStep;
    toPixel /= this->m_info.m_interlaceStep;

    /* Check the input */
    AssertRange(fromPixel, toPixel);

    double maxv = m_data[fromPixel];
    for (long i = fromPixel + 1; i <= toPixel; ++i) {
        maxv = std::max(maxv, m_data[i]);
    }
    return maxv;
}

double CSpectrum::MinValue(long fromPixel, long toPixel) const
{
    fromPixel /= this->m_info.m_interlaceStep;
    toPixel /= this->m_info.m_interlaceStep;

    /* Check the input */
    AssertRange(fromPixel, toPixel);

    double minv = m_data[fromPixel];
    for (long i = fromPixel + 1; i <= toPixel; ++i)
    {
        minv = std::min(minv, m_data[i]);
    }
    return minv;
}

double CSpectrum::AverageValue(long fromPixel, long toPixel) const
{
    fromPixel /= this->m_info.m_interlaceStep;
    toPixel /= this->m_info.m_interlaceStep;

    /* Check the input */
    AssertRange(fromPixel, toPixel);

    double avg = m_data[fromPixel];
    for (long i = fromPixel + 1; i <= toPixel; ++i) {
        avg += m_data[i];
    }
    return (avg / (double)(toPixel - fromPixel + 1));
}

int CSpectrum::Add(const CSpectrum &spec)
{
    if (m_length != spec.m_length)
    {
        return 1;
    }

    m_length = spec.m_length;
    m_info.m_numSpec += spec.m_info.m_numSpec;

    CDateTime localCopy = spec.m_info.m_startTime;
    if (localCopy < m_info.m_startTime)
        m_info.m_startTime = localCopy;

    localCopy = spec.m_info.m_stopTime;
    if (m_info.m_stopTime < localCopy)
    {
        m_info.m_stopTime = localCopy;
    }

    return PixelwiseOperation(spec, &CSpectrum::Plus);
}
int CSpectrum::Add(const double value)
{
    return PixelwiseOperation(value, &CSpectrum::Plus);
}

int CSpectrum::Sub(const CSpectrum &spec)
{
    return PixelwiseOperation(spec, &CSpectrum::Minus);
}
int CSpectrum::Sub(const double value)
{
    return PixelwiseOperation(value, &CSpectrum::Minus);
}

int CSpectrum::Mult(const CSpectrum &spec)
{
    return PixelwiseOperation(spec, &CSpectrum::Multiply);
}
int CSpectrum::Mult(const double value)
{
    return PixelwiseOperation(value, &CSpectrum::Multiply);
}

int CSpectrum::Div(const CSpectrum &spec)
{
    return PixelwiseOperation(spec, &CSpectrum::Divide);
}
int CSpectrum::Div(const double value)
{
    return PixelwiseOperation(value, &CSpectrum::Divide);
}

// performes the supplied operation on all pixels in the two spectra
int CSpectrum::PixelwiseOperation(const CSpectrum &spec, double f(double, double))
{
    if (spec.m_length != m_length)
    {
        return 1;
    }

    for (long i = 0; i < m_length; ++i)
    {
        m_data[i] = f(m_data[i], spec.m_data[i]);
    }

    return 0;
}
//  Performs the supplied function to every pixel in the current spectrum with the supplied constant 
int CSpectrum::PixelwiseOperation(const double value, double f(double, double))
{
    for (long i = 0; i < m_length; ++i)
    {
        m_data[i] = f(m_data[i], value);
    }

    return 0;
}

double CSpectrum::GetOffset() const
{
    double avg;

    if (m_info.m_startChannel > 20)
    {
        return 0; // no idea...
    }

    // The covered pixels, where the offset is calculated
    int from = 2 / m_info.m_interlaceStep;
    int to = 20 / m_info.m_interlaceStep;

    // the offset is the average of the eighteen first pixels minus the highest one
    avg = AverageValue(from, to) * (to - from + 1);
    avg -= MaxValue(from, to);
    avg /= (to - from);

    return avg;
}

bool CSpectrum::IsDark() const
{
    // take the highest part of the (normal) spectrum
    int low = 1130 / m_info.m_interlaceStep;
    int high = 1158 / m_info.m_interlaceStep;

    // check the ranges.
    if (low > m_length)
    {
        int diff = high - low;
        high = m_length - 1;
        low = m_length - 1 - diff;
    }
    else if (high > m_length)
    {
        high = m_length;
    }

    double average = this->AverageValue(low, high);
    double offset = this->GetOffset();

    if (offset > 1e-2)
    {
        if (fabs(average - offset) < 4)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        // offset == 0 means that we could not calculate the offset
        double maxV = this->MaxValue();
        double minV = this->MinValue();
        if (maxV - minV > 20)
        {
            return false;
        }
        else
        {
            return true;
        }
    }
}

void CSpectrum::Clear()
{
    memset(m_data, 0, MAX_SPECTRUM_LENGTH * sizeof(double));
    m_length = 0;
    // uchar
    m_info.m_channel = m_info.m_flag = 0;
    // ushort
    m_info.m_startChannel = 0;
    // float
    m_info.m_compass = m_info.m_scanAngle = m_info.m_scanAngle2 = m_info.m_peakIntensity = 0;
    // long
    m_info.m_exposureTime = m_info.m_numSpec = 0;
    // CString
    m_info.m_device = "";
    m_info.m_name = "";
    // short
    m_info.m_scanIndex = m_info.m_scanSpecNum = 0;
    // GPS
    m_info.m_gps.m_altitude = 0;
    m_info.m_gps.m_latitude = m_info.m_gps.m_longitude = 0;
    // Time
    m_info.m_startTime = CDateTime();
    m_info.m_stopTime = CDateTime();
}

int	CSpectrum::Split(CSpectrum *spec[MAX_CHANNEL_NUM]) const
{
    int i;

    // If the spectrum is collected using a single channel, do nothing....
    if (Channel() < 128)
    {
        return 0;
    }

    // The number of spectra that this spectrum can be separated into
    int NSpectra = (Channel() - 127);

    // Check for illegal numbers
    if (NSpectra <= 0 || NSpectra > MAX_CHANNEL_NUM)
    {
        return 0;
    }

    // Copy the spectrum information
    for (i = 0; i < NSpectra; ++i)
    {
        spec[i]->m_info = m_info;
        spec[i]->m_info.m_interlaceStep = NSpectra;
        spec[i]->m_info.m_channel = (unsigned char)(i + 16 * (spec[i]->m_info.m_interlaceStep - 1));
        spec[i]->m_length = 0;
    }

    // Which spectrum to start with, the master or the slave channel
    //	This depends on the start-channel if a partial spectrum has been
    //	read out. By default, the odd data-points belong to the master
    //	channel and the even to the 1:st slave channel. This since the first
    //	data-point in the spectrum always is 0, and then follows a data-point
    //	from the master, one from the slave, one from the master etc....
    // TODO!!! Test this with the tripple spectormeter

    int specIndex = m_info.m_startChannel % NSpectra;

    for (i = 1; i < m_length; ++i)
    {
        if (specIndex >= 0 && specIndex < MAX_CHANNEL_NUM)
        {
            // Set the pixel data of the correct spectrum
            spec[specIndex]->m_data[spec[specIndex]->m_length] = m_data[i];

            // Change the length of the spectrum
            spec[specIndex]->m_length++;
        }
        else
        {
            // ShowMessage("CSpectrum::Split was called with an illegal start-channel");
        }

        // Take the next spectrum to update
        specIndex += 1;
        specIndex %= NSpectra;
    }

    return NSpectra;
}

int CSpectrum::GetInterlaceSteps(int channel, int &interlaceSteps)
{
    // if the spectrum is a mix of several spectra
    if (channel >= 129)
    {
        interlaceSteps = channel - 127;
        return -1;
    }

    // special case, channel = 128 is same as channel = 0
    if (channel == 128)
    {
        channel = 0;
    }

    // If the spectrum is a single spectrum
    interlaceSteps = (channel / 16) + 1; // 16->31 means interlace=2, 32->47 means interlace=3 etc.
    return (channel % 16); // the remainder tells the channel number
}

/** Interpolate the spectrum originating from the channel number 'channel' */
bool CSpectrum::InterpolateSpectrum()
{
    double	data[MAX_SPECTRUM_LENGTH];
    memset(data, 0, MAX_SPECTRUM_LENGTH * sizeof(double));
    int step = 2, start = 0;	// start is the first data-point we know in the spectrum

    // If this is not an partial spectrum, then return false
    if (m_info.m_channel < MAX_CHANNEL_NUM)
    {
        return false;
    }

    // Get the channel number of this spectrum
    switch (m_info.m_channel)
    {
    case 16:	start = 0; step = 2; break;
    case 17:	start = 1; step = 2; break;
    default:	return false;
    }

    // Get the length of this spectrum
    int newLength = m_length * step;

    // Copy the data we have
    for (int k = 0; k < m_length; ++k)
    {
        data[step*k + start] = m_data[k];
    }

    // Interpolate the data we don't have
    if (start == 1)
        data[0] = data[1];
    for (int k = start + 1; k < (newLength - 1); k += step) {
        data[k] = (data[k - 1] + data[k + 1]) * 0.5;
    }
    if (start == 0)
        data[newLength - 1] = data[newLength - 2];

    // Get the data back
    memcpy(m_data, data, newLength * sizeof(double));
    m_length = m_length * step;

    // Correct the channel number
    switch (m_info.m_channel) {
    case 16:	m_info.m_channel = 0; break;
    case 17:	m_info.m_channel = 1; break;
    case 18:	m_info.m_channel = 2; break;
    }

    return true;
}

/** Interpolate the spectrum originating from the channel number 'channel' */
bool CSpectrum::InterpolateSpectrum(CSpectrum &spec) const
{
    double	data[MAX_SPECTRUM_LENGTH];
    memset(data, 0, MAX_SPECTRUM_LENGTH * sizeof(double));
    int step = 2, start = 0;	// start is the first data-point we know in the spectrum

    // If this is not an partial spectrum, then return false
    if (m_info.m_channel < MAX_CHANNEL_NUM)
        return false;

    // Get the channel number of this spectrum
    switch (m_info.m_channel) {
    case 16:	start = 0; step = 2; break;
    case 17:	start = 1; step = 2; break;
    default:	return false;
    }

    // Get the 'new' length of this spectrum
    int newLength = m_length * step;

    // Copy the data we have
    for (int k = 0; k < m_length; ++k) {
        data[step*k + start] = m_data[k];
    }

    // Interpolate the data we don't have
    if (start == 1)
        data[0] = data[1];
    for (int k = start + 1; k < (newLength - 1); k += step) {
        data[k] = (data[k - 1] + data[k + 1]) * 0.5;
    }
    if (start == 0)
        data[newLength - 1] = data[newLength - 2];

    // Get the data back
    spec = *this;
    spec.m_length = m_length * step;
    memcpy(spec.m_data, data, newLength * sizeof(double));

    // Correct the channel number
    switch (m_info.m_channel) {
    case 16:	spec.m_info.m_channel = 0; break;
    case 17:	spec.m_info.m_channel = 1; break;
    case 18:	spec.m_info.m_channel = 2; break;
    }

    return true;
}

double GetMaximumSaturationRatioOfSpectrum(const CSpectrum& spectrum, double maximumIntensity)
{
    if (spectrum.m_info.m_numSpec > 0)
    {
        return spectrum.MaxValue(0, spectrum.m_length) / (maximumIntensity * spectrum.m_info.m_numSpec);
    }
    else
    {
        return spectrum.MaxValue(0, spectrum.m_length) / maximumIntensity;
    }
}

double GetMaximumSaturationRatioOfSpectrum(
    const CSpectrum& spectrum,
    const SpectrometerModel& model)
{
    return GetMaximumSaturationRatioOfSpectrum(spectrum, model.maximumIntensity);
}

double GetMaximumSaturationRatioOfSpectrum(const CSpectrum& spectrum)
{
    auto model = CSpectrometerDatabase::GetInstance().GetModel(spectrum.m_info.m_specModelName);
    return GetMaximumSaturationRatioOfSpectrum(spectrum, model.maximumIntensity);
}
