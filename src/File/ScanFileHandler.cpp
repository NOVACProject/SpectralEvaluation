#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/File/SpectrumIO.h>
#include <SpectralEvaluation/Spectra/SpectrometerModel.h>

#ifdef _MSC_VER
#pragma warning (push, 4)
#endif

namespace novac
{

CScanFileHandler::CScanFileHandler()
{
}

bool CScanFileHandler::CheckScanFile(const std::string& fileName)
{
    CSpectrumIO reader;
    const std::string strings[] = { std::string("sky"), std::string("zenith"), std::string("dark"), std::string("offset"), std::string("dark_cur"), std::string("darkcur") };
    int indices[] = { -1, -1, -1, -1, -1, -1 };
    bool error = false;

    m_fileName = fileName;

    // Count the number of spectra in the .pak-file
    m_specNum = reader.ScanSpectrumFile(m_fileName, strings, 6, indices);
    if (m_specNum == 0)
    {
        return false;
    }

    // Read in the spectra into the buffer, if the file is not too long
    if (m_specNum < 200)
    {
        m_spectrumBuffer.resize(m_specNum);
        FILE* f = fopen(m_fileName.c_str(), "rb");
        if (f != NULL)
        {
            CSpectrum tempSpec;
            for (unsigned int k = 0; k < m_specNum; ++k)
            {
                if (true == reader.ReadNextSpectrum(f, tempSpec))
                {
                    m_spectrumBuffer[k] = tempSpec;
                }
                else
                {
                    printf("Could not read spectrum from file: %s", fileName.c_str());
                    this->m_lastError = reader.m_lastError;
                    fclose(f);
                    return false;
                }
            }
            fclose(f);
        }
        m_spectrumBufferNum = m_specNum;
    }
    else
    {
        // The file's too large, don't store it in memory!
        m_spectrumBufferNum = 0;
        m_spectrumBuffer.clear();
    }

    // --------------- read the sky spectrum ----------------------
    int indexOfSkySpectrum = 0;
    if (indices[0] != -1)
    {
        // The spectrum named 'sky' is to be read as the sky spectrum.
        indexOfSkySpectrum = indices[0];
    }
    else if (indices[1] != -1)
    {
        // The spectrum named 'zenith' is to be read as the sky spectrum.
        indexOfSkySpectrum = indices[1];
    }
    else
    {
        // no spectrum named 'sky' or 'zenith'. Just assume that the first spectrum is the sky spectrum
        indexOfSkySpectrum = 0;
    }

    {
        std::unique_ptr<CSpectrum> spectrum = std::make_unique<CSpectrum>();
        if (true != reader.ReadSpectrum(m_fileName, indexOfSkySpectrum, *spectrum))
        {
            error = true;
        }
        else
        {
            this->m_sky = std::move(spectrum);
        }
    }

    if (error)
    {
        printf("Could not read sky-spectrum in file: %s", fileName.c_str());
        this->m_lastError = reader.m_lastError;
        return false;
    }

    // --------------- read the dark spectrum ----------------------
    int indexOfDarkSpectrum = 0;
    if (indices[2] != -1)
    {
        // If there is a dark-spectrum specified, then read it!
        indexOfDarkSpectrum = indices[2];
    }
    else if (indices[3] == -1 && indices[4] == -1)
    {
        // If there's no spectrum called 'dark' and no spectrum called 'dark_cur'
        //	and also no spectrum called 'offset', then we assume that something is wrong
        //	in the names and use the second spectrum in the scan as dark
        indexOfDarkSpectrum = 1;
    }

    {
        std::unique_ptr<CSpectrum> spectrum = std::make_unique<CSpectrum>();
        if (true != reader.ReadSpectrum(m_fileName, indexOfDarkSpectrum, *spectrum))
        {
            error = true;
        }
        else
        {
            this->m_dark = std::move(spectrum);
        }
    }

    if (error)
    {
        printf("Could not read dark-spectrum in file: %s", fileName.c_str());
        this->m_lastError = reader.m_lastError;
        return false;
    }

    // --------------- read the offset spectrum (if any) ----------------------
    if (indices[3] != -1)
    {
        std::unique_ptr<CSpectrum> spectrum = std::make_unique<CSpectrum>();
        if (true != reader.ReadSpectrum(m_fileName, indices[3], *spectrum))
        {
            printf("Could not read offset-spectrum in file: %s", fileName.c_str());
            this->m_lastError = reader.m_lastError;
            return false;
        }

        this->m_offset = std::move(spectrum);
    }

    // --------------- read the dark-current spectrum (if any) ----------------------
    if (indices[4] != -1)
    {
        std::unique_ptr<CSpectrum> spectrum = std::make_unique<CSpectrum>();
        if (true != reader.ReadSpectrum(m_fileName, indices[4], *spectrum))
        {
            printf("Could not read offset-spectrum in file: %s", fileName.c_str());
            this->m_lastError = reader.m_lastError;
            return false;
        }
        this->m_darkCurrent = std::move(spectrum);
    }
    if (indices[5] != -1)
    {
        std::unique_ptr<CSpectrum> spectrum = std::make_unique<CSpectrum>();
        if (true != reader.ReadSpectrum(m_fileName, indices[5], *spectrum))
        {
            printf("Could not read offset-spectrum in file: %s", fileName.c_str());
            this->m_lastError = reader.m_lastError;
            return false;
        }
        this->m_offset = std::move(spectrum);
    }

    // set the start and stop time of the measurement
    {
        std::unique_ptr<CSpectrum> spectrum = std::make_unique<CSpectrum>();
        if (true == reader.ReadSpectrum(m_fileName, 0, *spectrum))
        {
            this->m_startTime = spectrum->m_info.m_startTime;
            this->m_stopTime = spectrum->m_info.m_stopTime;

            // get the serial number of the spectrometer
            m_device = std::string(spectrum->m_info.m_device.c_str());

            // get the channel of the spectrometer
            m_channel = spectrum->m_info.m_channel;
        }
    }

    // set the number of spectra that we have read from the file so far
    if (m_sky != nullptr && m_sky->ScanIndex() == 0 && m_dark != nullptr && m_dark->ScanIndex() == 1)
    {
        m_specReadSoFarNum = 2;
    }
    else
    {
        m_specReadSoFarNum = 0;
    }

    // This CScanFileHandler object has been initialized
    m_initialized = true;

    return true;
}

int CScanFileHandler::GetNextSpectrum(CSpectrum& spec)
{
    CSpectrumIO reader;

    if (m_spectrumBufferNum == (unsigned int)m_specNum)
    {
        // We've read in the spectra into the buffer, just read it from there
        // instead of reading from the file itself.
        if (m_specReadSoFarNum >= m_spectrumBufferNum)
        {
            this->m_lastError = CSpectrumIO::ERROR_SPECTRUM_NOT_FOUND;
            ++m_specReadSoFarNum; // <-- go to the next spectum
            return 0;
        }
        else
        {
            spec = m_spectrumBuffer.at(m_specReadSoFarNum);
        }
    }
    else
    {
        // read the next spectrum in the file
        if (true != reader.ReadSpectrum(m_fileName, m_specReadSoFarNum, spec))
        {
            // if there was an error reading the spectrum, set the error-flag
            this->m_lastError = reader.m_lastError;
            ++m_specReadSoFarNum; // <-- go to the next spectum
            return 0;
        }
    }

    ++m_specReadSoFarNum;

    UpdateStartAndStopTimeOfScan(spec);

    // Extract the spectrometer-model from the serial-number of the spectrometer
    SpectrometerModel spectrometer = CSpectrometerDatabase::GetInstance().GuessModelFromSerial(spec.m_info.m_device);

    spec.m_info.m_average = spectrometer.averagesSpectra;
    spec.m_info.m_specModelName = spectrometer.modelName;

    return 1;
}

int CScanFileHandler::GetSpectrum(CSpectrum& spec, long specNo)
{
    if ((unsigned int)specNo < m_spectrumBufferNum)
    {
        // We've read in the spectra into the buffer, just read it from there
        // instead of reading from the file itself.
        spec = m_spectrumBuffer.at(specNo);
    }
    else
    {
        // read the desired spectrum from file
        CSpectrumIO reader;
        if (true != reader.ReadSpectrum(m_fileName, specNo, spec))
        {
            this->m_lastError = reader.m_lastError;
            return 0;
        }
    }

    UpdateStartAndStopTimeOfScan(spec);

    return 1;
}

void CScanFileHandler::UpdateStartAndStopTimeOfScan(novac::CSpectrum& spec)
{
    if (this->m_stopTime < spec.m_info.m_stopTime)
    {
        this->m_stopTime = spec.m_info.m_stopTime;
    }
    if (spec.m_info.m_startTime < this->m_startTime)
    {
        this->m_startTime = spec.m_info.m_startTime;
    }
}

int CopySpectrumIfNotNull(const std::unique_ptr<CSpectrum>& possibleNullSource, CSpectrum& destination)
{
    if (possibleNullSource == nullptr)
    {
        return 1;
    }

    destination = *possibleNullSource;
    return 0;
}

int CScanFileHandler::GetDark(CSpectrum& spec) const
{
    return CopySpectrumIfNotNull(this->m_dark, spec);
}

int CScanFileHandler::GetSky(CSpectrum& spec) const
{
    return CopySpectrumIfNotNull(this->m_sky, spec);
}

int CScanFileHandler::GetOffset(CSpectrum& spec) const
{
    return CopySpectrumIfNotNull(this->m_offset, spec);
}

int CScanFileHandler::GetDarkCurrent(CSpectrum& spec) const
{
    return CopySpectrumIfNotNull(this->m_darkCurrent, spec);
}

const CGPSData CScanFileHandler::GetGPS() const
{
    if (this->m_dark != nullptr)
    {
        return this->m_dark->GPS();
    }

    CGPSData defaultValue;
    return defaultValue;
}

double CScanFileHandler::GetCompass() const
{
    if (this->m_dark != nullptr)
    {
        return this->m_dark->Compass();
    }

    return 0.0;
}

void  CScanFileHandler::ResetCounter()
{
    m_specReadSoFarNum = 0;

    if (m_sky != nullptr && m_sky->ScanIndex() == 0)
    {
        m_specReadSoFarNum = 1;
    }

    if (m_dark != nullptr && (unsigned int)m_dark->ScanIndex() == m_specReadSoFarNum)
    {
        m_specReadSoFarNum += 1;
    }

    if (m_offset != nullptr && ((unsigned int)m_offset->ScanIndex() == m_specReadSoFarNum))
    {
        m_specReadSoFarNum += 1;
    }

    if (m_darkCurrent != nullptr && (unsigned int)m_darkCurrent->ScanIndex() == m_specReadSoFarNum)
    {
        m_specReadSoFarNum += 1;
    }
}

int CScanFileHandler::GetSpectrumNumInFile() const
{
    return 	m_specNum;
}

int	CScanFileHandler::GetInterlaceSteps() const
{
    if (!m_initialized || this->m_sky == nullptr)
    {
        return -1;
    }

    return this->m_sky->m_info.m_interlaceStep;
}

int	CScanFileHandler::GetSpectrumLength() const
{
    if (!m_initialized || this->m_sky == nullptr)
    {
        return -1;
    }

    return this->m_sky->m_length;
}

int	CScanFileHandler::GetStartChannel() const
{
    if (!m_initialized || this->m_sky == nullptr)
    {
        return -1;
    }

    return this->m_sky->m_info.m_startChannel;
}
}

#ifdef _MSC_VER
#pragma warning (pop)
#endif
