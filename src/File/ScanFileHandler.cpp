#include <SpectralEvaluation/File/ScanFileHandler.h>
#include <SpectralEvaluation/Spectra/SpectrometerModel.h>

#ifdef _MSC_VER
#pragma warning (push, 4)
#endif

using namespace SpectrumIO;

namespace FileHandler
{

CScanFileHandler::CScanFileHandler()
{
    m_specReadSoFarNum = 0;
    m_initialized = false;
    m_channel = 0;
    m_specNum = 0;

    m_spectrumBufferNum = 0;

    m_fHasDark = true;
    m_fHasSky = true;
    m_fHasOffset = false;
    m_fHasDarkCurrent = false;
}

CScanFileHandler::~CScanFileHandler()
{
    m_spectrumBuffer.clear();
}

bool CScanFileHandler::CheckScanFile(const std::string& fileName)
{
    CSpectrumIO reader;
    const std::string strings[] = { std::string("sky"), std::string("zenith"), std::string("dark"), std::string("offset"), std::string("dark_cur"), std::string("darkcur") };
    int indices[] = { -1, -1, -1, -1, -1, -1 };
    bool error = false;
    CSpectrum tempSpec;

    m_fileName = fileName;

    // Count the number of spectra in the .pak-file
    m_specNum = reader.ScanSpectrumFile(m_fileName, strings, 6, indices);

    // Error checking on the contents of the file
    if (m_specNum <= 0)
    {
        return false;
    }

    // Read in the spectra into the buffer, if the file is not too long
    if (m_specNum < 200)
    {
        m_spectrumBuffer.resize(m_specNum);
        FILE *f = fopen(m_fileName.c_str(), "rb");
        if (f != NULL)
        {
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
    if (indices[0] != -1) {
        if (true != reader.ReadSpectrum(m_fileName, indices[0], m_sky)) {
            error = true;
            m_fHasSky = false;
        }
    }
    else if (indices[1] != -1) {
        if (true != reader.ReadSpectrum(m_fileName, indices[1], m_sky)) {
            error = true;
            m_fHasSky = false;
        }
    }
    else if (true != reader.ReadSpectrum(m_fileName, 0, m_sky)) {
        error = true;
        m_fHasSky = false;
    }
    if (error) {
        printf("Could not read sky-spectrum in file: %s", fileName.c_str());
        this->m_lastError = reader.m_lastError;
        return false;
    }

    // --------------- read the dark spectrum ----------------------
    if (indices[2] != -1) {
        // If there is a dark-spectrum specified, then read it!
        if (true != reader.ReadSpectrum(m_fileName, indices[2], m_dark)) {
            error = true;
        }
    }
    else if (indices[3] == -1 && indices[4] == -1) {
        // If there's no spectrum called 'dark' and no spectrum called 'dark_cur'
        //	and also no spectrum called 'offset', then we assume that something is wrong
        //	in the names and use the second spectrum in the scan as dark
        if (true != reader.ReadSpectrum(m_fileName, 1, m_dark)) {
            error = true;
        }
    }
    if (error) {
        m_fHasDark = false;
        printf("Could not read dark-spectrum in file: %s", fileName.c_str());
        this->m_lastError = reader.m_lastError;
        return false;
    }

    // --------------- read the offset spectrum (if any) ----------------------
    if (indices[3] != -1) {
        if (true != reader.ReadSpectrum(m_fileName, indices[3], m_offset)) {
            printf("Could not read offset-spectrum in file: %s", fileName.c_str());
            this->m_lastError = reader.m_lastError;
            return false;
        }
        m_fHasOffset = true;
    }

    // --------------- read the dark-current spectrum (if any) ----------------------
    if (indices[4] != -1) {
        if (true != reader.ReadSpectrum(m_fileName, indices[4], m_darkCurrent)) {
            printf("Could not read offset-spectrum in file: %s", fileName.c_str());
            this->m_lastError = reader.m_lastError;
            return false;
        }
        m_fHasDarkCurrent = true;
    }
    if (indices[5] != -1) {
        if (true != reader.ReadSpectrum(m_fileName, indices[5], m_darkCurrent)) {
            printf("Could not read offset-spectrum in file: %s", fileName.c_str());
            this->m_lastError = reader.m_lastError;
            return false;
        }
        m_fHasDarkCurrent = true;
    }
    // set the start and stop time of the measurement
    if (true == reader.ReadSpectrum(m_fileName, 0, tempSpec)) {
        this->m_startTime = tempSpec.m_info.m_startTime;
        this->m_stopTime = tempSpec.m_info.m_stopTime;

        // get the serial number of the spectrometer
        m_device = std::string(tempSpec.m_info.m_device.c_str());

        // get the channel of the spectrometer
        m_channel = tempSpec.m_info.m_channel;
    }

    // set the number of spectra that we have read from the file so far
    if (m_sky.ScanIndex() == 0 && m_dark.ScanIndex() == 1)
        m_specReadSoFarNum = 2;
    else
        m_specReadSoFarNum = 0;

    // This CScanFileHandler object has been initialized
    m_initialized = true;

    return true;
}

/** Returns the next spectrum in the scan */
int CScanFileHandler::GetNextSpectrum(CSpectrum &spec) {
    CSpectrumIO reader;

    if (m_spectrumBufferNum == (unsigned int)m_specNum) {
        // We've read in the spectra into the buffer, just read it from there
        // instead of reading from the file itself.
        if (m_specReadSoFarNum >= m_spectrumBufferNum) {
            this->m_lastError = SpectrumIO::CSpectrumIO::ERROR_SPECTRUM_NOT_FOUND;
            ++m_specReadSoFarNum; // <-- go to the next spectum
            return 0;
        }
        else {
            spec = m_spectrumBuffer.at(m_specReadSoFarNum);
        }
    }
    else {
        // read the next spectrum in the file
        if (true != reader.ReadSpectrum(m_fileName, m_specReadSoFarNum, spec)) {
            // if there was an error reading the spectrum, set the error-flag
            this->m_lastError = reader.m_lastError;
            ++m_specReadSoFarNum; // <-- go to the next spectum
            return 0;
        }
    }

    ++m_specReadSoFarNum;

    // set the start and stop time of the measurement
    if (this->m_stopTime < spec.m_info.m_stopTime)
        this->m_stopTime = spec.m_info.m_stopTime;
    if (spec.m_info.m_startTime < this->m_startTime)
        this->m_startTime = spec.m_info.m_startTime;

    // Extract the spectrometer-model from the serial-number of the spectrometer
    spec.m_info.m_specModelName = CSpectrometerDatabase::GetInstance().GuessModelFromSerial(spec.m_info.m_device).modelName;

    return 1;
}

int CScanFileHandler::GetSpectrum(CSpectrum &spec, long specNo)
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

    // set the start and stop time of the measurement
    if (this->m_stopTime < spec.m_info.m_stopTime)
        this->m_stopTime = spec.m_info.m_stopTime;
    if (spec.m_info.m_startTime < this->m_startTime)
        this->m_startTime = spec.m_info.m_startTime;

    return 1;
}
/** Gets the dark spectrum of the scan */
int CScanFileHandler::GetDark(CSpectrum &spec) const {
    spec = m_dark;

    if (m_fHasDark)
        return 0;
    else
        return 1;
}

/** Gets the sky spectrum of the scan */
int CScanFileHandler::GetSky(CSpectrum &spec) const {
    spec = m_sky;

    if (m_fHasSky)
        return 0;
    else
        return 1;
}

/** Gets the offset spectrum of the scan - if any */
int CScanFileHandler::GetOffset(CSpectrum &spec) const {
    spec = m_offset;

    if (m_fHasOffset)
        return 0;
    else
        return 1;
}

/** Gets the dark-current spectrum of the scan - if any */
int CScanFileHandler::GetDarkCurrent(CSpectrum &spec) const {
    spec = m_darkCurrent;

    if (m_fHasDarkCurrent)
        return 0;
    else
        return 1;
}

/** Retrieves GPS-information from the spectrum files */
const CGPSData &CScanFileHandler::GetGPS() const {
    return m_dark.GPS();
}

/** Retrieves compass-information from the spectrum files */
double CScanFileHandler::GetCompass() const {
    return m_dark.Compass();
}


/** Resets the m_specReadSoFarNum to 0 */
void  CScanFileHandler::ResetCounter() {
    m_specReadSoFarNum = 0;

    if (m_sky.ScanIndex() == 0)
        m_specReadSoFarNum = 1;

    if ((unsigned int)m_dark.ScanIndex() == m_specReadSoFarNum)
        m_specReadSoFarNum += 1;

    if ((unsigned int)m_offset.ScanIndex() == m_specReadSoFarNum || (unsigned int)m_darkCurrent.ScanIndex() == m_specReadSoFarNum)
        m_specReadSoFarNum += 1;

    if ((unsigned int)m_offset.ScanIndex() == m_specReadSoFarNum || (unsigned int)m_darkCurrent.ScanIndex() == m_specReadSoFarNum)
        m_specReadSoFarNum += 1;
}

int CScanFileHandler::GetSpectrumNumInFile() const 
{
    return 	m_specNum;
}
/** Returns the interlace steps for the spectra in this scan-file.
            @return the interlace steps for the spectra in this scan.
            @return -1 if the function 'CheckScanFile' has not been called */
int	CScanFileHandler::GetInterlaceSteps() const {
    if (!m_initialized)
        return -1;

    return this->m_sky.m_info.m_interlaceStep;
}

/** Returns the length of the spectra in this scan-file.
            @return the spectrum-length for the spectra in this scan.
            @return -1 if the function 'CheckScanFile' has not been called */
int	CScanFileHandler::GetSpectrumLength() const {
    if (!m_initialized)
        return -1;

    return this->m_sky.m_length;
}

/** Returns the start-channel for the spectra in this scan-file.
        This is the pixel on the detector for which corresponds to the first
            datapoint in the spectra (normally 0).
            @return the start-channel for the spectra in this scan.
            @return -1 if the function 'CheckScanFile' has not been called */
int	CScanFileHandler::GetStartChannel() const {
    if (!m_initialized)
        return -1;

    return this->m_sky.m_info.m_startChannel;
}

}

#ifdef _MSC_VER
#pragma warning (pop)
#endif
