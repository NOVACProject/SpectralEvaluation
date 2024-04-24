#include <SpectralEvaluation/File/SpectrumIO.h>
#include <SpectralEvaluation/File/MKPack.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>
#include <SpectralEvaluation/Spectra/SpectrometerModel.h>
#include <SpectralEvaluation/StringUtils.h>
#include <cstring>
#include <algorithm>


#ifdef _MSC_VER
#pragma warning (push, 4)
#endif

#undef min
#undef max

namespace novac
{
    CSpectrumIO::CSpectrumIO()
    {
        this->m_lastError = ERROR_NO_ERROR;
        this->m_buffer.resize(16384);
        this->m_outbuf.resize(16384);
    }

    int CSpectrumIO::CountSpectra(const std::string& fileName)
    {
        std::uint32_t specNum = 0;

        FILE* f = fopen(fileName.c_str(), "rb");

        if (f == nullptr)
        {
            printf("Could not open spectrum file: %s\n", fileName.c_str());
            m_lastError = ERROR_COULD_NOT_OPEN_FILE;
            return(1);
        }

        while (1)
        {
            int headerSize;
            struct MKZYhdr MKZY;
            int ret = ReadNextSpectrumHeader(f, MKZY, headerSize);
            if (ret == 1)
            {
                break;
            }
            if (ret == 2)
            {
                continue;
            }

            // Seek our way into the next spectrum...
            if (0 != fseek(f, std::min(MKZY.size, (std::uint16_t)(4 * MAX_SPECTRUM_LENGTH)), SEEK_CUR))
            {
                break;
            }

            // Make sure we're at the right place, if not rewind again and search for the next
            // occurence of the "MKZY" string, which signals the start of a 'new' spectrum.
            char textBuffer[4];
            fread(textBuffer, 1, 4, f);
            if (nullptr == strstr(textBuffer, "MKZY"))
            {
                // rewind
                if (0 != fseek(f, -std::min(MKZY.size, (std::uint16_t)(4 * MAX_SPECTRUM_LENGTH)), SEEK_CUR))
                {
                    break;
                }
            }
            else
            {
                if (0 != fseek(f, -4, SEEK_CUR))
                {
                    break;
                }
            }

            ++specNum;
        }

        // signals that there's no error in the file.
        m_lastError = ERROR_SPECTRUM_NOT_FOUND;

        fclose(f);

        return specNum;
    }

    std::uint32_t CSpectrumIO::ScanSpectrumFile(const std::string& fileName, const std::string* specNamesToLookFor, int numSpecNames, int* indices)
    {
        std::uint32_t specNum = 0;

        FILE* f = fopen(fileName.c_str(), "rb");

        if (f == nullptr)
        {
            printf("Could not open spectrum file: %s\n", fileName.c_str());
            m_lastError = ERROR_COULD_NOT_OPEN_FILE;
            return 0;
        }

        while (1)
        {
            int headerSize;
            struct MKZYhdr MKZY;
            int ret = ReadNextSpectrumHeader(f, MKZY, headerSize);

            if (ret == 1)
            {
                break;
            }
            if (ret == 2)
            {
                continue;
            }

            /** Look in the buffer */
            // 1. Clean the spectrum name from special characters...
            std::string specName{ MKZY.name };
            CleanString(specName);
            Trim(specName, " \t");
            const size_t size1 = specName.size();

            for (int nameIndex = 0; nameIndex < numSpecNames; ++nameIndex)
            {
                // first of all, the strings must have equal size to be equal...
                if (size1 != specNamesToLookFor[nameIndex].size())
                {
                    continue;
                }

                if (EqualsIgnoringCase(specNamesToLookFor[nameIndex], specName))
                {
                    indices[nameIndex] = specNum;
                    continue;
                }
            }

            char textBuffer[4];

            // Seek our way into the next spectrum...
            if (0 != fseek(f, std::min(MKZY.size, (std::uint16_t)(4 * MAX_SPECTRUM_LENGTH)), SEEK_CUR))
                break;

            // Make sure we're at the right place, if not rewind again and search for the next
            // occurence of the "MKZY" string, which signals the start of a 'new' spectrum.
            fread(textBuffer, 1, 4, f);
            if (nullptr == strstr(textBuffer, "MKZY")) {
                // rewind
                if (0 != fseek(f, -std::min(MKZY.size, (std::uint16_t)(4 * MAX_SPECTRUM_LENGTH)), SEEK_CUR))
                {
                    break;
                }
            }
            else
            {
                if (0 != fseek(f, -4, SEEK_CUR))
                {
                    break;
                }
            }

            ++specNum;
        }

        // signals that there's no error in the file.
        this->m_lastError = ERROR_SPECTRUM_NOT_FOUND;

        fclose(f);

        return specNum;
    }

    bool CSpectrumIO::ReadSpectrum(const std::string& fileName, const int spectrumNumber, CSpectrum& spec, char* headerBuffer /* = nullptr*/, int headerBufferSize /* = 0*/, int* headerSize /* = nullptr*/)
    {
        long currentSpectrumNumber = 0;

        FILE* f = fopen(fileName.c_str(), "rb");

        if (f == nullptr)
        {
            printf("Could not open spectrum file: %s\n", fileName.c_str());
            m_lastError = ERROR_COULD_NOT_OPEN_FILE;
            return false;
        }

        while (1)
        {
            int ret;
            struct MKZYhdr MKZY;
            if (headerBuffer != nullptr)
            {
                ret = ReadNextSpectrumHeader(f, MKZY, *headerSize, &spec, headerBuffer, headerBufferSize);
            }
            else
            {
                int hdrSize;
                ret = ReadNextSpectrumHeader(f, MKZY, hdrSize, &spec);
            }
            if (ret == 1)
            {
                break;
            }
            if (ret == 2)
            {
                continue;
            }

            if (currentSpectrumNumber != spectrumNumber)
            {
                char textBuffer[4];

                // Seek our way into the next spectrum...
                if (0 != fseek(f, std::min(MKZY.size, (std::uint16_t)(4 * MAX_SPECTRUM_LENGTH)), SEEK_CUR))
                {
                    break;
                }

                // Make sure we're at the right place, if not rewind again and search for the next
                // occurence of the "MKZY" string, which signals the start of a 'new' spectrum.
                if (fread(textBuffer, 1, 4, f) < 4)
                {
                    break;
                }
                if (nullptr == strstr(textBuffer, "MKZY"))
                {
                    // rewind
                    if (0 != fseek(f, -std::min(MKZY.size, (std::uint16_t)(4 * MAX_SPECTRUM_LENGTH)), SEEK_CUR))
                    {
                        break;
                    }
                }
                else
                {
                    if (0 != fseek(f, -4, SEEK_CUR))
                    {
                        break;
                    }
                }

                ++currentSpectrumNumber;
                continue;
            }
            else
            {
                // read the spectrum from the file
                if (MKZY.size > m_buffer.size())
                {
                    // compressed data is too long. We cannot read the full spectrum.
                    m_lastError = ERROR_SPECTRUM_TOO_LARGE;
                    fclose(f);
                    return false;
                }

                if (fread(m_buffer.data(), 1, MKZY.size, f) < MKZY.size) //read compressed info
                {
                    printf("Error EOF! in %s\n", fileName.c_str());
                    fclose(f);
                    m_lastError = ERROR_EOF;
                    return false;
                }

                if (MKZY.pixels > m_buffer.size())
                {
                    // The spectrum is longer than what the buffer can handle. Trying to
                    // uncompress the whole spectrum will result in a buffer overflow.
                    // this spectrum cannot be read - return.
                    m_lastError = ERROR_SPECTRUM_TOO_LARGE;
                    fclose(f);
                    return false;
                }

                MKPack mkPack;
                const long outlen = mkPack.UnPack(m_buffer.data(), MKZY.pixels, m_outbuf.data()); //uncompress info(compressed buffer,num of sampling points, uncompressedinfo)

                // validate that the decompression was ok - Added 2006.02.13 by MJ
                if (outlen < 0)
                {
                    m_lastError = ERROR_DECOMPRESS;
                    fclose(f);
                    return false;
                }

                // validate that the spectrum is not too large - Added 2006.02.13 by MJ
                if (outlen > MAX_SPECTRUM_LENGTH)
                {
                    m_lastError = ERROR_SPECTRUM_TOO_LARGE;
                    fclose(f);
                    return false;
                }

                // calculate the checksum
                std::uint32_t chk = 0;
                for (long j = 0; j < outlen && j < MAX_SPECTRUM_LENGTH; j++)
                {
                    chk += m_outbuf[j];
                    spec.m_data[j] = m_outbuf[j];
                }
                const std::uint16_t* p = (std::uint16_t*)&chk;
                const std::uint16_t checksum = p[0] + p[1];

                if (checksum != MKZY.checksum)
                {
                    printf("Checksum mismatch %04x!=x%04x\n", checksum, MKZY.checksum);

                    m_lastError = ERROR_CHECKSUM_MISMATCH;
                    fclose(f);
                    return false;
                }

                // Get the maximum intensity
                if (MKZY.pixels > 0)
                {
                    spec.m_info.m_peakIntensity = (float)spec.MaxValue();
                    spec.m_info.m_offset = (float)spec.GetOffset();
                }

                fclose(f);

                return true;
            }
        }
        fclose(f);

        this->m_lastError = ERROR_SPECTRUM_NOT_FOUND;
        return false; // spectrum not found
    }

    bool CSpectrumIO::FindSpectrumNumber(FILE* f, int spectrumNumber)
    {
        std::string errorMessage; // a string used for error messages
        long curSpecNum = 0;

        if (f == nullptr)
        {
            return false;
        }

        // first rewind the file
        rewind(f);

        while (curSpecNum <= spectrumNumber)
        {
            int c;

            // find the next 'MKZY' - string in the spectrum file
            while ((c = getc(f)) != (int)'M')
            {
                if (c == EOF)
                {
                    return false;
                }
            }

            if (getc(f) != (int)'K')
                continue;

            if (getc(f) != (int)'Z')
                continue;

            if (getc(f) != (int)'Y')
                continue;

            // we've found a 'new' 'MKZY'-string, call this a spectrum and increase
            //	the spectrum counter...
            ++curSpecNum;
        }

        // we've found the spectrum we're looking for, now rewind past the 'MKZY'-string
        if (0 != fseek(f, -4, SEEK_CUR))
            return false;

        // signals that there's no error in the file.
        m_lastError = ERROR_NO_ERROR;

        return true;
    }

    bool CSpectrumIO::ReadNextSpectrum(FILE* f, CSpectrum& spec)
    {
        int tmp;
        return ReadNextSpectrum(f, spec, tmp);
    }

    bool CSpectrumIO::ReadNextSpectrum(FILE* f, CSpectrum& spec, int& headerSize, char* headerBuffer, int headerBufferSize)
    {
        std::uint16_t* p = nullptr;

        struct MKZYhdr MKZY;
        int ret = ReadNextSpectrumHeader(f, MKZY, headerSize, &spec, headerBuffer, headerBufferSize);
        if (ret != 0)
        {
            return false;
        }

        // read the spectrum from the file
        if (MKZY.size > m_buffer.size())
        {
            // compressed data is too long. We cannot read the full spectrum.
            this->m_lastError = ERROR_SPECTRUM_TOO_LARGE;
            return false;
        }

        if (fread(m_buffer.data(), 1, MKZY.size, f) < MKZY.size) //read compressed info
        {
            printf("Error EOF! in pak-file\n");

            m_lastError = ERROR_EOF;
            return false;
        }

        if (MKZY.pixels > m_outbuf.size())
        {
            // The spectrum is longer than what the buffer can handle. Trying to
            // uncompress the whole spectrum will result in a buffer overflow.
            // this spectrum cannot be read - return.
            m_lastError = ERROR_SPECTRUM_TOO_LARGE;
            return false;
        }

        // We've managed to read the spectrum header, write that information
        //	to the supplied spectrum data-structure
        spec.m_info.m_device = std::string(MKZY.instrumentname);
        Trim(spec.m_info.m_device, " ");  // remove spaces in the beginning or the end
        spec.m_info.m_name = std::string(MKZY.name);

        // Decompress the spectrum itself
        MKPack mkPack;
        const long outlen = mkPack.UnPack(m_buffer.data(), MKZY.pixels, m_outbuf.data()); //uncompress info(compressed buffer,num of sampling points, uncompressedinfo)

        // validate that the decompression was ok - Added 2006.02.13 by MJ
        if (outlen < 0)
        {
            this->m_lastError = ERROR_DECOMPRESS;
            return false;
        }
        // validate that the spectrum is not too large - Added 2006.02.13 by MJ
        if (outlen > MAX_SPECTRUM_LENGTH)
        {
            this->m_lastError = ERROR_SPECTRUM_TOO_LARGE;
            return false;
        }

        // calculate the checksum
        std::uint32_t chk = 0;
        for (long j = 0; j < outlen && j < MAX_SPECTRUM_LENGTH; j++)
        {
            chk += m_outbuf[j];
        }
        p = (std::uint16_t*)&chk;

        const std::uint16_t checksum = p[0] + p[1];
        if (checksum != MKZY.checksum)
        {
            printf("Checksum mismatch %04x!=x%04x\n", checksum, MKZY.checksum);

            this->m_lastError = ERROR_CHECKSUM_MISMATCH;
            return false;
        }

        // copy the spectrum
        for (long j = 0; j < outlen && j < MAX_SPECTRUM_LENGTH; j++)
        {
            spec.m_data[j] = m_outbuf[j];
        }

        // Get the maximum intensity
        if (MKZY.pixels > 0)
        {
            spec.m_info.m_peakIntensity = (float)spec.MaxValue();
            spec.m_info.m_offset = (float)spec.GetOffset();
        }

        return true;
    }

    void CSpectrumIO::ParseTime(const std::uint32_t t, CDateTime& time) const
    {
        time.hour = (unsigned char)(t / 1000000);
        time.minute = (unsigned char)((t - time.hour * 1000000) / 10000);
        time.second = (unsigned char)((t - time.hour * 1000000 - time.minute * 10000) / 100);
        time.millisecond = 10 * ((std::uint16_t)(t % 100));
    }

    void CSpectrumIO::WriteTime(std::uint32_t& t, const CDateTime& time) const {
        t = time.hour * 1000000 + time.minute * 10000 + time.second * 100 + time.millisecond / 10;
    }

    void CSpectrumIO::ParseDate(const std::uint32_t d, CDateTime& day) const {
        day.day = (unsigned char)(d / 10000);                  // the day
        day.month = (unsigned char)((d - day.day * 10000) / 100);  // the month
        day.year = (std::uint16_t)(d % 100);                  // the year

        if (day.year < 100)
            day.year += 2000; // assume the 21:st century (should be ok for another 95 years)
    }

    // Write the date in Manne's format: ddmmyy 
    void CSpectrumIO::WriteDate(std::uint32_t& d, const CDateTime& day) const {
        if (day.year < 100)
            d = day.day * 10000 + day.month * 100 + day.year;
        else
            d = day.day * 10000 + day.month * 100 + day.year - (day.year / 100) * 100;
    }

    int CSpectrumIO::AddSpectrumToFile(const std::string& fileName, const CSpectrum& spectrum, const char* headerBuffer, int headerSize, bool overwrite)
    {
        // Test the input-data
        if (spectrum.m_length <= 0)
        {
            return 1;
        }

        // ---- start by converting the spectrum into 'long'
        std::vector<long> spec(spectrum.m_length);
        for (long i = 0; i < spectrum.m_length; ++i)
        {
            spec[i] = (long)spectrum.m_data[i];
        }

        // ---- create the proper header information ---- 
        struct MKZYhdr MKZY;

        // calculate checksum
        std::uint32_t checksum = 0;
        for (long i = 0; i < spectrum.m_length; ++i)
        {
            checksum += spec[i];
        }
        const std::uint16_t* p = (std::uint16_t*)&checksum;
        MKZY.checksum = p[0] + p[1];

        // the spectrum should be stored as just the difference
        //  between each two pixels (delta compression)
        //  except for the first pixel (of course)
        long last = spec[0];
        for (long i = 1; i < spectrum.m_length; i++)
        {
            long tmp = spec[i];
            spec[i] = tmp - last;
            last = tmp;
        }

        // Compress the spectrum..
        std::vector<std::uint16_t> sbuf(16384);
        memset(sbuf.data(), 0, 16384);
        MKPack mkPack;
        const std::uint16_t outsiz = mkPack.mk_compress(spec.data(), (unsigned char*)sbuf.data(), (std::uint16_t)spectrum.m_length);
        const CSpectrumInfo& info = spectrum.m_info;

        MKZY.ident[0] = 'M';
        MKZY.ident[1] = 'K';
        MKZY.ident[2] = 'Z';
        MKZY.ident[3] = 'Y';

        MKZY.altitude = (short)spectrum.Altitude();
        MKZY.channel = info.m_channel;
        MKZY.compassdir = (short)(info.m_compass * 10.0f);

        MKZY.ADC[0] = (std::uint16_t)(info.m_batteryVoltage * 100.0f);

        MKZY.coneangle = (char)info.m_coneAngle;
        WriteDate(MKZY.date, info.m_startTime);
        MKZY.exptime = (short)info.m_exposureTime;
        MKZY.flag = info.m_flag;
        MKZY.hdrsize = sizeof(struct MKZYhdr);
        MKZY.hdrversion = hdr_version;
        sprintf(MKZY.instrumentname, "%.16s", spectrum.m_info.m_device.c_str());
        MKZY.lat = spectrum.Latitude();
        MKZY.lon = spectrum.Longitude();
        MKZY.measurecnt = (char)info.m_scanSpecNum;
        MKZY.measureidx = (char)info.m_scanIndex;
        sprintf(MKZY.name, "%.12s", spectrum.m_info.m_name.c_str());
        MKZY.pixels = (std::uint16_t)spectrum.m_length;
        MKZY.size = outsiz;
        MKZY.startc = info.m_startChannel;
        MKZY.scans = (std::uint16_t)info.m_numSpec;
        WriteTime(MKZY.starttime, info.m_startTime);
        WriteTime(MKZY.stoptime, info.m_stopTime);
        MKZY.temperature = info.m_temperature;
        MKZY.tiltX = (short)info.m_roll;		// <-- The leaning in the direction perpendicular to the scanner
        MKZY.tiltY = (short)info.m_pitch;		// <-- The leaning in the direction of the scanner
        MKZY.viewangle = (std::uint16_t)info.m_scanAngle;
        MKZY.viewangle2 = (std::uint16_t)info.m_scanAngle2;

        FILE* f = nullptr;

        if (overwrite)
        {
            f = fopen(fileName.c_str(), "wb");
        }
        else
        {
            f = fopen(fileName.c_str(), "r+b");
            if (f == nullptr) // this will happen if the file does not exist...
            {
                f = fopen(fileName.c_str(), "w+b");
            }
        }
        if (f == nullptr)
        {
            return 1;
        }

        if (0 == fseek(f, 0, SEEK_END)) {
            // Write the header
            if (headerBuffer != nullptr && headerSize != 0) {
                fwrite(headerBuffer, headerSize, 1, f);
            }
            else {
                fwrite(&MKZY, sizeof(struct MKZYhdr), 1, f);
            }

            // Write the spectrum data
            fwrite(sbuf.data(), outsiz, 1, f);
        }
        fclose(f);

        return 0;
    }

    int CSpectrumIO::ReadNextSpectrumHeader(FILE* f, struct MKZYhdr& MKZYHeader, int& headerSize, CSpectrum* spec, char* headerBuffer, int headerBufferSize)
    {
        memset(&MKZYHeader, 0, sizeof(MKZYHeader));  // clear header information
        MKZYHeader.measureidx = -1;   // this is for compatibility reasons, if the file does not contain spectrum number, we'll know about it

        if (fread(&MKZYHeader, 1, 8, f) < 8) {        // reads MKZY.ident and MKZY.hdrsize
            return 1;  // could not read header info, break /* TODO - this is not a good way to quit */
        }

        // check that we are actually at the beginning of the header - Added 2006.02.13 by MJ
        if (strncmp(MKZYHeader.ident, "MKZY", 4 * sizeof(char)))
        {
            return 2;
        }

        std::uint32_t local_headersize = MKZYHeader.hdrsize;
        headerSize = MKZYHeader.hdrsize;

        if (sizeof(MKZYHeader) < local_headersize)
        {
            local_headersize = sizeof(MKZYHeader);
        }

        /** If the file contains a smaller header than the program version can read (MKZY.hdrsize < sizeof(MKZY))
                    only read the actual header in the file.
                If the file contains a bigger header than the program version can read (sizeof(MKZY) < MKZY.hdrsize)
                    only read what we can understand. */

        if (fread((char*)&MKZYHeader + 8, 1, local_headersize - 8, f) < local_headersize - 8) // read the rest of the header
            return 1;

        // If the user wants the header in binary format, copy it
        if (headerBuffer != nullptr && headerBufferSize > (int)local_headersize) {
            memset(headerBuffer, 0, headerBufferSize);
            memcpy(headerBuffer, &MKZYHeader, local_headersize);
        }


        // Calculate how much of the information in the file that we could not read because the program version is too old.
        const long sizdiff = MKZYHeader.hdrsize - sizeof(MKZYHeader);
        if (sizdiff > 0)
        {
            // If the user want the whole header, read it. Otherwise jump formwards
            if (headerBuffer != nullptr && headerBufferSize > MKZYHeader.hdrsize)
            {
                if (fread(headerBuffer + local_headersize, 1, sizdiff, f) < (std::uint32_t)sizdiff)
                {
                    m_lastError = ERROR_SPECTRUM_NOT_FOUND;
                }
                return false;
            }
            else
            {
                fseek(f, sizdiff, SEEK_CUR);       // NOTE -- BUG CORRECTED 2006.02.14 BY MJ - was "fseek(f,sizdiff-8,SEEK_CUR);"
            }
        }

        if (spec != nullptr)
        {
            // clear the spectrum
            memset(spec->m_data, 0, MAX_SPECTRUM_LENGTH * sizeof(double));

            CSpectrumInfo* info = &spec->m_info;
            // save the spectrum information in the CSpectrum data structure
            spec->m_length = std::max(std::min(MKZYHeader.pixels, (std::uint16_t)(MAX_SPECTRUM_LENGTH)), (std::uint16_t)(0));
            info->m_startChannel = MKZYHeader.startc;
            info->m_numSpec = MKZYHeader.scans;
            info->m_exposureTime = (MKZYHeader.exptime > 0) ? MKZYHeader.exptime : -MKZYHeader.exptime;
            info->m_gps.m_longitude = MKZYHeader.lon;
            info->m_gps.m_latitude = MKZYHeader.lat;
            info->m_gps.m_altitude = MKZYHeader.altitude;
            info->m_channel = MKZYHeader.channel;
            CSpectrum::GetInterlaceSteps(info->m_channel, info->m_interlaceStep);
            info->m_scanAngle = MKZYHeader.viewangle;
            if (info->m_scanAngle > 180.0)
            {
                info->m_scanAngle -= 360.0; // map 270 -> -90
            }
            info->m_scanAngle2 = (float)MKZYHeader.viewangle2;
            info->m_coneAngle = MKZYHeader.coneangle;
            info->m_compass = (float)MKZYHeader.compassdir / 10.0f;
            if (info->m_compass > 360.0 || info->m_compass < 0)
            {
                printf("Spectrum has compass angle outside of the expected [0, 360] degree range.\n");
            }
            info->m_batteryVoltage = (float)MKZYHeader.ADC[0] / 100.0f;
            info->m_temperature = MKZYHeader.temperature;

            info->m_scanIndex = MKZYHeader.measureidx;
            info->m_scanSpecNum = MKZYHeader.measurecnt;
            info->m_flag = MKZYHeader.flag;

            ParseTime(MKZYHeader.starttime, info->m_startTime);
            ParseTime(MKZYHeader.stoptime, info->m_stopTime);
            ParseDate(MKZYHeader.date, info->m_startTime);
            ParseDate(MKZYHeader.date, info->m_stopTime);

            info->m_device = std::string(MKZYHeader.instrumentname);
            Trim(info->m_device, " "); // remove spaces in the beginning or the end
            info->m_name = std::string(MKZYHeader.name);
        }

        return 0;
    }
}

#ifdef _MSC_VER
#pragma warning (pop)
#endif
