#pragma once

#include <cstdint>

namespace novac
{
#define hdr_version 5

    /* The structure used for saving spectra to disk.
        The types were changed on 2018-01-13 by Mattias from 'unsigned short', 'unsigned long' etc
        to 'std::uint16_t', 'std::uint32_t' etc to better handle different types of architechtures.
        This header is written to/read from file as it is and it is very important that the layout
        of the structure doesn't change */
    struct MKZYhdr
    {
        char ident[4];                 // "MKZY"
        std::uint16_t hdrsize;         // this is the size in bytes of the header
        std::uint16_t hdrversion;      // version of the header
        std::uint16_t size;            // the number of bytes with compressed data
        std::uint16_t checksum;        // checksum for the uncompressed data
        char name[12];                 // the name of this specific measurement
        char instrumentname[16];       // the name of the instrument
        std::uint16_t startc;          // the startchannel for the first data-point
        std::uint16_t pixels;          // number of pixels saved in the data-field
        std::int16_t viewangle;        // the viewing angle of the instrument
        std::uint16_t scans;           // total number of scans added
        std::int16_t exptime;          // exposure time, negative if set automatic
        std::uint8_t channel;          // channel of the spectrometer, typically 0
        std::uint8_t flag;             // for further use, currently contains the
                                       // status of the solenoid(s) in bit 0 and 1
        std::uint32_t date;            // date.
        std::uint32_t starttime;       // time when the scanning was started.
        std::uint32_t stoptime;        // time when the scanning was finished.
        double lat;                    // GPS latitude in degrees
        double lon;                    // GPS longitude in degrees
        std::int16_t altitude;         // new in version 2
        char measureidx;               // new in version 2, nr between 0 and measurecnt-1
        char measurecnt;               //new in version 2
                                       // number of MEAS= lines in cfg.txt
        std::int16_t viewangle2;       //new in version 3, direction of 2nd motor
        std::int16_t compassdir;       //new in version 3, given in cfg.txt
        std::int16_t tiltX;            //new in version 3, given in cfg.txt
        std::int16_t tiltY;            //new in version 3, given in cfg.txt
        float temperature;             //new in version 3, given in cfg.txt
        char coneangle;                //new in version 4, given in cfg.txt
        std::uint16_t ADC[8];          //new in version 5
    };

#define headsiz 12

    /** <b>MKPack</b> is a class for reading and writing spectra in
        Manne Kihlman's Pak-format. */
    class MKPack
    {
    public:
        MKPack();
        ~MKPack();

        /** Compress the given input-buffer to the given output-buffer
            @param in - the uncompressed spectral data
            @param ut - will on return contain the compressed spectrum
            @param size - the number of data-points in the uncompressed spectrum
            @return - the number of data-points in the compressed spectrum */
        std::uint16_t mk_compress(long *in, std::uint8_t *ut, std::uint16_t size);

        /** Uncompress the given input-buffer to the given output-buffer
                @param inpek - the input-buffer
                @param kvar - the number of data-points to be written to the output-buffer
                @param ut - will in return contain the uncompressed spectrum
                @return -  the number of data-points in the uncompressed spectrum ???? */
        long UnPack(std::uint8_t *inpek, long kvar, long *ut);

    private:
        void SetBit(std::uint8_t *pek, long bit);
        void ClearBit(std::uint8_t *pek, long bit);

        void WriteBits(std::int16_t a, std::int16_t curr, long *inpek, std::uint8_t *utpek, long bitnr);
        std::int16_t BitsPrec(long i);

        void PackSeg(std::uint8_t *utpek, long *kvar);

        long m_bitnr;
        long *m_strt;
    };
}