#include <SpectralEvaluation/File/MKPack.h>
#include <SpectralEvaluation/Spectra/Spectrum.h>

namespace novac
{

MKPack::MKPack()
{
}

void MKPack::SetBit(std::uint8_t* pek, long bit)
{
    std::uint16_t ut = 0x80;
    pek[bit >> 3] |= (ut >> (bit & 7));
}

void MKPack::ClearBit(std::uint8_t* pek, long bit)
{
    std::uint8_t ut = 0x80;
    pek[bit >> 3] &= ~(ut >> (bit & 7));
}

void MKPack::WriteBits(std::int16_t a, std::int16_t curr, long* inpek, std::uint8_t* utpek, long bitnr)
{
    std::int16_t jj, j;
    long utwrd;
    long kk;
    std::uint16_t utlng = 0x80;

    utwrd = (std::uint32_t)((a << 5) | (curr & 0x1f));
    kk = 1L << (headsiz - 1);
    for (j = 0; j < headsiz; j++)
    {
        if (utwrd & kk) utpek[(bitnr >> 3)] |= (utlng >> (bitnr & 7));
        bitnr++;
        kk = kk >> 1;
    }

    for (jj = 0; jj < a; jj++)						 /* spara undan alla */
    {
        kk = (1L << (curr - 1));
        utwrd = *inpek++;
        for (j = 0; j < curr; j++)
        {
            if (utwrd & kk) utpek[(bitnr >> 3)] |= (utlng >> (bitnr & 7));
            bitnr++;
            kk = kk >> 1;
        }
    }
}

std::int16_t MKPack::BitsPrec(long i)
{
    std::int16_t j = 1;
    if (!i)
        return(0);

    if (i < 0) {
        if (i == -1)
            return(1);
        while (i != -1) {
            j++;
            i = (i >> 1);
        }
    }
    else {
        while (i) {
            j++;
            i = (i >> 1);
        }
    }
    return(j);
}

void MKPack::PackSeg(std::uint8_t* utpek, long* kvar)
{
    std::int16_t len[33];
    long j;
    long* incpy;
    std::int16_t curr, i, a;

    for (j = 0; j < 33; j++)
        len[j] = 0;

    incpy = m_strt;

    i = BitsPrec(*incpy++);
    curr = i;
    a = 0;
    do {
        a++;
        for (j = 0; j < curr; j++) {
            if (i > j)
                len[j] = 0;
            else {
                len[j]++;
                if (len[j] * (curr - j) > headsiz * 2) {
                    a -= len[j];
                    goto Fixat;
                }
            }
        }
        i = BitsPrec(*incpy++);
        if (i > curr)
        {
            /* i har blivit för stort. Vi ska då titta bakåt så att
                 vi tjänar in plats bakåt också på att öka bitantalet */
            if (a * (i - curr) > headsiz) goto Fixat;

            /* gå till fixat om det inte lönar sig att öka
                 bitantalet på den föregående gruppen */

            while (curr != i) /* det lönade sig att byta */
            {
                len[curr] = a;
                curr++; /* öka bitantalet */
            }
        }
    } while (a < *kvar && a < 127);
Fixat:

    WriteBits(a, curr, m_strt, utpek, m_bitnr);
    *kvar -= a;
    m_strt += a; /* öka m_strt */
    m_bitnr += a * curr + headsiz;
}

std::uint16_t MKPack::mk_compress(long* in, std::uint8_t* ut, std::uint16_t size)
{
    std::uint16_t outsize;

    m_strt = in;
    long kvar = size;

    m_bitnr = 0;
    do {
        PackSeg(ut, &kvar);
    } while (kvar > 0);

    outsize = (std::uint16_t)(m_bitnr + 7) >> 3;
    return(outsize);
}

long MKPack::UnPack(std::uint8_t* inpek, long kvar, long* ut)
{
    long* utpek = NULL;
    std::int16_t len, curr;
    std::int16_t j, jj;
    long a;
    std::uint16_t lentofile = 0;
    long bit = 0;

    // validate the input data - Added 2006.02.13 by MJ
    if (kvar > MAX_SPECTRUM_LENGTH)
        return -1;
    if (ut == NULL || inpek == NULL)
        return -1;


    utpek = ut;
    lentofile = 0;
    while (kvar > 0)
    {
        len = 0;
        for (j = 0; j < 7; j++)
        {
            len += len;
            len |= inpek[(bit >> 3)] >> (7 - (bit & 0x7)) & 1;
            bit++;
        }
        curr = 0;
        for (j = 0; j < 5; j++)
        {
            curr += curr;
            curr |= inpek[(bit >> 3)] >> (7 - (bit & 0x7)) & 1;
            bit++;
        }
        if (curr)
        {
            for (jj = 0; jj < len; jj++)
            {
                a = inpek[(bit >> 3)] >> (7 - (bit & 0x7)) & 1;
                if (a)
                    a = -1;
                bit++;
                for (j = 1; j < curr; j++)
                {
                    a += a;
                    a |= inpek[(bit >> 3)] >> (7 - (bit & 0x7)) & 1;
                    bit++;
                }
                *utpek++ = a;
            }
        }
        else
            for (jj = 0; jj < len; jj++)
                *utpek++ = 0;

        kvar -= len;
        lentofile += len;
    }
    for (jj = 1; jj < lentofile; jj++)
    {
        ut[jj] += ut[jj - 1];
    }

    return(lentofile);
}
}

