#pragma once

#ifndef DATETIME_H
#define DATETIME_H

class CDateTime
{
public:
    CDateTime() = default;
    ~CDateTime() = default;

    CDateTime(const CDateTime &t2);
    CDateTime(int y, int mo, int d, int h, int mi, int s);

    // ---------------- DATA ------------------
    /** The year (0 - 16384) */
    unsigned short year = 0;

    /** The month (1-12) */
    unsigned char month = 0;

    /** The day (1-31) */
    unsigned char day = 0;

    /** The hour (0-23) */
    unsigned char hour = 0;

    /** The minute (0-59) */
    unsigned char minute = 0;

    /** The second (0-59) */
    unsigned char second = 0;

    /** The milli second (0-999) */
    unsigned short millisecond = 0;

    // --------------- OPERATORS ------------------
    bool operator<(const CDateTime& t2) const;
    bool operator==(const CDateTime& t2) const;
    CDateTime &operator=(const CDateTime& t2);

    // --------------- METHODS --------------------

    /** Sets the time of this CDateTime-object to the current (local PC) time */
    void SetToNow();

    /** Increments the current time with the supplied number of seconds */
    void Increment(int seconds);

    /** Calculates the difference, in seconds, between two times.
            If t2 is later than t1, then the result will be negative. */
    static double Difference(const CDateTime &t1, const CDateTime &t2);

};


#endif