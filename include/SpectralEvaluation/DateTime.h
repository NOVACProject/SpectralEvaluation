#pragma once

#ifndef DATETIME_H
#define DATETIME_H

#include <iosfwd>
#include <string>

namespace novac
{
class CDateTime
{
public:
    CDateTime() = default;
    ~CDateTime() = default;

    CDateTime(int y, int mo, int d, int h, int mi, int s);
    CDateTime(int y, int mo, int d, int h, int mi, int s, int ms);

    // Object is copyable
    CDateTime(const CDateTime& t2);
    CDateTime& operator=(const CDateTime& t2);

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
    bool operator<=(const CDateTime& t2) const;
    bool operator>(const CDateTime& t2) const;
    bool operator>=(const CDateTime& t2) const;

    friend std::ostream& operator<<(std::ostream& os, const CDateTime& dt);

    std::string ToString() const;

    // --------------- METHODS --------------------

    /** Sets the time of this CDateTime-object to the current (local PC) time */
    void SetToNow();

    /** Sets the time of this CDateTime-object to the current (UTC) time */
    void SetToNowUTC();

    /** Increments the current time with the supplied number of seconds */
    void Increment(int seconds);

    /** Decrements the current time with the supplied number of seconds */
    void Decrement(int seconds);

    /** Increments the current time to the same time the day before */
    void IncrementOneDay();

    /** Decrements the current time to the same time the day before */
    void DecrementOneDay();

    /** Calculates the difference, in seconds, between two times.
            If t2 is later than t1, then the result will be negative.
            If either of the two is an invalid time (e.g. month=0 or month=14) then zero is returned. */
    static double Difference(const CDateTime& t1, const CDateTime& t2);

    /** Attempts to parse the date found in the string 'dateStr', the resulting
        date is (on successful parsing) returned in 't'

        - The date MUST be in either of the formats YYYY.MM.DD / YYYY:MM:DD or YYYYMMDD

        - This can also manage functional expressions such as;
            TODAY(n) where n is an integer. This sets t to be n-days after today (before if n is negative)

        NOTE only the date is parsed not the time!!!!
        @return true on success else false.
    */
    static bool ParseDate(const char* dateStr, CDateTime& t);

    /** Calculates and returns the number of seconds since midnight of this timestamp.
    *   Ignores the date component and will only use the hours, minutes, seconds and milliseconds */
    double SecondsSinceMidnight() const;
};


/** Takes a given year and month and returns the number of days in that month.
    The month ranges from 1 to 12. Any illegal values in the month will return 0. */
int DaysInMonth(int year, int month);

/** Takes a given date and calculates the day of the year.
    An illegal day will return 0. */
int DayNr(const unsigned short day[3]);

/** Takes a given date and calculates the day of the year.
    An illegal day will return 0. */
int DayNr(const CDateTime& day);

/** Returns the Julian Day corresponding to the given date and time of day. */
double JulianDay(const CDateTime& utcTime);

/** Splits a time, given in seconds since midnight to hour, minutes and seconds */
void SplitToHourMinuteSecond(const int time, int& hours, int& minutes, int& seconds);

}

#endif
