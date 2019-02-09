#include "DateTime.h"

#include <time.h>

CDateTime::CDateTime(const CDateTime &t2) {
    year = t2.year;
    month = t2.month;
    day = t2.day;
    hour = t2.hour;
    minute = t2.minute;
    second = t2.second;
    millisecond = t2.millisecond;
}

CDateTime::CDateTime(int y, int mo, int d, int h, int mi, int s) {
    year = y;
    month = mo;
    day = d;
    hour = h;
    minute = mi;
    second = s;
    millisecond = 0;
}

bool CDateTime::operator <(const CDateTime &t2) const
{
    if (this->year < t2.year)
        return true;
    if (this->year > t2.year)
        return false;

    if (this->month < t2.month)
        return true;
    if (this->month > t2.month)
        return false;

    if (this->day < t2.day)
        return true;
    if (day > t2.day)
        return false;

    if (this->hour < t2.hour)
        return true;
    if (this->hour > t2.hour)
        return false;

    if (this->minute < t2.minute)
        return true;
    if (this->minute > t2.minute)
        return false;

    if (this->second < t2.second)
        return true;
    if (this->second > t2.second)
        return false;

    if (this->millisecond < t2.millisecond)
        return true;
    if (this->millisecond > t2.millisecond)
        return false;

    // equal
    return false;
}

bool CDateTime::operator==(const CDateTime& t2) const
{
    if (this->millisecond != t2.millisecond)
        return false;
    if (this->second != t2.second)
        return false;
    if (this->minute != t2.minute)
        return false;
    if (this->hour != t2.hour)
        return false;
    if (this->day != t2.day)
        return false;
    if (this->month != t2.month)
        return false;
    if (this->year != t2.year)
        return false;
    return true;
}

CDateTime	& CDateTime::operator=(const CDateTime& t2) {
    this->year = t2.year;
    this->month = t2.month;
    this->day = t2.day;
    this->hour = t2.hour;
    this->minute = t2.minute;
    this->second = t2.second;
    return *this;
}

void CDateTime::SetToNow() {
    struct tm *tim;
    time_t t;
    time(&t);
    tim = localtime(&t);

    this->year = 1900 + tim->tm_year;
    this->month = 1 + tim->tm_mon;
    this->day = tim->tm_mday;
    this->hour = tim->tm_hour;
    this->minute = tim->tm_min;
    this->second = tim->tm_sec;
}

/** Calculates the difference, in seconds, between two times.
        If t2 is later than t1, then the result will be negative. */
double CDateTime::Difference(const CDateTime &t1, const CDateTime &t2) {
    struct tm tid1, tid2;
    tid1.tm_year = t1.year - 1900;
    tid1.tm_mon = t1.month - 1;
    tid1.tm_mday = t1.day;
    tid1.tm_hour = t1.hour;
    tid1.tm_min = t1.minute;
    tid1.tm_sec = t1.second;
    tid1.tm_wday = 0;
    tid1.tm_yday = 0;
    tid1.tm_isdst = 0;

    tid2.tm_year = t2.year - 1900;
    tid2.tm_mon = t2.month - 1;
    tid2.tm_mday = t2.day;
    tid2.tm_hour = t2.hour;
    tid2.tm_min = t2.minute;
    tid2.tm_sec = t2.second;
    tid2.tm_wday = 0;
    tid2.tm_yday = 0;
    tid2.tm_isdst = 0;

    time_t t_1 = mktime(&tid1);
    time_t t_2 = mktime(&tid2);

    return difftime(t_1, t_2);
}

/** Helper function, Takes a given year and month and returns the number of days in that month. */
int	DaysInMonth(int year, int month)
{
    static int nDays[] = { 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 };

    // detect non-existing months.
    if (month < 1 || month > 12)
        return 0;

    // If the month is not february, then it's easy!!!
    if (month != 2)
        return nDays[month - 1];

    // If february, then check for leap-years
    if (year % 4 != 0)
        return 28; // not a leap-year

    if (year % 400 == 0) // every year dividable by 400 is a leap-year
        return 29;

    if (year % 100 == 0) // years diviable by 4 and by 100 are not leap-years
        return 28;
    else
        return 29;		// years dividable by 4 and not by 100 are leap-years
}

void CDateTime::Increment(int secs) {
    this->second += secs;
    if (this->second < 60)
        return;

    this->minute += this->second / 60;
    this->second %= 60;
    if (this->minute < 60)
        return;

    this->hour += this->minute / 60;
    this->minute %= 60;
    if (this->hour < 24)
        return;

    int daysInMonth = DaysInMonth(this->year, this->month);
    this->day += this->hour / 24;
    this->hour %= 24;
    if (this->day < daysInMonth)
        return;

    while (this->day / daysInMonth > 1) {
        this->month += 1;
        if (this->month > 12) {
            this->month = 1;
            this->year += 1;
        }
        this->day -= daysInMonth;
        daysInMonth = DaysInMonth(this->year, this->month);
    }
}
