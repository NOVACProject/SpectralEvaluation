#include <SpectralEvaluation/DateTime.h>
#include <time.h>
#include <cstring>
#include <string>
#include <SpectralEvaluation/StringUtils.h>

namespace novac
{

    CDateTime::CDateTime(const CDateTime& t2)
        : year(t2.year), month(t2.month), day(t2.day), hour(t2.hour), minute(t2.minute), second(t2.second), millisecond(t2.millisecond)
    {
    }

    CDateTime::CDateTime(int y, int mo, int d, int h, int mi, int s)
        : year((unsigned short)y),
        month((unsigned char)mo),
        day((unsigned char)d),
        hour((unsigned char)h),
        minute((unsigned char)mi),
        second((unsigned char)s),
        millisecond(0)
    {
    }

    bool CDateTime::operator <(const CDateTime& t2) const
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

    bool CDateTime::operator <=(const CDateTime& t2) const {
        if (*this < t2)
            return true;
        if (*this == t2)
            return true;

        return false;
    }

    bool CDateTime::operator >(const CDateTime& t2) const {
        if (*this < t2)
            return false;
        if (*this == t2)
            return false;

        return true;
    }

    bool CDateTime::operator >=(const CDateTime& t2) const {
        if (*this < t2)
            return false;

        return true;
    }

    CDateTime& CDateTime::operator=(const CDateTime& t2) {
        this->year = t2.year;
        this->month = t2.month;
        this->day = t2.day;
        this->hour = t2.hour;
        this->minute = t2.minute;
        this->second = t2.second;
        return *this;
    }

    void CDateTime::SetToNow() {
        struct tm* tim;
        time_t t;
        time(&t);
        tim = localtime(&t);

        this->year = (unsigned short)(1900 + tim->tm_year);
        this->month = (unsigned char)(1 + tim->tm_mon);
        this->day = (unsigned char)tim->tm_mday;
        this->hour = (unsigned char)tim->tm_hour;
        this->minute = (unsigned char)tim->tm_min;
        this->second = (unsigned char)tim->tm_sec;
    }

    void CDateTime::SetToNowUTC() {
        struct tm* tim;
        time_t t;
        time(&t);
        tim = gmtime(&t);

        this->year = (unsigned short)(1900 + tim->tm_year);
        this->month = (unsigned char)(1 + tim->tm_mon);
        this->day = (unsigned char)tim->tm_mday;
        this->hour = (unsigned char)tim->tm_hour;
        this->minute = (unsigned char)tim->tm_min;
        this->second = (unsigned char)tim->tm_sec;
    }

    double CDateTime::Difference(const CDateTime& t1, const CDateTime& t2) {
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
        // calculate the number of seconds since midnight this is...
        long nSecsSinceMidnight = second + 60 * minute + 3600 * hour;

        if (secs <= (86400 - nSecsSinceMidnight)) {
            // this is a small change, only change within the same day
            nSecsSinceMidnight += secs;

            this->hour = (unsigned char)(nSecsSinceMidnight / 3600);
            this->minute = (unsigned char)((nSecsSinceMidnight - 3600 * hour) / 60);
            this->second = (unsigned char)((nSecsSinceMidnight % 60));

            return;
        }
        else {
            // go to the next before
            IncrementOneDay();

            while (secs > 86400) {
                DecrementOneDay();
                secs -= 86400;
            }

            // we should increment less than 86400 seconds (less than one day)
            nSecsSinceMidnight = secs - (86400 - hour * 3600 - minute * 60 - second);

            this->hour = (unsigned char)(nSecsSinceMidnight / 3600);
            this->minute = (unsigned char)((nSecsSinceMidnight - 3600 * hour) / 60);
            this->second = (unsigned char)((nSecsSinceMidnight % 60));

            return;
        }
    }

    /** Decrements the current time with the supplied number of seconds */
    void CDateTime::Decrement(int secs) {
        // calculate the number of seconds since midnight this is...
        long nSecsSinceMidnight = second + 60 * minute + 3600 * hour;

        if (secs <= nSecsSinceMidnight) {
            // this is a small change, only change within the same day
            nSecsSinceMidnight -= secs;

            this->hour = (unsigned char)(nSecsSinceMidnight / 3600);
            this->minute = (unsigned char)((nSecsSinceMidnight - 3600 * hour) / 60);
            this->second = (unsigned char)((nSecsSinceMidnight % 60));

            return;
        }
        else {
            // go to the day before
            DecrementOneDay();

            while (secs > 86400) {
                DecrementOneDay();
                secs -= 86400;
            }

            // we should decrement less than 86400 seconds (less than one day)
            nSecsSinceMidnight = nSecsSinceMidnight + 86400 - secs;

            this->hour = (unsigned char)(nSecsSinceMidnight / 3600);
            this->minute = (unsigned char)((nSecsSinceMidnight - 3600 * hour) / 60);
            this->second = (unsigned char)((nSecsSinceMidnight % 60));

            return;
        }
    }

    /** Decrements the current time to the same time the day before */
    void CDateTime::IncrementOneDay() {
        long daysInThisMonth = DaysInMonth(this->year, this->month);

        if (day < daysInThisMonth) {
            day += 1;
            return;
        }
        else {
            if (this->month < 12) {
                day = (unsigned char)DaysInMonth(this->year, this->month);
                month += 1;
            }
            else {
                day = 01;
                month = 01;
                year += 1;
            }
        }
    }

    /** Decrements the current time to the same time the day before */
    void CDateTime::DecrementOneDay() {
        if (day > 1) {
            day -= 1;
            return;
        }
        else {
            if (this->month > 1) {
                day = (unsigned char)DaysInMonth(this->year, this->month - 1);
                month -= 1;
            }
            else {
                day = 31;
                month = 12;
                year -= 1;
            }
        }
    }

    double CDateTime::SecondsSinceMidnight() const
    {
        return 0.001 * millisecond + second + 60 * minute + 3600 * hour;
    }

    bool CDateTime::ParseDate(const char* dateStr, CDateTime& t)
    {
        char* pt = NULL, * pt2 = NULL;
        char str[256]; // local copy
        snprintf(str, 255, "%s", (const char*)dateStr);
        int y, m, d;
        y = m = d = 0;

        // look if this is a function...
        if (EqualsIgnoringCase(dateStr, "TODAY(", 6))
        {
            int numberOfDays = 0;
            t.SetToNow();
            sprintf(str, "%s", (const char*)(dateStr + 7));
            char* right = strchr(str, ')');
            if (NULL != right) {
                right[0] = '\0';
            }
            if (sscanf(str, "%d", &numberOfDays)) {
                if (numberOfDays < 0) {
                    for (int k = 0; k > numberOfDays; --k) {
                        t.DecrementOneDay();
                    }
                }
                else if (numberOfDays > 0) {
                    for (int k = 0; k < numberOfDays; ++k) {
                        t.IncrementOneDay();
                    }
                }
            }
            return true;
        }

        // guess that the year-month-day might be separated with dots
        pt = strstr(str, ".");
        if (pt != NULL) {
            pt2 = strstr(pt + 1, ".");
            if (pt2 == NULL) {
                return false;
            }
            else {
                if (3 == sscanf(str, "%d.%d.%d", &y, &m, &d)) {
                    t.year = (unsigned short)y;
                    t.month = (unsigned char)m;
                    t.day = (unsigned char)d;
                    return true;
                }
                else {
                    return false;
                }
            }
        }

        // guess that the year-month-day might be separated with colons
        pt = strstr(str, ":");
        if (pt != NULL) {
            pt2 = strstr(pt + 1, ":");
            if (pt2 == NULL) {
                return false;
            }
            else {
                if (3 == sscanf(str, "%d:%d:%d", &y, &m, &d)) {
                    t.year = (unsigned short)y;
                    t.month = (unsigned char)m;
                    t.day = (unsigned char)d;
                    return true;
                }
                else {
                    return false;
                }
            }
        }

        // guess that the year-month-day might be separated with nothing
        int N = (int)strlen(str);
        if (N != 8)
            return false;

        if (0 == sscanf(&str[N - 2], "%d", &d))
            return false;
        str[N - 2] = '\0';
        if (0 == sscanf(&str[N - 4], "%d", &m))
            return false;
        str[N - 4] = '\0';
        if (0 == sscanf(&str[0], "%d", &y))
            return false;

        t.year = (unsigned short)y;
        t.month = (unsigned char)m;
        t.day = (unsigned char)d;
        return true;
    }

    int	DayNr(const unsigned short day[3]) {
        CDateTime d;
        d.year = day[2];
        d.month = (unsigned char)day[1];
        d.day = (unsigned char)day[0];
        return novac::DayNr(d);
    }

    int	DayNr(const CDateTime& day) {
        // Check errors in input
        if (day.month <= 0 || day.month > 12 || day.day < 1 || day.day > DaysInMonth(day.year, day.month))
            return 0;

        int dayNr = day.day; // the daynumber

        int m = day.month;
        while (m > 1) {
            dayNr += DaysInMonth(day.year, m - 1);
            --m;
        }

        return dayNr;
    }

    double JulianDay(const CDateTime& utcTime) {
        int N, J;
        int b = 0;
        double Hd, JD;

        if (utcTime.year < 1901 || utcTime.year > 2100 || utcTime.month < 1 || utcTime.month > 12 || utcTime.day < 1 || utcTime.day > 31)
            return 0.0;
        if (utcTime.hour < 0 || utcTime.hour > 23 || utcTime.minute < 0 || utcTime.minute > 59 || utcTime.second < 0 || utcTime.second > 59)
            return 0.0;

        N = 4713 + utcTime.year - 1;
        J = N * 365 + N / 4 - 10 - 3;

        if (N % 4 == 1 || N % 4 == 2 || N % 4 == 3)
            ++J;

        switch (utcTime.month)
        {
        case(1):  b = utcTime.day - 1;		break;
        case(2):  b = utcTime.day + 30;	break;
        case(3):  b = utcTime.day + 58;	break;
        case(4):  b = utcTime.day + 89;	break;
        case(5):  b = utcTime.day + 119;	break;
        case(6):  b = utcTime.day + 150;	break;
        case(7):  b = utcTime.day + 180;	break;
        case(8):  b = utcTime.day + 211;	break;
        case(9):  b = utcTime.day + 242;	break;
        case(10): b = utcTime.day + 272;	break;
        case(11): b = utcTime.day + 303;	break;
        case(12): b = utcTime.day + 333;	break;
        }
        if (utcTime.year % 4 == 0)
            ++b;

        Hd = (utcTime.hour - 12.0) / 24.0 + utcTime.minute / (60.0 * 24.0) + utcTime.second / (3600.0 * 24.0);						/*CONVERSION HORAR TO DECIMAL SYSTEM*/

        JD = J + Hd + b;
        return JD;
    }

    void SplitToHourMinuteSecond(const int time, int& hours, int& minutes, int& seconds)
    {
        hours = (int)(time / 3600);
        minutes = (time - hours * 3600) / 60;
        seconds = time % 60;
    }

}
