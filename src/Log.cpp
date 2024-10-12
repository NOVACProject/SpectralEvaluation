#include <SpectralEvaluation/Log.h>
#include <iostream>
#include <sstream>

novac::LogContext novac::LogContext::With(std::string name, std::string value)
{
    novac::LogContext c(*this);
    c.properties.push_back(std::pair<std::string, std::string>(name, value));
    return c;
}

novac::LogContext novac::LogContext::With(std::string name, int value)
{
    std::stringstream s;
    s << value;
    return this->With(name, s.str());
}

novac::LogContext novac::LogContext::With(std::string name, double value)
{
    std::stringstream s;
    s << value;
    return this->With(name, s.str());
}

static std::ostream& operator << (std::ostream& out, const novac::LogContext& c)
{
    for (const auto & p : c.properties)
    {
        out << "[" << p.first << "=" << p.second << "] ";
    }

    return out;
}

static void Output(const char* level, std::string message)
{
    std::cout << level << message << std::endl;
}

static void Output(const char* level, const novac::LogContext& c, std::string message)
{
    std::cout << level << c << message << std::endl;
}

void novac::ConsoleLog::Debug(const std::string& message)
{
    Output("[Debug] ", message);
}

void novac::ConsoleLog::Debug(const LogContext& c, const std::string& message)
{
    Output("[Debug] ", c, message);
}

void novac::ConsoleLog::Information(const std::string& message)
{
    Output("[Info] ", message);
}

void novac::ConsoleLog::Information(const LogContext& c, const std::string& message)
{
    Output("[Info] ", c, message);
}

void novac::ConsoleLog::Error(const std::string& message)
{
    Output("[Error] ", message);
}

void novac::ConsoleLog::Error(const LogContext& c, const std::string& message)
{
    Output("[Error] ", c, message);
}
