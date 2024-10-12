#include <SpectralEvaluation/Log.h>
#include <iostream>
#include <sstream>

namespace novac
{

LogContext::LogContext(std::string name, std::string value)
{
    this->properties.push_back(std::pair<std::string, std::string>(name, value));
}

LogContext LogContext::With(std::string name, std::string value)
{
    LogContext c(*this);
    c.properties.push_back(std::pair<std::string, std::string>(name, value));
    return c;
}

LogContext LogContext::With(std::string name, int value)
{
    std::stringstream s;
    s << value;
    return this->With(name, s.str());
}

LogContext LogContext::With(std::string name, double value)
{
    std::stringstream s;
    s << value;
    return this->With(name, s.str());
}

std::ostream& operator << (std::ostream& out, const LogContext& c)
{
    for (const auto& p : c.properties)
    {
        out << "[" << p.first << "=" << p.second << "] ";
    }

    return out;
}

static void Output(const char* level, std::string message)
{
    std::cout << level << message << std::endl;
}

static void Output(const char* level, const LogContext& c, std::string message)
{
    std::cout << level << c << message << std::endl;
}

void ConsoleLog::Debug(const std::string& message)
{
    Output("[Debug] ", message);
}

void ConsoleLog::Debug(const LogContext& c, const std::string& message)
{
    Output("[Debug] ", c, message);
}

void ConsoleLog::Information(const std::string& message)
{
    Output("[Info] ", message);
}

void ConsoleLog::Information(const LogContext& c, const std::string& message)
{
    Output("[Info] ", c, message);
}

void ConsoleLog::Error(const std::string& message)
{
    Output("[Error] ", message);
}

void ConsoleLog::Error(const LogContext& c, const std::string& message)
{
    Output("[Error] ", c, message);
}

}
