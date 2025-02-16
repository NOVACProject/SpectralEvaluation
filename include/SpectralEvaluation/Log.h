#pragma once

#include <vector>
#include <string>
#include <iosfwd>

namespace novac
{

class CDateTime;

class LogContext
{
public:
    LogContext() = default;

    LogContext(std::string name, std::string value);

    ~LogContext() = default;

    friend std::ostream& operator<<(std::ostream& os, const LogContext& c);

    std::vector<std::pair<std::string, std::string>> properties;

    LogContext With(std::string name, std::string value);
    LogContext With(std::string name, int value);
    LogContext With(std::string name, double value);
    LogContext WithTimestamp(const CDateTime& value);

    /// List of recommended and commonly used 'name' parameters

    // The name of a file
    static const std::string FileName;

    // The name of a directory
    static const std::string Directory;

    // The name of a fit window
    static const std::string FitWindow;

    // Device serial number
    static const std::string Device;

    // Device model
    static const std::string DeviceModel;

    // A timestamp
    static const std::string Time;
};

/** Abstract logger base class. */
class ILogger
{
public:
    virtual void Debug(const std::string& message) = 0;
    virtual void Debug(const LogContext& c, const std::string& message) = 0;

    virtual void Information(const std::string& message) = 0;
    virtual void Information(const LogContext& c, const std::string& message) = 0;

    virtual void Error(const std::string& message) = 0;
    virtual void Error(const LogContext& c, const std::string& message) = 0;
};

/** The simplest form of logging, using the console */
class ConsoleLog : public ILogger
{
public:

    virtual void Debug(const std::string& message) override;
    virtual void Debug(const LogContext& c, const std::string& message) override;

    virtual void Information(const std::string& message) override;
    virtual void Information(const LogContext& c, const std::string& message) override;

    virtual void Error(const std::string& message) override;
    virtual void Error(const LogContext& c, const std::string& message) override;
};

}