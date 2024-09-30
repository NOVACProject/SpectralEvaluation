#pragma once

#include <string>

namespace novac
{

/** Abstract logger base class. */
class ILogger
{
public:
    virtual void Debug(const std::string& message) = 0;

    virtual void Information(const std::string& message) = 0;

    virtual void Error(const std::string& message) = 0;
};

/** The simplest form of logging, using the console */
class ConsoleLog : public ILogger
{
public:

    // Inherited via ILogger
    virtual void Debug(const std::string& message) override;

    virtual void Information(const std::string& message) override;

    virtual void Error(const std::string& message) override;

};



}