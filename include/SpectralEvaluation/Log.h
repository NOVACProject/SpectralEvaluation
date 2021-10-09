#pragma once

namespace novac
{

/** Simple interface for a basic logging class. */
class SimpleLog
{
public:

    virtual void Information(const char* message) = 0;
};

/** The simplest form of logging, using the console */
class ConsoleLog : public SimpleLog
{
public:

    virtual void Information(const char* message) override;
};



}