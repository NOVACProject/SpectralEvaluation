#include <SpectralEvaluation/Log.h>
#include <iostream>

using namespace novac;


void novac::ConsoleLog::Debug(const std::string& message)
{
    std::cout << "[Debug]" << message << std::endl;
}

void novac::ConsoleLog::Information(const std::string& message)
{
    std::cout << "[Info]" << message << std::endl;
}

void novac::ConsoleLog::Error(const std::string& message)
{
    std::cout << "[Error]" << message << std::endl;
}
