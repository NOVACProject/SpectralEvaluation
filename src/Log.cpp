#include <SpectralEvaluation/Log.h>
#include <iostream>

using namespace novac;


void ConsoleLog::Information(const char* message)
{
    std::cout << message << std::endl;
}
