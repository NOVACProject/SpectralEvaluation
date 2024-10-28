#include <SpectralEvaluation/NullableValue.h>
#include <iostream>

namespace novac
{

std::ostream& operator<<(std::ostream& os, const Nullable<double>& value)
{
    char buffer[128];
    sprintf(buffer, "%lf", value.Value());
    os << std::string(buffer);
    return os;
}

std::ostream& operator<<(std::ostream& os, const Nullable<int>& value)
{
    char buffer[128];
    sprintf(buffer, "%d", value.Value());
    os << std::string(buffer);
    return os;
}


}