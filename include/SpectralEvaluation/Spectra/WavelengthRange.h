#pragma once

namespace Evaluation
{

struct WavelengthRange
{
    WavelengthRange(double from, double to)
    {
        this->low = (from < to) ? from : to;
        this->high = (from < to) ? to : from;
    }
    double low = 0.0;
    double high = 0.0;

    inline bool Empty() const { return high - low < 0.001; }
};

}
