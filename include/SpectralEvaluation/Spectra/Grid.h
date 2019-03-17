#pragma once

struct UniformGrid
{
    double minValue = 0.0;
    double maxValue = 0.0;
    size_t length = 0U;

    double Resolution() const { return (maxValue - minValue) / (double)(length - 1); }

    inline double At(size_t idx) const
    {
        if (idx < length)
        {
            return minValue + (maxValue - minValue) * (double)idx / (double)(length - 1);
        }
        return 0.0; // invalid index
    }

    void Generate(std::vector<double>& values)
    {
        values.resize(this->length);
        for (size_t ii = 0; ii < length; ++ii)
        {
            values[ii] = minValue + (maxValue - minValue) * (double)ii / (length - 1);
        }
    }
};