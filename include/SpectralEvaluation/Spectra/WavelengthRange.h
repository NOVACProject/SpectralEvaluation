#pragma once

namespace novac
{
    /** Basic representation of a wavelength range. */
    struct WavelengthRange
    {
        WavelengthRange(double from, double to)
        {
            this->low = (from < to) ? from : to;
            this->high = (from < to) ? to : from;
        }

        // The lower edge of the range. Should always be smaller than 'high'.
        double low = 0.0;

        // The higher edge of the range. Should always be larger than 'low'.
        double high = 0.0;

        inline bool Empty() const { return high - low < 0.001; }

        // @return true if this WavelengthRange equals the other with the provided margin.
        bool Equals(const WavelengthRange& other, double margin) const
        {
            return std::abs(other.low - this->low) < margin && std::abs(other.high - this->high) < margin;
        }
    };

}
