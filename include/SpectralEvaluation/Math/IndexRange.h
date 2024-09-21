#pragma once

#include <stddef.h>

namespace novac
{

/** Basic representation of a range of indices. */
struct IndexRange
{
    IndexRange(size_t low, size_t high)
    {
        this->from = (low < high) ? low : high;
        this->to = (low < high) ? high : low;
    }

    size_t from = 0;

    size_t to = 0;

    inline bool Empty() const { return to <= from; }

    inline size_t Length() const { return to - from; }
};

}
