#pragma once

#include <stddef.h>

namespace novac
{

// IndexRange is a basic representation of a range of indices.
// Notice that there is no default constructor here, the values must be set.
struct IndexRange
{
    IndexRange(size_t low, size_t high)
        : from((low < high) ? low : high), to((low < high) ? high : low)
    {}

    size_t from = 0;

    size_t to = 0;

    inline bool Empty() const { return to <= from; }

    inline size_t Length() const { return to - from; }

    bool operator==(const IndexRange& other) const { return this->from == other.from && this->to == other.to; }
};

}
