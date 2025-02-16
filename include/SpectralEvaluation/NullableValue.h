#pragma once

#include <iosfwd>
#include <stdexcept>

namespace novac
{

/// Nullable represents a value which may or may not have been set.
/// Until the value has been set, it should not be used and doing so results in an exception.
template<class T> struct Nullable {

private:
    T m_value = T();
    bool m_isSet = false;

    void VerifyIsSet() const
    {
        if (!m_isSet)
        {
            throw std::logic_error("Attempted to use a nullable value before it has been set.");
        }
    }

public:
    Nullable() = default;
    // Creates a new Nullable with the provided value. This value is then set and Value() returns the set value.
    Nullable(T value) { m_value = value; m_isSet = true; }

    Nullable(const Nullable& other) = default;
    Nullable(Nullable&& other) = default;

    Nullable& operator=(const Nullable& other) = default;
    Nullable& operator=(Nullable&& other) = default;

    Nullable& operator+(double value) { static_assert(std::is_arithmetic<T>::value, "Not an arithmetic type"); VerifyIsSet(); return this->m_value + value; }
    Nullable& operator-(double value) { static_assert(std::is_arithmetic<T>::value, "Not an arithmetic type"); VerifyIsSet(); return this->m_value - value; }
    Nullable& operator*(double value) { static_assert(std::is_arithmetic<T>::value, "Not an arithmetic type"); VerifyIsSet(); return this->m_value * value; }
    Nullable& operator/(double value) { static_assert(std::is_arithmetic<T>::value, "Not an arithmetic type"); VerifyIsSet(); return this->m_value / value; }


    Nullable& operator+=(double value) { static_assert(std::is_arithmetic<T>::value, "Not an arithmetic type"); VerifyIsSet(); m_value += value; return *this; }
    Nullable& operator-=(double value) { static_assert(std::is_arithmetic<T>::value, "Not an arithmetic type"); VerifyIsSet(); m_value -= value; return *this; }
    Nullable& operator*=(double value) { static_assert(std::is_arithmetic<T>::value, "Not an arithmetic type"); VerifyIsSet(); m_value *= value; return *this; }

    // Value() retrieves the current value.
    // Throws std::logic_error if the value has not yet been set.
    T Value() const
    {
        VerifyIsSet();
        return m_value;
    }

    T ValueOrDefault(T defaultValue) const
    {
        if (!m_isSet)
        {
            return defaultValue;
        }
        return m_value;
    }

    bool HasValue() const { return m_isSet; }

    void Set(T value) { m_value = value; m_isSet = true; }
};

std::ostream& operator<<(std::ostream& os, const Nullable<double>& value);
std::ostream& operator<<(std::ostream& os, const Nullable<int>& value);


}