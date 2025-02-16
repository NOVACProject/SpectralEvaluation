#pragma once

#include <SpectralEvaluation/Definitions.h>

namespace novac
{

/** Simple structure for remembering the scaling unit of a reference or absorption cross section */
enum class CrossSectionUnit
{
    unknown = 0,

    // An absorption cross section, scaled to the absorption of 2.5e15 molecules / cm2
    //  A reference in this unit means that a resulting DOAS fit will produce a column in ppmm.
    ppmm = 1,

    // cm2 per molecule, a typical unit for absorption cross sections. 
    //  A reference in this unit means that a resulting DOAS fit will produce a column in molecules / cm2.
    cm2_molecule = 2,

    // no unit at all.
    none = 99
};

// Representation of an angle, in degrees.
// Separate struct for this allows us to not mix up radians and degrees.
struct angle_degrees_t
{
public:
    angle_degrees_t() = default;
    angle_degrees_t(const angle_degrees_t& other) = default;
    angle_degrees_t(angle_degrees_t&& other) = default;
    angle_degrees_t(double v) : value(v) { }

    angle_degrees_t& operator=(const angle_degrees_t& other) { this->value = other.value; return *this; }
    angle_degrees_t& operator=(angle_degrees_t&& other) noexcept { this->value = other.value; return *this; }

    double value = 0.0;

    // Returns the angle converted to radians.
    double InRadians() const { return this->value * DEGREETORAD; }

    // Retuns the angle in degrees
    double InDegrees() const { return this->value; }
};

}