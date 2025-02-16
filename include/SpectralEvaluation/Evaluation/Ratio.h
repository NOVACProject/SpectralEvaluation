#pragma once

#include <string>

namespace novac
{
// This struct is used to represent a ratio between two species (e.g. BrO to SO2)
struct Ratio
{
    // The estimated quotient between the two columns.
    double ratio = 0.0;

    // An estimation of the error in the quotient. Calculated from the retrieved columns and their Doas fit errors
    double error = 0.0;

    // The evaluated column of the minor specie (BrO).
    double minorResult = 0.0;

    // The error in the evaluated column of the minor specie (BrO)
    double minorError = 0.0;

    // The name of the minor specie (e.g. "BrO")
    std::string minorSpecieName = "";

    // The evaluated column of the major specie (SO2)
    double majorResult = 0.0;

    // The error in the evaluated column of the major specie (SO2)
    double majorError = 0.0;

    // The name of the major specie (e.g. "SO2")
    std::string majorSpecieName = "";
};
}
