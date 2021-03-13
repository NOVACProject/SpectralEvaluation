#pragma once

#include <string>

namespace novac
{
    // This struct is used to represent a ratio between two species
    struct Ratio
    {
        double ratio; //<- The estimated quotient.
        double error; //<- An estimation of the error in the quotient.

        double minorResult; //<- The result of the minor evaluation (Bro).
        double minorError; //<- The estimated column error of minorResult.
        std::string minorSpecieName; //<- The name of the minor specie (Bro?)

        double majorResult; //<- The result of the major evaluation (SO2).
        double majorError; //<- The estimated column error of majorResult.
        std::string majorSpecieName; //<- The name of the minor specie (Bro?)
    };
}
