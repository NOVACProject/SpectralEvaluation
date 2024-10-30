#pragma once

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

}