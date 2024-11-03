#pragma once

#include <string>

namespace novac
{

/** The class <b>Molecule</b> is used to store information on a single
    trace gas species, such as it's molecular weight.
    This is also able to use these numbers to convert between columns in number
    of molecules and in mass.
*/
enum class StandardMolecule
{
    SO2,
    O3,
    BrO,
    NO2,
    HCHO
};

std::string ToString(StandardMolecule molecule);

const double AVOGADROS_NUMBER = 6.02214179e23;

/** The class <b>Molecule</b> is used to pass the properties of an individual
    molecule to different parts of the program. The main usage (for now) is
    passing on the molecular weight to convert from molec/cm2 to kg/m2.
*/

struct Molecule
{
public:

    Molecule() = default;

    Molecule(StandardMolecule molec);

    // ----------------------------------------------------------------------
    // ---------------------- PUBLIC DATA -----------------------------------
    // ----------------------------------------------------------------------

    /** The trivial name of this gas molecule */
    std::string name = "SO2";

    /** The molecular weight of the molecule
        In grams per mol (g/mol) or (u/molecule).
        Defaults to the molecular weight of SO2. */
    double molecularWeight = 64.0638;

    // ----------------------------------------------------------------------
    // --------------------- PUBLIC METHODS ---------------------------------
    // ----------------------------------------------------------------------

    /** Takes a column in molec/cm2 and converts it to kg/m2 */
    double Convert_MolecCm2_to_kgM2(double molec_cm2) const;

};

} // namespace novac