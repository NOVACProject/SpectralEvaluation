#include <SpectralEvaluation/Molecule.h>

namespace novac
{

Molecule::Molecule(StandardMolecule molec)
{
    switch (molec)
    {
    case StandardMolecule::SO2:
        name = "SO2";
        molecularWeight = 64.0638;
        break;
    case StandardMolecule::O3:
        name = "O3";
        molecularWeight = 47.9982;
        break;
    case StandardMolecule::BrO:
        name = "BrO";
        molecularWeight = 95.8980;
        break;
    case StandardMolecule::NO2:
        name = "NO2";
        molecularWeight = 46.0055; break;
    case StandardMolecule::HCHO:
        name = "HCHO";
        molecularWeight = 30.0206; break;
    default:
        name = "SO2";
        molecularWeight = 64.0638; break;
    }
}

double Molecule::Convert_MolecCm2_to_kgM2(double molec_cm2) const
{
    return molec_cm2 * (10.0) * this->molecularWeight / AVOGADROS_NUMBER;
}

}  // namespace novac
