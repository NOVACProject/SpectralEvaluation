#include <SpectralEvaluation/Molecule.h>

namespace novac
{

std::string ToString(StandardMolecule molecule)
{
    switch (molecule)
    {
    case novac::StandardMolecule::SO2:
        return "SO2";
    case novac::StandardMolecule::O3:
        return "O3";
    case novac::StandardMolecule::BrO:
        return "BrO";
    case novac::StandardMolecule::NO2:
        return "NO2";
    case novac::StandardMolecule::HCHO:
        return "HCHO";
    default:
        return "SO2"; // default assumption
    }
}

Molecule::Molecule(StandardMolecule molec)
{
    switch (molec)
    {
    case StandardMolecule::SO2:
        molecularWeight = 64.0638;
        break;
    case StandardMolecule::O3:
        molecularWeight = 47.9982;
        break;
    case StandardMolecule::BrO:
        molecularWeight = 95.8980;
        break;
    case StandardMolecule::NO2:
        molecularWeight = 46.0055; break;
    case StandardMolecule::HCHO:
        molecularWeight = 30.0206; break;
    default:
        molecularWeight = 64.0638; break;
    }

    name = ToString(molec);
}

double Molecule::Convert_MolecCm2_to_kgM2(double molec_cm2) const
{
    return molec_cm2 * (10.0) * this->molecularWeight / AVOGADROS_NUMBER;
}

}  // namespace novac
