#pragma once

namespace novac
{

// Refractive index of air, as calculated by the formula from Edlén [Edlén, 1966]
inline double RefractiveIndexOfAir_Edlen(double nmVac)
{
    const double sigma_vac_2 = 1e-2 / (nmVac * nmVac);
    return 1.0 + 1e-8 * (8342.13 + 2406030.0 / (130.0 - sigma_vac_2) + 15997.0 / (38.9 - sigma_vac_2));
}

// Refractive index of air, as calculated by the formula from 
inline double RefractiveIndexOfAir_Morton(double nmVac)
{
    const double sigma_vac_2 = 1.0 / (nmVac * nmVac);
    return 1.0 + 0.0000834254 + 0.02406147 / (130.0 - sigma_vac_2) + 0.00015998 / (38.9 - sigma_vac_2);
}

// Converts a wavelength in vacuum to a wavelength in air using Edléns formula.
inline double NmVacuumToNmAir(double nmVac)
{
    return nmVac / RefractiveIndexOfAir_Edlen(nmVac);
}

} // namespace novac
