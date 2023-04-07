#include "adsorption_potentials.h"

/// Boltzmann constant, [m2 kg /s2 /K]
double BOLTZMANN_CONSTANT = 1.38064852e-23;

/// Mathematical constant, [-]
double PI = 3.14159265359;

//  spacing between basal planes, Angstroms
double BASAL_SPACING = 3.35;

double DRA(double z, double eps0, double z0, double beta)
{
    return eps0 * pow(log(z0 / z), (1.0 / beta));
}

double LEE(double z, double eps_k, double sigma_ff, Adsorbent adsorbent)
{
    if (!adsorbent.SolidDiameter || !adsorbent.SolidAtomicDensity)
    {
        throw std::invalid_argument("Missing adsorbent diameter and atomic density");
    }
    double epsilon_fs = eps_k * BOLTZMANN_CONSTANT;
    double sigma_fs = (adsorbent.SolidDiameter + sigma_ff) / 2.0;
    double z_dummy = z + (adsorbent.SolidDiameter / 2.0);
    double summation = pow(sigma_fs, 4.0) / pow(z_dummy + 0.0 * adsorbent.SolidDiameter, 4.0) +
                       pow(sigma_fs, 4.0) / pow(z_dummy + 1.0 * adsorbent.SolidDiameter, 4.0) +
                       pow(sigma_fs, 4.0) / pow(z_dummy + 2.0 * adsorbent.SolidDiameter, 4.0) +
                       pow(sigma_fs, 4.0) / pow(z_dummy + 3.0 * adsorbent.SolidDiameter, 4.0);

    double PSI1 = 4.0 * PI * adsorbent.SolidAtomicDensity * epsilon_fs * sigma_fs * sigma_fs;
    double PSI2 = pow(sigma_fs, 10.0) / (5.0 * pow(z_dummy, 10.0)) - 0.5 * summation;

    return PSI1 * PSI2;
}

/**
 * @brief Calculate the adsorption potential energy using the STEELE equation.
 *
 * This function calculates the adsorption potential energy using the STEELE equation,
 * which is a model used in the field of adsorption science to estimate the potential
 * energy of adsorption for a gas molecule on a solid adsorbent surface.
 *
 * @param z Charge of the adsorbed gas molecule.
 * @param eps_k Dielectric constant of the adsorbent material.
 * @param sigma_ff Lennard-jonnes diameter of the adsorbed gas molecule.
 * @param adsorbent Solid adsorbent material.
 *
 * @return Adsorption potential energy.
 */
double STEELE(double z, double eps_k, double sigma_ff, Adsorbent adsorbent)
{
    double epsilon_fs = eps_k * BOLTZMANN_CONSTANT;
    double sigma_fs = (adsorbent.SolidDiameter + sigma_ff) / 2.0; // Angstrom
    double division = sigma_fs / z;
    double PSI1 = 2.0 * PI * adsorbent.SolidAtomicDensity * epsilon_fs * sigma_fs * sigma_fs * BASAL_SPACING;
    double PSI2 = 0.4 * pow(division, 10.0) - pow(division, 4.0) - (pow(sigma_fs, 4.0) / (3.0 * BASAL_SPACING * pow(z + 0.61 * BASAL_SPACING, 3.0)));

    return PSI1 * PSI2;
}
