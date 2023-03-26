#include "adsorption_potentials.h"

/// Boltzmann constant, [m2 kg /s2 /K]
double BOLTZMANN_CONSTANT = 1.38064852e-23;

/// Mathematical constant, [-]
double PI = 3.14159265359;

double DRA(double z, double eps0, double z0, double beta)
{
    return eps0 * pow(log(z0 / z), (1.0 / beta));
}

double LEE(double z, double eps_k, double sigma_ff, Adsorbent adsorbent)
{
    double epsilon_fs = eps_k * BOLTZMANN_CONSTANT;
    double sigma_fs = (adsorbent.SolidDiameter + sigma_ff) / 2.0;
    double z_dummy = z + (adsorbent.SolidDiameter / 2.0);
    double summation = pow(sigma_fs, 4.0) / pow(z_dummy + 0.0 * adsorbent.SolidDiameter, 4.0) +
                       pow(sigma_fs, 4.0) / pow(z_dummy + 1.0 * adsorbent.SolidDiameter, 4.0) +
                       pow(sigma_fs, 4.0) / pow(z_dummy + 2.0 * adsorbent.SolidDiameter, 4.0) +
                       pow(sigma_fs, 4.0) / pow(z_dummy + 3.0 * adsorbent.SolidDiameter, 4.0);

    double PSI1 = 4.0 * PI * adsorbent.SolidAtomicDensity * epsilon_fs * sigma_fs * sigma_fs;
    double PSI2 = std::pow(sigma_fs, 10.0) / (5.0 * std::pow(z_dummy, 10.0)) - 0.5 * summation;

    return PSI1 * PSI2;
}

double STEELE(double z, double eps_k, double sigma_ff, Adsorbent adsorbent)
{
    double epsilon_fs = eps_k * BOLTZMANN_CONSTANT;
    double sigma_fs = (adsorbent.SolidDiameter + sigma_ff) / 2.0; // Angstrom
    double d_delta = 3.35;                                        // Angstroms
    double rho_s = 0.114;                                         // Angstroms-3
    double division = sigma_fs / z;
    double PSI1 = 2.0 * PI * rho_s * epsilon_fs * sigma_fs * sigma_fs * d_delta;
    double PSI2 = 0.4 * std::pow(division, 10.0) - std::pow(division, 4.0) - std::pow(sigma_fs, 4.0) / (3.0 * d_delta * std::pow(z + 0.61 * d_delta, 3.0));

    return PSI1 * PSI2;
}
