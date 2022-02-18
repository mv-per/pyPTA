#pragma once
#include <cmath>
#include <iostream>
using namespace std;

const double PI = 3.14159265359;
const double KB = 1.38064852e-23;

double DRA(double z, double eps0, double z0, double beta)
{
	return eps0 * pow(log(z0 / z), (1.0 / beta));
}

double LEE(double z, double eps_k, double sigma_ff, adsorbent adsorbent_properties)
{
	double epsilon_fs = eps_k * KB;
	double sigma_fs = (adsorbent_properties.sigma_ss + sigma_ff) / 2.0;
	double z_dummy = z + (adsorbent_properties.sigma_ss / 2.0);
	double summation = pow(sigma_fs, 4.0) / pow(z_dummy + 0.0 * adsorbent_properties.sigma_ss, 4.0) +
					   pow(sigma_fs, 4.0) / pow(z_dummy + 1.0 * adsorbent_properties.sigma_ss, 4.0) +
					   pow(sigma_fs, 4.0) / pow(z_dummy + 2.0 * adsorbent_properties.sigma_ss, 4.0) +
					   pow(sigma_fs, 4.0) / pow(z_dummy + 3.0 * adsorbent_properties.sigma_ss, 4.0);

	double PSI1 = 4.0 * PI * adsorbent_properties.rho_atoms * epsilon_fs * sigma_fs * sigma_fs;
	double PSI2 = pow(sigma_fs, 10.0) / (5.0 * pow(z_dummy, 10.0)) - 0.5 * summation;

	return PSI1 * PSI2;
}

double STEELE(double z, double eps_k, double sigma_ff, double sigma_ss)
{
	double epsilon_fs = eps_k * KB;
	double sigma_fs = (sigma_ss + sigma_ff) / 2.0; // Angstrom
	double DDelta = 3.35;						   // Angstroms
	double rho_s = 0.114;						   // Angstroms-3
	double PSI1 = 2.0 * PI * rho_s * epsilon_fs * sigma_fs * sigma_fs * DDelta;
	double division = sigma_fs / z;
	double PSI2 = 0.4 * pow(division, 10.0) - pow(division, 4.0) - pow(sigma_fs, 4.0) / (3.0 * DDelta * pow(z + 0.61 * DDelta, 3.0));
	//
	return PSI1 * PSI2;
}
