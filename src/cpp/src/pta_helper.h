#ifndef PTA_HELPER_H
#define PTA_HELPER_H

#include "data_classes.h"
#include "equations_of_state/eos.h"
#include <vector>
#include <functional>

typedef std::function<double(double)> call_potential;

typedef std::function<double(std::size_t, double, Fluid)> call_potential_mix;

typedef std::function<double(double, double, std::vector<double>, call_mono_eos, call_potential)> call_mono_get_load;

typedef std::function<std::vector<double>(std::vector<double>, double, double, std::vector<std::vector<double>>, std::vector<Fluid>)> call_mix_get_load;

struct Equilibrium
{
	double Pz;
	std::vector<double> x;
};

/**
 * Returns the first derivative of the equilibrium function
 *
 *
 * @param bulk_fugacities Vector of fugacities for each fluid in mixture.
 * @param f_eps Vector of Adsorption potential strengths for each fluid in mixture.
 * @param adsorbed_compressibility Vector of Compressibility factor for each fluid in the adsorbed mixture.
 * @param adsorbed_pressure Pressure of the adsorbed phase at a specific position.
 * @return Sum of the adsorbed component fractions minus 1.
 */
double mix_dfx(std::vector<double> bulk_fugacities,
			   std::vector<double> f_eps,
			   std::vector<double> adsorbed_compressibility,
			   double adsorbed_pressure);

/**
 * Returns the second derivative of the equilibrium function
 *
 *
 * @param bulk_fugacities Vector of fugacities for each fluid in mixture.
 * @param f_eps Vector of Adsorption potential strengths for each fluid in mixture.
 * @param adsorbed_compressibility Vector of Compressibility factor for each fluid in the adsorbed mixture.
 * @param adsorbed_pressure Pressure of the adsorbed phase at a specific position.
 * @return Sum of the adsorbed component fractions minus 1.
 */
double mix_ddfx(std::vector<double> fb, std::vector<double> f_eps, std::vector<double> adsorbed_compressibility, double Pz);

/**
 * Verifies the fractions of the adsorbed phase at a specific position.
 *
 *
 * @param bulk_fugacities Vector of fugacities for each fluid in mixture.
 * @param f_eps Vector of Adsorption potential strengths for each fluid in mixture.
 * @param adsorbed_compressibility Vector of Compressibility factor for each fluid in the adsorbed mixture.
 * @param adsorbed_pressure Pressure of the adsorbed phase at a specific position.
 * @return Sum of the adsorbed component fractions minus 1.
 */
double mix_fx(std::vector<double> bulk_fugacities,
			  std::vector<double> f_eps,
			  std::vector<double> adsorbed_compressibility,
			  double Pz);

double GetAdsorptionPotentialEnergy(double eps, double T, std::string potential);

double GetCalculatedAdsorbedDensity(double adsorbed_density, double bulk_density, std::string isotherm_type);

double integrate(std::vector<double> integral_vector, double step);

/**
 * Returns a new set of adsorbed fractions for a specific Pz
 *
 *
 * @param bulk_fugacities Vector of fugacities for each fluid in mixture.
 * @param f_eps Vector of Adsorption potential strengths for each fluid in mixture.
 * @param adsorbed_compressibility Vector of Compressibility factor for each fluid in the adsorbed mixture.
 * @param adsorbed_pressure Pressure of the adsorbed phase at a specific position.
 * @return Vector containing a new set of adsorbed fractions.
 */
std::vector<double> get_new_adsorbed_fractions(std::vector<double> bulk_fugacities,
											   std::vector<double> f_eps,
											   std::vector<double> adsorbed_compressibility,
											   double Pz);
#endif