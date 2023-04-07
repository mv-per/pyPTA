#ifndef PTA_MIXTURE_H
#define PTA_MIXTURE_H

#include "pta_helper.h"
#include "helpers.h"
#include "equations_of_state/eos.h"
#include "optimization_algorithms/brent.h"
#include "adsorption_potentials.h"
#include "pta_solver.h"

#include <stdio.h>
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <cfloat>

class MixturePTA
{

public:
	std::string equation_of_state;
	std::string isotherm_type;
	std::string potential;
	std::size_t num_of_layers;
	~MixturePTA() {}

	/**
	 * MixturePTA
	 *
	 * @param potential_ potential to be used in the calculations
	 * @param eos Name of the Equation of State to be used in the calculations
	 * @param isotherm_type_ Type of the isotherm to be used in calculations (excess | absolute).
	 * @param num_of_layers_ Number of layers from the solid to Gibbs phase.
	 */
	MixturePTA(std::string potential_, std::string eos, std::string isotherm_type_, std::size_t num_of_layers_);

	std::vector<double> GetLoading(std::vector<double> composition, double P, double T, std::vector<std::vector<double>> potential_params, std::vector<Fluid> fluids);

	void SetAdsorbent(Adsorbent adsorbent);

private:
	bool AdsorbentConfigured = false;
	Adsorbent adsorbent;
	double min_int = 0.7;

	std::function<mix_eos(std::vector<double>, double, double)>
	get_equation_of_state_mixture(std::vector<Fluid> fluids);

	std::function<std::vector<double>(std::vector<double>, double, double, std::vector<std::vector<double>>, std::vector<Fluid>)> GetLoadingFunction(void);

	std::function<double(std::size_t, double, Fluid)>
	GetAdsorptionPotentialInvoker(std::vector<std::vector<double>> potential_params);

	std::vector<double>
	get_mixture_lj_loading(std::vector<double> bulk_composition, double bulk_pressure, double Temperature, std::vector<std::vector<double>> potential_parameters, std::vector<Fluid> fluids);

	std::vector<double>
	get_mixture_dr_loading(std::vector<double> bulk_composition, double bulk_pressure, double Temperature, std::vector<std::vector<double>> potential_parameters, std::vector<Fluid> fluids);

	double get_equilibrium_difference(double Pz,
									  std::vector<double> adsorbed_composition,
									  std::vector<double> bulk_fugacity,
									  std::vector<double> f_eps,
									  double T,
									  call_mix_eos eos);

	struct Equilibrium find_equilibrium_properties(double Pz_old, std::vector<double> xold, std::vector<double> fb,
												   std::vector<double> f_eps, double T, call_mix_eos eos);
};

#endif