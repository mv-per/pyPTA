#ifndef PTA_HELPER_H
#define PTA_HELPER_H

#include "pr77.h"
#include "data_classes.h"

struct Equilibrium
{
	double Pz;
	std::vector<double> x;
};

typedef std::function<double(double, double, std::vector<double>, std::function<mono_eos(double, double)>, std::function<double(double)>, std::size_t num_of_layers)> call_mono_get_load;

typedef std::function<mono_eos(double, double)> call_mono_eos;

typedef std::function<mix_eos(std::vector<double>, double, double)> call_mix_eos;
typedef std::function<std::vector<double>(std::vector<double>, double, double, std::vector<std::vector<double>>, std::vector<Fluid>)> call_mix_get_load;

typedef std::function<double(double)> call_potential;
typedef std::function<double(std::size_t, double)> call_potential_mix;

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
			  double Pz)
{
	double sum = 0.0;
	for (size_t i = 0; i < f_eps.size(); i++)
	{
		sum += (bulk_fugacities[i] * f_eps[i]) / (adsorbed_compressibility[i] * Pz);
	}
	return sum - 1;
}

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
			   double adsorbed_pressure)
{
	double sum = 0.0;
	for (size_t i = 0; i < f_eps.size(); i++)
	{
		sum += (bulk_fugacities[i] * f_eps[i]) / (adsorbed_compressibility[i] * adsorbed_pressure * adsorbed_pressure);
	}

	return -sum;
}

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
double mix_ddfx(std::vector<double> fb, std::vector<double> f_eps, std::vector<double> adsorbed_compressibility, double Pz)
{
	double sum = 0.0;
	for (size_t i = 0; i < f_eps.size(); i++)
	{
		sum += (2.0 * fb[i] * f_eps[i]) / (adsorbed_compressibility[i] * Pz * Pz * Pz);
	}

	return sum;
}

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
											   double Pz)
{
	size_t n = f_eps.size();
	std::vector<double> xnew(n, 0.0);
	double sum = 0.0;
	for (std::size_t i = 0; i < n; i++)
	{
		xnew[i] = bulk_fugacities[i] * f_eps[i] / (adsorbed_compressibility[i] * Pz);
		sum += xnew[i];
	}

	// Check the factions consistency
	for (std::size_t i = 0; i < n; i++)
	{
		xnew[i] = xnew[i] / sum;
	}
	return xnew;
}

// std::function<mono_eos(double, double, Fluid)> get_mono_equation_of_state(std::string eos)
// {

// 	if (eos == "pr77")
// 	{
// 		return [](double P, double T, Fluid fluid)
// 		{
// 			return pr77().get_mono_fluid_properties(P, T, fluid);
// 		};
// 	}
// 	else
// 	{
// 		throw std::invalid_argument("Equation of State not found/defined.");
// 	}
// }

double equilibrium(double p_,
				   double bulk_fugacity,
				   double f_eps,
				   double T,
				   call_mono_eos eos)
{
	if (std::isinf(p_) == 1 || std::isnan(p_) == 1)
	{
		return 10000.;
	}
	struct mono_eos adsorbed = eos(p_, T);
	return bulk_fugacity * f_eps - adsorbed.fug;
}

double find_pz(double Pz, double bulk_fugacity, double f_eps, double temperature, call_mono_eos eos)
{
	static std::function<double(double)> get_equilibrium = [&](double x)
	{
		return equilibrium(x, bulk_fugacity, f_eps, temperature, eos);
	};

	// static auto get_equilibrium = [&](double x)
	// {
	// 	return equilibrium(x, bulk_fugacity, f_eps, temperature, eos);
	// };

	double (*equilibrium)(double) = [](double p)
	{ return get_equilibrium(p); };
	return brent_zeroin(equilibrium, Pz, 1.0e-12);
}

double integrate(std::vector<double> integral_vector, double step)
{
	// if (this->integration_method == "Simpson")
	// {
	// 	return SIMPS(integral_vector, step, num_of_layers);
	// }

	// At this time, only simpson rule is implemented
	return simpson_rule(integral_vector, step);
}

double get_density_value(double adsorbed_density, double bulk_density, std::string isotherm_type)
{
	if (isotherm_type == "excess")
	{
		return adsorbed_density - bulk_density;
	}
	else if (isotherm_type == "absolute")
	{
		return adsorbed_density;
	}
	else
	{
		throw std::invalid_argument("Isotherm type not found/defined. Available types: excess, absolute");
	}
}

double get_f_eps(double eps, double T, std::string potential)
{
	if (potential == "DRA")
	{
		return std::exp(eps / (8.314 * T));
	}
	else if (potential == "STELLE" || potential == "LEE")
	{
		return std::exp(-eps / (1.38064852e-23 * T));
	}
	else
	{
		throw std::invalid_argument("Isotherm potential not found/defined.");
	}
}

// std::function<double(double)> get_adsorption_potential(std::vector<double> potential_params, Fluid fluid, std::string potential, Adsorbent adsorbent_properties = Adsorbent())
// {
// 	if (potential == "DRA")
// 	{
// 		return [=](double z)
// 		{
// 			return DRA(z, potential_params[0], potential_params[1], potential_params[2]);
// 		};
// 	}
// 	else if (potential == "STELLE")
// 	{
// 		return [=](double z)
// 		{
// 			return STEELE(z, potential_params[0], fluid.lj_diameter, adsorbent_properties);
// 		};
// 	}
// 	else if (potential == "LEE")
// 	{
// 		return [=](double z)
// 		{
// 			return LEE(z, potential_params[0], fluid.lj_diameter, adsorbent_properties);
// 		};
// 	}
// 	else
// 	{
// 		throw std::invalid_argument("Isotherm potential not found/defined.");
// 	}
// }

// double equilibrium(double x,
// 				   double fb,
// 				   double f_eps,
// 				   double T,
// 				   Fluid fluid,
// 				   std::string eos)
// {
// 	if (std::isinf(x) == 1 || std::isnan(x) == 1)
// 	{
// 		return -1;
// 	}
// 	call_mono_eos caller = get_mono_equation_of_state(eos);
// 	struct mono_eos adsorbed = caller(x, T, fluid);
// 	return fb * f_eps - adsorbed.fug;
// }

// double get_dr_loading(double bulk_pressure, double temperature, std::vector<double> potential_params, Fluid fluid, *PurePTA MyInstance)

// {
// 	double Pz, delta, eps, f_eps;
// 	call_mono_eos eos_caller = get_mono_equation_of_state(MyInstance->equation_of_state);
// 	call_potential get_potential = get_adsorption_potential(potential_params, fluid, MyInstance->potential, MyInstance->adsorbent_properties);

// 	struct mono_eos bulk = eos_caller(bulk_pressure, temperature, fluid);

// 	// Initial estimates for the adsorbed region
// 	Pz = bulk_pressure;

// 	delta = potential_params[1] / (double)MyInstance->num_of_layers;

// 	std::vector<double>
// 		z = linespace(potential_params[1], delta, MyInstance->num_of_layers);
// 	std::vector<double>
// 		integral(MyInstance->num_of_layers, 0.0);

// 	for (std::size_t i = 0; i < z.size(); i++)
// 	{
// 		eps = get_potential(z[i]);
// 		f_eps = get_f_eps(eps, temperature);

// 		static auto MINFUNC = [&](double x)
// 		{
// 			return equilibrium(x, bulk.fug, f_eps, temperature, fluid, MyInstance->equation_of_state)
// 		};
// 		double (*minfun_pointer)(double) = [](double x)
// 		{ return MINFUNC(x); };

// 		Pz = brent_zeroin(minfun_pointer, Pz, DEFAULT_TOL);

// 		struct mono_eos adsorbed = eos_caller(Pz, temperature, fluid);
// 		integral[i] = get_density_value(adsorbed.dens, bulk.dens, MyInstance->isotherm_type) * 1e-3;
// 	}
// 	return integrate(integral, delta);
// }

// double get_lj_loading(double Pb, double T, std::vector<double> potential_params, Fluid fluid, PurePTA MyInstance)
// {
// 	double Pz, delta, eps, f_eps;
// 	call_mono_eos eos_caller = get_mono_equation_of_state(MyInstance->equation_of_state);
// 	call_potential get_potential = get_adsorption_potential(potential_params, fluid, MyInstance->potential, MyInstance->adsorbent_properties);

// 	struct mono_eos bulk = eos_caller(Pb, T, fluid);

// 	// Initial estimates for the adsorbed region
// 	Pz = Pb;
// 	double minimal_space = MyInstance->min_int * potential_params[1];

// 	std::vector<double>
// 		z = linespace(potential_params[1] / 2.0, minimal_space, MyInstance->num_of_layers);

// 	delta = z[1] - z[0];

// 	std::vector<double> integral(MyInstance->num_of_layers, 0.0);

// 	for (std::size_t i = 0; i < z.size(); i++)
// 	{
// 		eps = get_potential(z[i]) + get_potential(potential_params[1] - z[i]);
// 		f_eps = get_f_eps(eps, T);

// 		static auto MINFUNC = [&](double p)
// 		{
// 			return equilibrium(p, bulk.fug, f_eps, T, fluid, MyInstance->equation_of_state);
// 		};
// 		double (*minfun_pointer)(double) = [](double p)
// 		{ return MINFUNC(p); };
// 		Pz = brent_zeroin(minfun_pointer, Pz, 1.0e-16);

// 		struct mono_eos adsorbed = eos_caller(Pz, T, fluid);
// 		integral[i] = get_density_value(adsorbed.dens, bulk.dens);
// 	}
// 	return -potential_params[2] * integrate(integral, delta) * 1e-7;
// }

#endif