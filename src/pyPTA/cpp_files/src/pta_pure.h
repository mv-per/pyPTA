#ifndef PTA_PURE_H
#define PTA_PURE_H

#include <cmath>
#include <vector>
#include <stdexcept>
#include <functional>

#include "data_classes.h"
#include "global_helper.h"
#include "pta_helper.h"
#include "pr77.h"
#include "adsorption_potentials.h"

/**
 * Perform Pure PTA calculations
 *
 * Shapiro, A. A., & Stenby, E. H. (1998).
 * Potential theory of multicomponent adsorption.
 * Journal of Colloid and Interface Science, 201(2), 146-157.
 *
 */
class PurePTA
{

public:
	std::string potential;
	std::size_t num_of_layers;
	std::string equation_of_state;
	std::string isotherm_type;
	/**
	 * PurePTA
	 *
	 * @param potential_ Name of the potential to be used in the calculations
	 * @param eos Name of the Equation of State to be used in the calculations
	 * @param isotherm_type_ Type of the isotherm to be used in calculations (excess | absolute).
	 * @param num_of_layers_ Number of layers from the solid to Gibbs phase.
	 */
	PurePTA(std::string potential_, std::string eos, std::string isotherm_type_, std::size_t num_of_layers_)
	{
		this->potential = potential_;
		this->num_of_layers = num_of_layers_;
		this->equation_of_state = eos;
		this->isotherm_type = isotherm_type_;
		// Set values based on inputs
	};

	~PurePTA(){};

	/**
	 * Get the calculated loading in a specific pressure
	 *
	 * @param P Pressure of the fluid
	 * @param T Temperature of the Fluid
	 * @param potential_params Params of the Adsorption Potential for this fluid
	 * @param fluid Fluid properties.
	 * @return Calculated adsorbed loading
	 */
	double get_loading(double P, double T, std::vector<double> potential_params, Fluid fluid, std::size_t num_of_layers)
	{

		call_mono_get_load _get_loading = get_proper_loader();
		call_potential get_potential = get_adsorption_potential(potential_params, fluid, this->adsorbent_properties);
		call_mono_eos eos_caller = get_equation_of_state(fluid);
		return _get_loading(P, T, potential_params, eos_caller, get_potential, num_of_layers);
	}

	/**
	 * Get the calculated loading for different pressures
	 *
	 * @param P List of pressure points to get loading
	 * @param T Temperature of the Fluid
	 * @param potential_params Params of the Adsorption Potential for this fluid
	 * @param fluid Fluid properties.
	 * @return List of calculated loadings
	 */
	std::vector<double> get_multiple_loadings(std::vector<double> P, double T, std::vector<double> potential_params, Fluid fluid, std::size_t num_of_layers)
	{
		std::vector<double>
			_loading(P.size(), 0.0);

		for (std::size_t i = 0; i < P.size(); i++)
		{
			_loading[i] = this->get_loading(P[i], T, potential_params, fluid, num_of_layers);
		}

		return _loading;
	}

	/**
	 * Get the deviation between calculated experimental loadings
	 *
	 * @param P List of pressure points to get loading
	 * @param T Temperature of the Fluid
	 * @param potential_params Params of the Adsorption Potential for this fluid
	 * @param fluid Fluid properties.
	 * @return List of calculated loadings
	 */
	double get_absolute_deviation(std::string deviation_type,
								  std::vector<double> loading_exp, std::vector<double> P, double T, std::vector<double> potential_params, Fluid fluid, std::size_t num_of_layers)
	{
		Deviation deviation = get_deviation_function(deviation_type);

		double difference = 0.;

		std::vector<double> calc_loading = this->get_multiple_loadings(P, T, potential_params, fluid, num_of_layers);

		for (std::size_t i = 0; i < P.size(); i++)
		{
			difference += deviation(loading_exp[i], calc_loading[i]);
		}

		if (deviation_type.find("relative") != std::string::npos)
		{
			return 100. / loading_exp.size() * difference;
		}

		return difference;
	}

	void set_adsorbent_property(Adsorbent properties)
	{
		this->adsorbent_properties = properties;
	}

private:
	Adsorbent adsorbent_properties;

	std::function<double(double, double, std::vector<double>, std::function<mono_eos(double, double)>, std::function<double(double)>, std::size_t)> get_proper_loader()
	{
		if (this->potential == "DRA")
		{
			return [this](double P, double T, std::vector<double> params, call_mono_eos eos, call_potential get_potential, std::size_t num_of_layers)
			{ return this->get_dr_loading(P, T, params, eos, get_potential, num_of_layers); };
		}
		else if (this->potential == "STEELE" || this->potential == "LEE")
		{
			return [this](double P, double T, std::vector<double> params, call_mono_eos eos, call_potential get_potential, std::size_t num_of_layers)
			{ return this->get_lj_loading(P, T, params, eos, get_potential, num_of_layers); };
		}
		else
		{

			throw std::invalid_argument("Isotherm potential not found/defined.");
		}
	}

	std::function<mono_eos(double, double)> get_equation_of_state(Fluid fluid_)
	{

		if (this->equation_of_state == "pr77")
		{
			return [=](double P, double T)
			{
				return pr77().get_mono_fluid_properties(P, T, fluid_);
			};
		}
		else
		{
			throw std::invalid_argument("Equation of State not found/defined.");
		}
	}

	std::function<double(double)> get_adsorption_potential(std::vector<double> potential_params, Fluid fluid, Adsorbent solid_properties)
	{
		if (this->potential == "DRA")
		{
			return [=](double z)
			{
				return DRA(z, potential_params[0], potential_params[1], potential_params[2]);
			};
		}
		else if (this->potential == "STELLE")
		{
			return [=](double z)
			{
				return STEELE(z, potential_params[0], fluid.lj_diameter, solid_properties);
			};
		}
		else if (this->potential == "LEE")
		{
			return [=](double z)
			{
				return LEE(z, potential_params[0], fluid.lj_diameter, solid_properties);
			};
		}
		else
		{
			throw std::invalid_argument("Isotherm potential not found/defined.");
		}
	}

	double get_f_eps(double eps, double T)
	{
		if (this->potential == "DRA")
		{
			return exp(eps / (8.314 * T));
		}
		else if (this->potential == "STELLE" || this->potential == "LEE")
		{
			return exp(-eps / (1.38064852e-23 * T));
		}
		else
		{
			throw std::invalid_argument("Isotherm potential not found/defined.");
		}
	}

	double get_density_value(double adsorbed_density, double bulk_density)
	{
		if (this->isotherm_type == "excess")
		{
			return adsorbed_density - bulk_density;
		}
		else if (this->isotherm_type == "absolute")
		{
			return adsorbed_density;
		}
		else
		{
			throw std::invalid_argument("Isotherm type not found/defined. Available types: excess, absolute");
		}
	}

	double get_dr_loading(double bulk_pressure, double temperature, std::vector<double> potential_params, call_mono_eos eos, call_potential get_potential, std::size_t num_of_layers)
	{

		double Pz, f_eps, delta, eps;

		struct mono_eos bulk = eos(bulk_pressure, temperature);

		// Initial estimates for the adsorbed region
		Pz = bulk_pressure;

		delta = potential_params[1] / (double)num_of_layers;

		std::vector<double>
			z = linespace(potential_params[1], delta, num_of_layers);

		std::vector<double> integral(num_of_layers, 0.0);

		for (std::size_t i = 0; i < z.size(); i++)
		{
			eps = get_potential(z[i]);
			f_eps = this->get_f_eps(eps, temperature);

			Pz = find_pz(Pz, bulk.fug, f_eps, temperature, eos);

			struct mono_eos adsorbed = eos(Pz, temperature);
			integral[i] = this->get_density_value(adsorbed.dens, bulk.dens) * 1e-3;
		}
		return integrate(integral, delta);
	}

	double get_lj_loading(double Pb, double T, std::vector<double> potential_params, call_mono_eos eos, call_potential get_potential, std::size_t num_of_layers)
	{
		double Pz, f_eps, delta, eps;
		struct mono_eos bulk = eos(Pb, T);

		// Initial estimates for the adsorbed region
		Pz = Pb;
		double minimal_space = 0.7 * potential_params[1];

		std::vector<double>
			z = linespace(potential_params[1] / 2.0, minimal_space, num_of_layers);

		delta = z[1] - z[0];

		std::vector<double> integral(num_of_layers, 0.0);

		for (std::size_t i = 0; i < z.size(); i++)
		{
			eps = get_potential(z[i]) + get_potential(potential_params[1] - z[i]);
			f_eps = this->get_f_eps(eps, T);

			Pz = find_pz(Pz, bulk.fug, f_eps, T, eos);

			struct mono_eos adsorbed = eos(Pz, T);
			integral[i] = this->get_density_value(adsorbed.dens, bulk.dens);
		}
		return -potential_params[2] * integrate(integral, delta) * 1e-7;
	}
};

#endif
