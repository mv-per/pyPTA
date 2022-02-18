#include <cmath>
#include <vector>
#include <stdexcept>
#include <functional>

#include <iostream>

#include "pta_helper.h"
#include "pr77.h"
#include "ADSORPTION_POTENTIAL_FUNCTIONS.h"
#include "GLOBAL_FUNCTIONS.h"

typedef std::function<mono_eos(double, double)> call_eos;
typedef std::function<double(double)> call_potential;
typedef std::function<double(double, double, std::vector<double>, std::function<mono_eos(double, double)>)> call_get_load;

double absolute_deviation(double nexp, double nads)
{
	return fabs(nexp - nads);
}

class PurePTA
{

public:
	PurePTA(std::string potential_, std::string eos, std::string isotherm_type_, std::size_t num_of_layers_)
	{
		this->potential = potential_;
		this->num_of_layers = num_of_layers_;
		this->equation_of_state = eos;
		this->isotherm_type = isotherm_type_;
		// Set values based on inputs
		this->integral.resize(num_of_layers);
		this->_get_loading = get_proper_loader();
	};

	double get_loading(double P, double T, std::vector<double> potential_params, Fluid fluid)
	{
		this->eos_caller = get_equation_of_state(fluid);
		this->get_potential = get_adsorption_potential(potential_params);

		return _get_loading(P, T, potential_params, eos_caller);
	}

	vector<double> get_multiple_loadings(std::vector<double> P, double T, std::vector<double> potential_params, Fluid fluid)
	{
		this->eos_caller = get_equation_of_state(fluid);
		this->get_potential = get_adsorption_potential(potential_params);

		std::vector<double>
			_loading(P.size(), 0.0);

		for (std::size_t i = 0; i < P.size(); i++)
		{
			_loading[i] = this->_get_loading(P[i], T, potential_params, this->eos_caller);
		}

		return _loading;
	}

	double get_absolute_deviation(std::vector<double> loading_exp, std::vector<double> P, double T, std::vector<double> potential_params, Fluid fluid)
	{
		this->eos_caller = get_equation_of_state(fluid);
		this->get_potential = get_adsorption_potential(potential_params);

		double difference = 0.;

		for (std::size_t i = 0; i < P.size(); i++)
		{
			difference += absolute_deviation(loading_exp[i], this->_get_loading(P[i], T, potential_params, this->eos_caller));
		}

		return difference;
	}

	void set_adsorbent_property(adsorbent properties)
	{
		this->adsorbent_properties = properties;
	}

private:
	std::vector<double> integral;
	std::string integration_method = "Simpson";
	std::size_t num_of_layers;
	std::string potential;
	std::string equation_of_state;
	std::string isotherm_type;
	call_get_load _get_loading;
	call_potential get_potential;
	call_eos eos_caller;
	adsorbent adsorbent_properties;
	double Pz, f_eps, delta, eps;
	double PSI_z, PSI_lz;

	std::function<double(double, double, std::vector<double>, std::function<mono_eos(double, double)>)> get_proper_loader()
	{
		if (this->potential == "DRA")
		{
			return [this](double P, double T, std::vector<double> params, call_eos eos)
			{ return this->get_dr_loading(P, T, params, eos); };
		}
		else if (this->potential == "STEELE" || this->potential == "LEE")
		{
			return [this](double P, double T, std::vector<double> params, call_eos eos)
			{ return this->get_lj_loading(P, T, params, eos); };
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

	std::function<double(double)> get_adsorption_potential(std::vector<double> potential_params)
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
				return STEELE(z, potential_params[0], potential_params[1], potential_params[2]);
			};
		}
		else if (this->potential == "LEE")
		{
			return [=](double z)
			{
				return LEE(z, potential_params[0], potential_params[1], this->adsorbent_properties);
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

	double equilibrium(double p_,
					   double bulk_fugacity,
					   double f_eps,
					   double T,
					   call_eos eos)
	{
		if (isinf(p_) == 1)
		{
			return 10000.;
		}
		struct mono_eos adsorbed = eos(p_, T);
		return bulk_fugacity * f_eps - adsorbed.fug;
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

	double integrate(std::vector<double> integral_vector, double step)
	{
		// if (this->integration_method == "Simpson")
		// {
		// 	return SIMPS(integral_vector, step, num_of_layers);
		// }

		// At this time, only simpson rule is implemented
		return SIMPS(integral_vector, step, this->num_of_layers);
	}

	double get_dr_loading(double bulk_pressure, double temperature, std::vector<double> potential_params, call_eos eos)
	{
		struct mono_eos bulk = eos(bulk_pressure, temperature);

		// Initial estimates for the adsorbed region
		this->Pz = bulk_pressure;

		this->delta = potential_params[1] / (double)this->num_of_layers;

		std::vector<double>
			z = linespace(potential_params[1], delta, this->num_of_layers);

		for (std::size_t i = 0; i < z.size(); i++)
		{
			this->eps = get_potential(z[i]);
			this->f_eps = get_f_eps(eps, temperature);

			static auto MINFUNC = [=](double p)
			{
				return equilibrium(p, bulk.fug, this->f_eps, temperature, eos);
			};

			double (*minfun_pointer)(double) = [](double p)
			{ return MINFUNC(p); };
			this->Pz = brent_zeroin(minfun_pointer, this->Pz, 1.0e-16);

			struct mono_eos adsorbed = eos(Pz, temperature);
			this->integral[i] = get_density_value(adsorbed.dens, bulk.dens) * 1e-3;
		}
		return integrate(this->integral, this->delta);
	}

	double get_lj_loading(double Pb, double T, std::vector<double> potential_params, call_eos eos)
	{
		struct mono_eos bulk = eos(Pb, T);

		// Initial estimates for the adsorbed region
		this->Pz = Pb;
		double minimal_space = 0.7 * potential_params[2];

		vector<double>
			z = linespace(potential_params[1] / 2.0, minimal_space, this->num_of_layers);

		this->delta = z[1] - z[0];

		for (std::size_t i = 0; i < z.size(); i++)
		{
			this->eps = get_potential(z[i]) + get_potential(potential_params[1] - z[i]);
			this->f_eps = get_f_eps(eps, T);

			static auto MINFUNC = [&](double p)
			{
				return equilibrium(p, bulk.fug, this->f_eps, T, eos);
			};
			double (*minfun_pointer)(double) = [](double p)
			{ return MINFUNC(p); };
			this->Pz = brent_zeroin(minfun_pointer, this->Pz, 1.0e-16);

			struct mono_eos adsorbed = eos(Pz, T);
			this->integral[i] = get_density_value(adsorbed.dens, bulk.dens);
		}
		return -potential_params[3] * integrate(this->integral, this->delta) * 1e-7;
	}
};
