#include "pta_mixture.h"

/**
 * MixturePTA
 *
 * @param potential_ Name of the potential to be used in the calculations
 * @param eos Name of the Equation of State to be used in the calculations
 * @param isotherm_type_ Type of the isotherm to be used in calculations (excess | absolute).
 * @param num_of_layers_ Number of layers from the solid to Gibbs phase.
 */
MixturePTA::MixturePTA(std::string potential_, std::string eos, std::string isotherm_type_, std::size_t num_of_layers_)
{
	this->potential = potential_;
	this->num_of_layers = num_of_layers_;
	this->equation_of_state = eos;
	this->isotherm_type = isotherm_type_;
};

void MixturePTA::SetAdsorbent(Adsorbent adsorbent)
{
	this->adsorbent = adsorbent;
	this->AdsorbentConfigured = true;
}

std::vector<double> MixturePTA::GetLoading(std::vector<double> composition, double P, double T, std::vector<std::vector<double>> potential_params, std::vector<Fluid> fluids)
{
	call_mix_get_load _get_loading = this->GetLoadingFunction();
	return _get_loading(composition, P, T, potential_params, fluids);
}

std::function<mix_eos(std::vector<double>, double, double)> MixturePTA::get_equation_of_state_mixture(std::vector<Fluid> fluids)
{
	if (this->equation_of_state == "pr77")
	{
		return [=](std::vector<double> composition, double P, double T)
		{
			return pr77().get_mixture_fluid_properties(composition, P, T, fluids);
		};
	}
	else if (this->equation_of_state == "pr77-peneloux")
	{
		for(std::size_t i=0;i<fluids.size();i++){
			fluids[i]= CheckForFluidCriticalCompressibility(fluids[i]);
		}
		return [=](std::vector<double> composition, double P, double T)
		{
			return pr77_peneloux().get_mixture_fluid_properties(composition, P, T, fluids);
		};
	}
	else if (this->equation_of_state == "srk")
	{
		return [=](std::vector<double> composition, double P, double T)
		{
			return srk().get_mixture_fluid_properties(composition, P, T, fluids);
		};
	}
	else if (this->equation_of_state == "srk-peneloux")
	{
		for(std::size_t i=0;i<fluids.size();i++){
			fluids[i]= CheckForFluidCriticalCompressibility(fluids[i]);
		}
		return [=](std::vector<double> composition, double P, double T)
		{
			return srk_peneloux().get_mixture_fluid_properties(composition, P, T, fluids);
		};
	}
	else
	{
		throw std::invalid_argument("Equation of State not found/defined.");
	}
}

std::function<std::vector<double>(std::vector<double>, double, double, std::vector<std::vector<double>>, std::vector<Fluid>)> MixturePTA::GetLoadingFunction(void)
{
	if (this->potential == DRA_POTENTIAL)
	{
		return [this](std::vector<double> composition, double P, double T, std::vector<std::vector<double>> params, std::vector<Fluid> fluid_params)
		{ return this->get_mixture_dr_loading(composition, P, T, params, fluid_params); };
	}
	else if (this->potential == STEELE_POTENTIAL || this->potential == LEE_POTENTIAL)
	{
		if (!this->AdsorbentConfigured)
		{
			throw std::invalid_argument("Adsorbent properties are needed for LJ-based potentials and is not defined.\nUse `set_adsorbent(adsorbent)` to configure.");
		}
		return [this](std::vector<double> composition, double P, double T, std::vector<std::vector<double>> params, std::vector<Fluid> fluid_params)
		{ return this->get_mixture_lj_loading(composition, P, T, params, fluid_params); };
	}
	else
	{

		throw std::invalid_argument("Isotherm potential not found/defined.");
	}
}

std::function<double(std::size_t, double, Fluid)>
MixturePTA::GetAdsorptionPotentialInvoker(std::vector<std::vector<double>> potential_params)
{
	if (this->potential == DRA_POTENTIAL)
	{
		return [=](std::size_t i, double z, Fluid fluid)
		{
			return DRA(z, potential_params[i][0], potential_params[i][1], potential_params[i][2]);
		};
	}
	else if (this->potential == STEELE_POTENTIAL)
	{
		return [=](std::size_t i, double z, Fluid fluid)
		{
			return STEELE(z, potential_params[i][0], fluid.LennardJonnesDiameter, this->adsorbent);
		};
	}
	else if (this->potential == LEE_POTENTIAL)
	{
		return [=](std::size_t i, double z, Fluid fluid)
		{
			return LEE(z, potential_params[i][0], fluid.LennardJonnesDiameter, this->adsorbent);
		};
	}
	else
	{
		throw std::invalid_argument("Isotherm potential not found/defined.");
	}
}

std::vector<double>
MixturePTA::get_mixture_lj_loading(std::vector<double> bulk_composition, double bulk_pressure, double Temperature, std::vector<std::vector<double>> potential_parameters, std::vector<Fluid> fluids)
{

	std::size_t i;
	std::size_t j;
	double minimal_space, eps;
	std::size_t ncomp = bulk_composition.size();

	call_potential_mix get_potential = GetAdsorptionPotentialInvoker(potential_parameters);

	call_mix_eos eos = get_equation_of_state_mixture(fluids);

	struct mix_eos bulk = eos(bulk_composition, bulk_pressure, Temperature);

	// Initial Conditions
	double Pz = bulk_pressure;
	std::vector<double> adsorbed_composition = bulk_composition;

	std::vector<std::vector<double>> z(ncomp, std::vector<double>(this->num_of_layers, 0.0));
	std::vector<double> deltas(ncomp, 0.0);
	std::vector<double> f_eps(ncomp, 0.0);
	std::vector<double> adsorbed_loading(ncomp, 0.0);
	std::vector<std::vector<double>> mix_integral(ncomp, std::vector<double>(this->num_of_layers, 0.0));

	for (j = 0; j < ncomp; j++)
	{
		if (!fluids[j].LennardJonnesDiameter)
		{
			throw std::invalid_argument("Missing the Lennard-Jonnes diameter for the fluid " + fluids[j].Name);
		}

		minimal_space = this->min_int * fluids[j].LennardJonnesDiameter;

		z[j] = linespace(potential_parameters[j][1] / 2.0, minimal_space, this->num_of_layers);
		deltas[j] = z[j][1] - z[j][0];
	}

	for (i = 0; i < this->num_of_layers; i++)
	{
		for (j = 0; j < ncomp; j++)
		{
			eps = get_potential(j, z[j][i], fluids[j]) + get_potential(j, potential_parameters[j][1] - z[j][i], fluids[j]);
			f_eps[j] = GetAdsorptionPotentialEnergy(eps, Temperature, this->potential);

			// TODO: study the proper rate with parameter L and the adsorbent sizes
		}

		// Find pressure and composition of equilibrium
		struct Equilibrium Props = this->find_equilibrium_properties(Pz, adsorbed_composition, bulk.fug, f_eps, Temperature, eos);
		Pz = Props.Pz;
		adsorbed_composition = Props.x;

		struct mix_eos adsorbed = eos(adsorbed_composition, Pz, Temperature);

		for (j = 0; j < ncomp; j++)
		{
			mix_integral[j][i] = GetCalculatedAdsorbedDensity(adsorbed.dens * adsorbed_composition[j], bulk.dens * bulk_composition[j], this->isotherm_type);
		}
	}

	for (j = 0; j < ncomp; j++)
	{
		adsorbed_loading[j] = -potential_parameters[j][2] * integrate(mix_integral[j], deltas[j]) * 1.0e-7;
	}

	return adsorbed_loading;
}

std::vector<double>
MixturePTA::get_mixture_dr_loading(std::vector<double> bulk_composition, double bulk_pressure, double Temperature, std::vector<std::vector<double>> potential_parameters, std::vector<Fluid> fluids)
{
	std::size_t i, j;
	std::size_t ncomp = bulk_composition.size();
	double eps;

	call_potential_mix get_potential = GetAdsorptionPotentialInvoker(potential_parameters);

	call_mix_eos eos = get_equation_of_state_mixture(fluids);

	struct mix_eos bulk = eos(bulk_composition, bulk_pressure, Temperature);

	// Initial Conditions
	double Pz = bulk_pressure;
	std::vector<double> adsorbed_composition = bulk_composition;

	std::vector<std::vector<double>> z(ncomp, std::vector<double>(this->num_of_layers, 0.0));
	std::vector<double> deltas(ncomp, 0.0);

	for (j = 0; j < ncomp; j++)
	{
		deltas[j] = potential_parameters[j][1] / (double)this->num_of_layers;
		z[j] = linespace(potential_parameters[j][1], deltas[j], this->num_of_layers);
	}
	// Instantiate blank vectors
	std::vector<double> f_eps(ncomp, 0.0);
	std::vector<double> adsorbed_loading(ncomp, 0.0);
	std::vector<std::vector<double>> mix_integral(ncomp, std::vector<double>(this->num_of_layers, 0.0));

	for (i = 0; i < this->num_of_layers; i++)
	{
		for (j = 0; j < ncomp; j++)
		{
			eps = get_potential(j, z[j][i], fluids[j]);
			f_eps[j] = GetAdsorptionPotentialEnergy(eps, Temperature, this->potential);
		}

		// Find pressure and composition of equilibrium
		struct Equilibrium Props = this->find_equilibrium_properties(Pz, adsorbed_composition, bulk.fug, f_eps, Temperature, eos);
		Pz = Props.Pz;
		adsorbed_composition = Props.x;

		struct mix_eos adsorbed = eos(adsorbed_composition, Pz, Temperature);

		for (j = 0; j < ncomp; j++)
		{
			mix_integral[j][i] = GetCalculatedAdsorbedDensity(adsorbed.dens * adsorbed_composition[j], bulk.dens * bulk_composition[j], this->isotherm_type);
		}
	}

	for (j = 0; j < ncomp; j++)
	{
		adsorbed_loading[j] = integrate(mix_integral[j], deltas[j]) / 1000.0;
	}

	return adsorbed_loading;
}

double MixturePTA::get_equilibrium_difference(double Pz,
											  std::vector<double> adsorbed_composition,
											  std::vector<double> bulk_fugacity,
											  std::vector<double> f_eps,
											  double T,
											  call_mix_eos eos)
{

	if (std::isinf(Pz) == 1 || std::isnan(Pz) == 1 || Pz <= 0)
	{
		return 1.0;
	}

	double difference = 0.0;
	struct mix_eos adsorbed = eos(adsorbed_composition, Pz, T);

	// Find the Equilibrium difference
	for (size_t i = 0; i < adsorbed_composition.size(); i++)
	{
		difference += bulk_fugacity[i] * f_eps[i] - adsorbed.fug[i];
	}
	return difference;
}

struct Equilibrium MixturePTA::find_equilibrium_properties(double Pz_old, std::vector<double> xold, std::vector<double> fb,
														   std::vector<double> f_eps, double T, call_mix_eos eos)
{
	double xsum, composition_balancer;
	double Pz = Pz_old;
	std::size_t i, k;
	std::vector<double> new_adsorbed_composition = xold;
	struct mix_eos adsorbed = eos(new_adsorbed_composition, Pz, T);
	new_adsorbed_composition = get_new_adsorbed_fractions(fb, f_eps, adsorbed.phi, Pz);
	k = 0;

	do
	{

		// Set the Equilibrium lambda function
		auto min = [&](double x)
		{
			return get_equilibrium_difference(x, new_adsorbed_composition, fb, f_eps, T, eos);
		};

		// Find the Equilibrium pressure for the current composition
		Pz = brent_zeroin(min, Pz, 1.0e-8);

		// Get the fluid properties
		struct mix_eos adsorbed = eos(new_adsorbed_composition, Pz, T);

		// Finds the difference from the current and previous composition
		xsum = 0;
		for (i = 0; i < fb.size(); i++)
		{
			xsum += std::fabs(new_adsorbed_composition[i] / xold[i]) - 1.0;
		}

		xsum *= xsum;

		// Limits to 25 iterations
		if (k >= 25)
		{
			break;
		}

		// Update composition
		xold = new_adsorbed_composition;
		k++;
		new_adsorbed_composition = get_new_adsorbed_fractions(fb, f_eps, adsorbed.phi, Pz);

	} while (xsum >= 1.0e-12);
	// Keep updating the composition until it doesn't changes anymore

	struct Equilibrium data = {Pz, new_adsorbed_composition};
	return data;
}
