#include "pta_pure.h"
#include <iostream>

PurePTA::PurePTA(std::string potential, std::string equation_of_state, std::string isotherm_type, std::size_t num_of_layers)
{
    this->Potential = potential;
    this->NumberOfLayers = num_of_layers;
    this->EquationOfState = equation_of_state;
    this->IsothermType = isotherm_type;
    this->Loader = GetLoadingFunction();
};

/**
 * Get the calculated loading in a specific pressure
 *
 * @param P Pressure of the fluid
 * @param T Temperature of the Fluid
 * @param potential_params Params of the Adsorption Potential for this fluid
 * @param fluid Fluid properties.
 * @return Calculated adsorbed loading
 */
double PurePTA::GetLoading(double P, double T, std::vector<double> potential_params, Fluid fluid)
{
    this->fluid = fluid;
    call_potential get_potential = GetAdsorptionPotentialInvoker(potential_params, fluid, this->adsorbent);
    call_mono_eos eos_caller = GetEquationOfStateInvoker(fluid);
    double Loading = this->Loader(P, T, potential_params, eos_caller, get_potential);
    return Loading;
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
std::vector<double> PurePTA::GetLoadings(std::vector<double> P, double T, std::vector<double> potential_params, Fluid fluid)
{
    std::vector<double>
        loadings(P.size());

    this->fluid = fluid;
    call_potential get_potential = GetAdsorptionPotentialInvoker(potential_params, fluid, this->adsorbent);
    call_mono_eos eos_caller = GetEquationOfStateInvoker(fluid);

    for (std::size_t i = 0; i < P.size(); i++)
    {
        loadings[i] = this->Loader(P[i], T, potential_params, eos_caller, get_potential);
    }

    return loadings;
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
double PurePTA::GetDeviationRange(std::string deviation_type,
                                  std::vector<double> loading_exp, std::vector<double> P, double T, std::vector<double> potential_params, Fluid fluid)
{
    deviation_caller deviation = get_deviation_function(deviation_type);

    double difference = 0.;

    std::vector<double> calc_loading = this->GetLoadings(P, T, potential_params, fluid);

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

void PurePTA::SetAdsorbent(Adsorbent adsorbent)
{
    this->adsorbent = adsorbent;
    this->AdsorbentConfigured = true;
}

call_mono_get_load PurePTA::GetLoadingFunction()
{

    if (this->Potential == DRA_POTENTIAL)
    {
        return [this](double P, double T, std::vector<double> params, call_mono_eos eos, call_potential get_potential)
        { return this->GetDRAAdsorbedAmout(P, T, params, eos, get_potential); };
    }
    else if (this->Potential == STEELE_POTENTIAL || this->Potential == LEE_POTENTIAL)
    {
        return [this](double P, double T, std::vector<double> params, call_mono_eos eos, call_potential get_potential)
        { 
            if (!this->AdsorbentConfigured) {
                throw std::invalid_argument("Adsorbent properties are needed for LJ-based potentials and is not defined.\nUse `set_adsorbent(adsorbent)` to configure.");
            }
            return this->GetLJAdsorbedAmout(P, T, params, eos, get_potential); };
    }
    else
    {

        throw std::invalid_argument("Isotherm potential not found/defined.");
    }
}

std::function<mono_eos(double, double)> PurePTA::GetEquationOfStateInvoker(Fluid fluid_)
{

    if (this->EquationOfState == "pr77")
    {
        return [=](double P, double T)
        {
            return pr77().get_mono_fluid_properties(P, T, fluid_);
        };
    }
    else if (this->EquationOfState == "pr77-peneloux")
    {
        return [=](double P, double T)
        {
            return pr77_peneloux().get_mono_fluid_properties(P, T, fluid_);
        };
    }
    else if (this->EquationOfState == "srk")
    {
        return [=](double P, double T)
        {
            return srk().get_mono_fluid_properties(P, T, fluid_);
        };
    }
    else if (this->EquationOfState == "srk-peneloux")
    {
        return [=](double P, double T)
        {
            return srk_peneloux().get_mono_fluid_properties(P, T, fluid_);
        };
    }
    else
    {
        throw std::invalid_argument("Equation of State not found/defined.");
    }
}

std::function<double(double)> PurePTA::GetAdsorptionPotentialInvoker(std::vector<double> potential_params, Fluid fluid, Adsorbent solid_properties)
{
    if (this->Potential == DRA_POTENTIAL)
    {
        return [=](double z)
        {
            return DRA(z, potential_params[0], potential_params[1], potential_params[2]);
        };
    }
    else if (this->Potential == STEELE_POTENTIAL)
    {
        return [=](double z)
        {
            return STEELE(z, potential_params[0], fluid.LennardJonnesDiameter, solid_properties);
        };
    }
    else if (this->Potential == LEE_POTENTIAL)
    {
        return [=](double z)
        {
            return LEE(z, potential_params[0], fluid.LennardJonnesDiameter, solid_properties);
        };
    }
    else
    {
        throw std::invalid_argument("Isotherm potential not found/defined.");
    }
}

double PurePTA::GetDRAAdsorbedAmout(double bulk_pressure, double temperature, std::vector<double> potential_params, call_mono_eos eos, call_potential get_potential)
{

    double Pz, f_eps, delta, eps;

    struct mono_eos bulk = eos(bulk_pressure, temperature);
    PTASolver solver = PTASolver(bulk.fug, temperature, eos);

    // Initial estimates for the adsorbed region
    Pz = bulk_pressure;

    delta = potential_params[1] / (double)this->NumberOfLayers;

    std::vector<double>
        z = linespace(potential_params[1], delta, this->NumberOfLayers);

    std::vector<double> integral(this->NumberOfLayers, 0.0);

    for (std::size_t i = 0; i < z.size(); i++)
    {
        eps = get_potential(z[i]);
        f_eps = GetAdsorptionPotentialEnergy(eps, temperature, this->Potential);
        Pz = solver.findOptimizedPressure(Pz, f_eps);

        struct mono_eos adsorbed = eos(Pz, temperature);
        integral[i] = GetCalculatedAdsorbedDensity(adsorbed.dens, bulk.dens, this->IsothermType) * 1e-3;
    }
    return integrate(integral, delta);
}

double PurePTA::GetLJAdsorbedAmout(
    double BulkPressure,
    double Temperature,
    std::vector<double> PotentialParameters,
    call_mono_eos EquationOfStateInvoker,
    call_potential get_potential)
{
    double Pz, f_eps, delta, eps;

    struct mono_eos bulk = EquationOfStateInvoker(BulkPressure, Temperature);
    PTASolver solver = PTASolver(bulk.fug, Temperature, EquationOfStateInvoker);

    // Initial estimates for the adsorbed region
    Pz = BulkPressure;
    if (!this->fluid.LennardJonnesDiameter)
    {
        throw std::invalid_argument("Missing the Lennard-Jonnes diameter for the fluid " + this->fluid.Name);
    }
    double minimal_space = 0.7 * this->fluid.LennardJonnesDiameter;

    std::vector<double>
        z = linespace(PotentialParameters[1] / 2.0, minimal_space, this->NumberOfLayers);

    delta = z[1] - z[0];

    std::vector<double> integral(this->NumberOfLayers, 0.0);

    for (std::size_t i = 0; i < z.size(); i++)
    {
        eps = get_potential(z[i]) + get_potential(PotentialParameters[1] - z[i]);
        f_eps = GetAdsorptionPotentialEnergy(eps, Temperature, this->Potential);
        Pz = solver.findOptimizedPressure(Pz, f_eps);

        struct mono_eos adsorbed = EquationOfStateInvoker(Pz, Temperature);
        integral[i] = GetCalculatedAdsorbedDensity(adsorbed.dens, bulk.dens, this->IsothermType);
    }
    return -PotentialParameters[2] * integrate(integral, delta) * 1e-7;
}
