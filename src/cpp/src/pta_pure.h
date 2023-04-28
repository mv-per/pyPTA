#ifndef PTA_PURE_H
#define PTA_PURE_H

#include <cmath>
#include <vector>
#include <stdexcept>
#include <functional>

#include "data_classes.h"
#include "helpers.h"
#include "pta_helper.h"
#include "pta_solver.h"
#include "equations_of_state/eos.h"
#include "adsorption_potentials.h"
#include "optimization_algorithms/bisection.h"
#include "optimization_algorithms/brent.h"
#include "optimization_algorithms/fmin.h"

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
	call_mono_get_load Loader;
	Adsorbent adsorbent;
	Fluid fluid;

public:
	std::string Potential;
	std::size_t NumberOfLayers;
	std::string EquationOfState;
	std::string IsothermType;

	~PurePTA() {}
	/**
	 * PurePTA
	 *
	 * @param potential Name of the potential to be used in the calculations
	 * @param equation_of_state Name of the Equation of State to be used in the calculations
	 * @param isotherm_type Type of the isotherm to be used in calculations (excess | absolute).
	 * @param num_of_layers Number of layers from the solid to Gibbs phase.
	 */
	PurePTA(std::string potential, std::string equation_of_state, std::string isotherm_type, std::size_t num_of_layers);

	/**
	 * Get the calculated loading in a specific pressure
	 *
	 * @param P Pressure of the fluid
	 * @param T Temperature of the Fluid
	 * @param potential_params Params of the Adsorption Potential for this fluid
	 * @param fluid Fluid properties.
	 * @return Calculated adsorbed loading
	 */
	double GetLoading(double P, double T, std::vector<double> potential_params, Fluid fluid);
	/**
	 * Get the calculated pressure in a specific loading
	 *
	 * @param n loading of the fluid
	 * @param T Temperature of the Fluid
	 * @param potential_params Params of the Adsorption Potential for this fluid
	 * @param fluid Fluid properties.
	 * @return Calculated adsorbed loading
	 */
	double GetPressure(double n, double T, std::vector<double> potential_params, Fluid fluid, double P_estimate_);
	/**
	 * Get the calculated loading for different pressures
	 *
	 * @param P List of pressure points to get loading
	 * @param T Temperature of the Fluid
	 * @param potential_params Params of the Adsorption Potential for this fluid
	 * @param fluid Fluid properties.
	 * @return List of calculated loadings
	 */
	std::vector<double> GetLoadings(std::vector<double> P, double T, std::vector<double> potential_params, Fluid fluid);
	/**
	 * Get the deviation between calculated experimental loadings
	 *
	 * @param P List of pressure points to get loading
	 * @param T Temperature of the Fluid
	 * @param potential_params Params of the Adsorption Potential for this fluid
	 * @param fluid Fluid properties.
	 * @return List of calculated loadings
	 */
	double GetDeviationRange(std::string deviation_type,
							 std::vector<double> loading_exp, std::vector<double> P, double T, std::vector<double> potential_params, Fluid fluid);

	void SetAdsorbent(Adsorbent properties);

private:
	bool AdsorbentConfigured = false;

	call_mono_get_load GetLoadingFunction();

	std::function<mono_eos(double, double)> GetEquationOfStateInvoker(Fluid fluid_);

	std::function<double(double)> GetAdsorptionPotentialInvoker(std::vector<double> potential_params, Fluid fluid, Adsorbent solid_properties);
	/**
	 * Get the potential function based on its energy interation and temperature
	 *
	 */
	double GetDRAAdsorbedAmout(double bulk_pressure, double temperature, std::vector<double> potential_params, call_mono_eos eos, call_potential get_potential);
	double GetLJAdsorbedAmout(double Pb, double T, std::vector<double> potential_params, call_mono_eos eos, call_potential get_potential);
};

#endif
