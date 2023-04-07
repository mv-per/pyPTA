#include "helpers.h"
#include "pta_helper.h"
#include "equations_of_state/eos.h"

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

double integrate(std::vector<double> integral_vector, double step)
{
    // if (this->integration_method == "Simpson")
    // {
    // 	return SIMPS(integral_vector, step, num_of_layers);
    // }

    // At this time, only simpson rule is implemented
    return simpson_rule(integral_vector, step);
}

double GetCalculatedAdsorbedDensity(double adsorbed_density, double bulk_density, std::string isotherm_type)
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

double GetAdsorptionPotentialEnergy(double eps, double T, std::string potential)
{
    if (potential == DRA_POTENTIAL)
    {
        return std::exp(eps / (8.314 * T));
    }
    else if (potential == STEELE_POTENTIAL || potential == LEE_POTENTIAL)
    {
        return std::exp(-eps / (1.38064852e-23 * T));
    }
    else
    {
        throw std::invalid_argument("Isotherm potential not found/defined.");
    }
}

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