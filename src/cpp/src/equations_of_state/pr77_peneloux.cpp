

#include "pr77_peneloux.h"

/**
 * Original Peng-Robinson EOS
 *
 * Peng, D., & Robinson, D. B. (1977).
 * A rigorous method for predicting the critical properties of multicomponent systems from an equation of state.
 * AIChE Journal, 23(2), 137–144. doi:10.1002/aic.690230202
 *
 */

/**
 * Get the properties of a mono-component fluid.
 *
 * @param P Pressure of the fluid.
 * @param T Temperature of the fluid.
 * @param fluid fluid properties.
 * @return Struct containing the fluid properties.
 */
struct mono_eos pr77_peneloux::get_mono_fluid_properties(double P, double T, Fluid fluid)
{

    double VolumeShiftFactor = GetTemperatureDependentVolumeShiftFactor(T, fluid);
    return pr77(VolumeShiftFactor).get_mono_fluid_properties(P, T, fluid);
}

/**
 * Get the properties of a mixture.
 *
 * @param x Composition fractions for each fluid in the mixture.
 * @param P Pressure of the fluid.
 * @param T Temperature of the fluid.
 * @param fluids List fluids on the mixture.
 * @return Struct containing the fluid properties
 */
struct mix_eos pr77_peneloux::get_mixture_fluid_properties(std::vector<double> x, double P, double T, std::vector<Fluid> fluids)
{
    double MixtureVolumeShiftFactor = 0;
    for (std::size_t i = 0; i < x.size(); i++)
    {
        MixtureVolumeShiftFactor += GetTemperatureDependentVolumeShiftFactor(T, fluids[i]) * x[i];
    }
    return pr77(MixtureVolumeShiftFactor).get_mixture_fluid_properties(x, P, T, fluids);
}
