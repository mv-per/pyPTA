

#include "srk_peneloux.h"

/**
 * Get the properties of a mono-component fluid.
 *
 * @param P Pressure of the fluid.
 * @param T Temperature of the fluid.
 * @param fluid fluid properties.
 * @return Struct containing the fluid properties.
 */
struct mono_eos srk_peneloux::get_mono_fluid_properties(double P, double T, Fluid fluid)
{
    double VolumeShiftFactor = GetTemperatureDependentVolumeShiftFactor(T, fluid);
    return srk(VolumeShiftFactor).get_mono_fluid_properties(P, T, fluid);
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
struct mix_eos srk_peneloux::get_mixture_fluid_properties(std::vector<double> x, double P, double T, std::vector<Fluid> fluids)
{
    double MixtureVolumeShiftFactor = 0;
    for (std::size_t i = 0; i < x.size(); i++)
    {
        MixtureVolumeShiftFactor += GetTemperatureDependentVolumeShiftFactor(T, fluids[i]) * x[i];
    }
    return srk(MixtureVolumeShiftFactor).get_mixture_fluid_properties(x, P, T, fluids);
}
