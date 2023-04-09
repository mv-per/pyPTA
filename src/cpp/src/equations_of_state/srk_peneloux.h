#ifndef SRK_PENELOUX_H
#define SRK_PENELOUX_H

#include "srk.h"

class srk_peneloux
{

public:
	struct mono_eos get_mono_fluid_properties(double P, double T, Fluid fluid);
	struct mix_eos get_mixture_fluid_properties(std::vector<double> x, double P, double T, std::vector<Fluid> fluids);
};

#endif