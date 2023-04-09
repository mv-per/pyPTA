#ifndef PR77_PENELOUX_H
#define PR77_PENELOUX_H

#include <tuple>
#include <cmath>
#include <vector>
#include "../data_classes.h"
#include "../helpers.h"
#include "eos_helper.h"
#include "pr77.h"

class pr77_peneloux
{

public:
	struct mono_eos get_mono_fluid_properties(double P, double T, Fluid fluid);
	struct mix_eos get_mixture_fluid_properties(std::vector<double> x, double P, double T, std::vector<Fluid> fluids);
};

#endif