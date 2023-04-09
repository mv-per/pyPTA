#ifndef PR77_PENELOUX_H
#define PR77_PENELOUX_H

#include <tuple>
#include <cmath>
#include <vector>
#include "../data_classes.h"
#include "../helpers.h"
#include "eos_helper.h"

class pr77_peneloux
{

	std::size_t i, j;
	std::size_t ncomp;
	double Zmin, Zvalue, Zmax;
	double gibbsenergymin, gibbsenergymax, amix, bmix, A, B, sumat;
	double Z, dens, vol, Z_;
	std::tuple<double, double> get_critical_properties(double T, double Pc, double Tc, double w);
	std::tuple<double, double, double, double> get_gibbs_energy(double P, std::vector<double> Z, double A, double B, int min_or_max);

public:
	struct mono_eos get_mono_fluid_properties(double P, double T, Fluid fluid);
	struct mix_eos get_mixture_fluid_properties(std::vector<double> x, double P, double T, std::vector<Fluid> fluids);
};

#endif