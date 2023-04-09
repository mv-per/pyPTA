#ifndef SRK_H
#define SRK_H

#include <tuple>
#include <cmath>
#include <vector>
#include "../data_classes.h"
#include "../helpers.h"
#include "eos_helper.h"

class srk
{

	std::size_t i, j;
	std::size_t ncomp;
	double Zmin, Zvalue, Zmax;
	double gibbsenergymin, gibbsenergymax, amix, bmix, A, B, sumat;
	double Z, dens, vol, Z_;
	double VolumeShiftFactor = 0;
	std::tuple<double, double> get_critical_properties(double T, double Pc, double Tc, double w);
	std::tuple<double, double, double, double> get_gibbs_energy(double P, std::vector<double> Z, double A, double B, int min_or_max);

public:
	srk();
	srk(double VolumeShiftFactor);
	struct mono_eos get_mono_fluid_properties(double P, double T, Fluid fluid);
	struct mix_eos get_mixture_fluid_properties(std::vector<double> x, double P, double T, std::vector<Fluid> fluids);
};

#endif