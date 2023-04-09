#ifndef EOS_HELPER_HPP
#define EOS_HELPER_HPP

#include <string>
#include <vector>
#include <iostream>
#include "../data_classes.h"

/// Gas contsant, [m3 Pa /K / mol]
extern const double R;

struct mono_eos
{
	double fug, dens, phi, Z;
};

struct mix_eos
{
	std::vector<double> fug;
	double dens;
	std::vector<double> phi;
	double Z;
};

std::vector<double> find_z(double a0, double a1, double a2);
double gx(double X, double a0, double a1, double a2);
double dgx(double X, double a1, double a2);
double d2gx(double X, double a2);
double minvalue(double num1, double num2, double num3);
double maxvalue(double num1, double num2, double num3);

void CheckValidPressure(double P);

double CalculateMixturePenelouxVolumeTranslation(double vol, double T, std::vector<double> molar_fractions, std::vector<Fluid> fluids);
double CalculatePurePenelouxVolumeTranslation(double vol, double T, Fluid fluid);

double GetTemperatureDependentVolumeShiftFactor(double T, Fluid fluid);
Fluid CheckForFluidCriticalCompressibility(Fluid fluid);

#endif
