#ifndef EOS_HELPER_HPP
#define EOS_HELPER_HPP

#include <string>
#include <vector>

/// Gas contsant, [m3 Pa /K / mol]
extern double R;

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

std::vector<double> find_z(double A, double B);
double gx(double X, double a0, double a1, double a2);
double dgx(double X, double a1, double a2);
double d2gx(double X, double a2);
double minvalue(double num1, double num2, double num3);
double maxvalue(double num1, double num2, double num3);

#endif
