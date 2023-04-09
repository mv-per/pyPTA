

#include "pr77.h"

/**
 * Original Peng-Robinson EOS
 *
 * Peng, D., & Robinson, D. B. (1977).
 * A rigorous method for predicting the critical properties of multicomponent systems from an equation of state.
 * AIChE Journal, 23(2), 137–144. doi:10.1002/aic.690230202
 *
 */
pr77::pr77(){}

pr77::pr77(double VolumeShiftFactor){
    this->VolumeShiftFactor = VolumeShiftFactor;
}

std::tuple<double, double,double> pr77::GetQuadraticCoefficients(double A, double B){
    double a2 = B - 1.0;
    double a1 = A - 3.0 * B * B - 2.0 * B;
    double a0 = -A * B + B * B + B * B * B;
    return std::make_tuple(a0, a1, a2);
}

/**
 * Get the Peng-Robinson critical parameters (a,b)
 *
 * @param T Temperature of the fluid.
 * @param Pc Critical Pressure of the Fluid.
 * @param Tc Critical Temperature of the fluid.
 * @param w Acentric Factor of the fluid.
 * @return A Tuple containin a and b.
 */
std::tuple<double, double> pr77::get_critical_properties(double T, double Pc, double Tc, double w)
{
    double Tr = T / Tc;
    double kappa = 0.37464 + 1.54226 * w - 0.26992 * w * w;
    double alpha = std::pow(1 + kappa * (1 - sqrt(Tr)), 2.0);
    double a = 0.457235 * R * R * Tc * Tc * alpha / Pc;
    double b = 0.077796 * R * Tc / Pc;
    return std::make_tuple(a, b);
}

/**
 * Get the Gibbs energy from a compressibility factor
 *
 * @param P Pressure of the fluid.
 * @param Z List of compressibility factors.
 * @param A Non-dimensional A parameter
 * @param B Non-dimensional B parameter
 * @param min_or_max Flag to check the minimal (= 0) or maximum (=1) compressibility factor of the list
 * @return A Tuple containing the gibbs energy, compressibility factor, fugacity coefficient and fugacity.
 */
std::tuple<double, double, double, double> pr77::get_gibbs_energy(double P, std::vector<double> Z, double A, double B, int min_or_max)
{
    double phi, fug;
    if (min_or_max == 0)
    {
        Z_ = minvalue(Z[0], Z[1], Z[2]);
    }

    else if (min_or_max == 1)
    {
        Z_ = maxvalue(Z[0], Z[1], Z[2]);
    }
    phi = exp((Z_ - 1.0) - log(Z_ - B) - A / B / 2.8284 * log((Z_ + 2.414 * B) / (Z_ - 0.414 * B)));
    fug = phi * P;
    return std::make_tuple(log(fug), Z_, phi, fug);
}

/**
 * Get the properties of a mono-component fluid.
 *
 * @param P Pressure of the fluid.
 * @param T Temperature of the fluid.
 * @param fluid fluid properties.
 * @return Struct containing the fluid properties.
 */
struct mono_eos pr77::get_mono_fluid_properties(double P, double T, Fluid fluid)
{
    CheckValidPressure(P);

    double phi, fug, phi_min, phi_max, fug_min, fug_max;
    double a, b;

    std::tie(a, b) = get_critical_properties(T, fluid.CriticalPressure, fluid.CriticalTemperature, fluid.AccentricFactor);

    if (this->VolumeShiftFactor){
        b = b-this->VolumeShiftFactor;
    }

    double A = a * P / R / R / T / T;
    double B = b * P / R / T;

    // cubic equation parameters
    std::tie(a0, a1, a2) = GetQuadraticCoefficients(A, B);

    std::vector<double> Z = find_z(a0, a1, a2);

    // Calculate the proper Minimal Gibbs Energy
    std::tie(gibbsenergymin, Zmin, phi_min, fug_min) = get_gibbs_energy(P, Z, A, B, 0);
    std::tie(gibbsenergymax, Zmax, phi_max, fug_max) = get_gibbs_energy(P, Z, A, B, 1);

    // Select the proper values
    if (gibbsenergymin < gibbsenergymax)
    {
        Zvalue = Zmin;
        phi = phi_min;
        fug = fug_min;
    }
    else
    {
        Zvalue = Zmax;
        phi = phi_max;
        fug = fug_max;
    }
    double vol = Zvalue * R * T / P;
    if (this->VolumeShiftFactor){
        vol = vol - this->VolumeShiftFactor;
    }
    
    double dens = 1.0 / vol;

    struct mono_eos results = {fug, dens, phi, Zvalue};
    return results;
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
struct mix_eos pr77::get_mixture_fluid_properties(std::vector<double> x, double P, double T, std::vector<Fluid> fluids)
{
    CheckValidPressure(P);

    ncomp = x.size();

    
    std::vector<double> fug(ncomp, 0.0);
    std::vector<double> phi(ncomp, 0.0);
    std::vector<double> phimin(ncomp, 0.0);
    std::vector<double> phimax(ncomp, 0.0);
    std::vector<double> fugmax(ncomp, 0.0);
    std::vector<double> fugmin(ncomp, 0.0);
    std::vector<double> Tr(ncomp, 0.0);
    std::vector<double> alpha(ncomp, 0.0);
    std::vector<double> kappa(ncomp, 0.0);
    std::vector<double> a(ncomp, 0.0);
    std::vector<double> b(ncomp, 0.0);

    std::vector<std::vector<double>> aij(ncomp, std::vector<double>(ncomp, 0.0));
    std::vector<std::vector<double>> bij(ncomp, std::vector<double>(ncomp, 0.0));
    std::vector<std::vector<double>> cij(ncomp, std::vector<double>(ncomp, 0.0));
    std::vector<std::vector<double>> dij(ncomp, std::vector<double>(ncomp, 0.0));

    for (i = 0; i < ncomp; i++)
    {
        std::tie(a[i], b[i]) = get_critical_properties(T, fluids[i].CriticalPressure, fluids[i].CriticalTemperature, fluids[i].AccentricFactor);
    }

    // Mixing-rules
    for (i = 0; i < ncomp; i++)
    {
        for (j = 0; j < ncomp; j++)
        {
            aij[i][j] = sqrt(a[i] * a[j]) * (1.0 - cij[i][j]);
            bij[i][j] = (b[i] + b[j]) / 2.0 * (1.0 + dij[i][j]);
        }
    }

    for (i = 0; i < ncomp; i++)
    {
        for (j = 0; j < ncomp; j++)
        {
            amix += x[i] * x[j] * aij[i][j];
            bmix += x[i] * x[j] * bij[i][j];
        }
    }

    A = amix * P / R / R / T / T;
    B = bmix * P / R / T;

    std::tie(a0, a1, a2) = GetQuadraticCoefficients(A, B);

    std::vector<double> Z = find_z(a0, a1, a2);

    Zmax = maxvalue(Z[0], Z[1], Z[2]);
    Zmin = minvalue(Z[0], Z[1], Z[2]);

    for (i = 0; i < ncomp; i++)
    {
        sumat = 0.;
        for (j = 0; j < ncomp; j++)
        {
            sumat += x[j] * aij[i][j];
        }
        phimin[i] = exp(b[i] / bmix * (Zmin - 1.0) - log(Zmin - B) + ((-A) / B / 2.8284) * (2.0 * sumat / amix - b[i] / bmix) * log((Zmin + 2.414 * B) / (Zmin - 0.414 * B)));
        phimax[i] = exp(b[i] / bmix * (Zmax - 1.0) - log(Zmax - B) + ((-A) / B / 2.8284) * (2.0 * sumat / amix - b[i] / bmix) * log((Zmax + 2.414 * B) / (Zmax - 0.414 * B)));
        fugmin[i] = phimin[i] * P * x[i];
        fugmax[i] = phimax[i] * P * x[i];
    }

    for (i = 0; i < ncomp; i++)
    {
        gibbsenergymin += x[i] * log(fugmin[i]);
        gibbsenergymax += x[i] * log(fugmax[i]);
    }

    if (gibbsenergymin < gibbsenergymax)
    {
        Zvalue = Zmin;
        for (i = 0; i < ncomp; i++)
        {
            phi[i] = phimin[i];
            fug[i] = fugmin[i];
        }
    }
    else
    {
        Zvalue = Zmax;
        for (i = 0; i < ncomp; i++)
        {
            phi[i] = phimax[i];
            fug[i] = fugmax[i];
        }
    }
    vol = Zvalue * R * T / P;
    dens = 1.0 / vol;

    struct mix_eos results = {fug, dens, phi, Zvalue};

    return results;
}
