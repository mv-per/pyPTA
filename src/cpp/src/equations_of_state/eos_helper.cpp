#include "math.h"
#include "eos_helper.h"

#include "eos.h"

const double R = 8.314462618;

std::vector<double>
find_z(double a0, double a1, double a2)
{

    double amax, rr, xinf1, xold, xnew, x1, x2, x3, gx0, gx1, gx2, a, b, c;
    std::size_t iter;

    amax = maxvalue(fabs(a2), fabs(a1), fabs(a0));
    rr = 1.0 + amax;
    xinf1 = -a2 / 3.0;

    if (gx(xinf1, a0, a1, a2) > 0.0)
    {
        xold = -rr;
    }
    else
    {
        xold = rr;
    }

    iter = 0;
    while (iter < 50)
    {
        gx0 = gx(xold, a0, a1, a2);
        gx1 = dgx(xold, a1, a2);
        gx2 = d2gx(xold, a2);
        xnew = xold - gx0 * gx1 / (gx1 * gx1 - 1.0 / 2.0 * gx0 * gx2); // kepler method

        if ((fabs(xnew - xold) < 1.0e-10) || (gx0 == 0))
        {
            break;
        }
        xold = xnew;
        iter++;
    }

    x1 = xnew;
    a = 1.0;
    b = a * x1 + a2;
    c = b * x1 + a1;

    double discrim = b * b - 4.0 * a * c;
    if (discrim < 0.0)
    {
        x2 = x1;
        x3 = x1;
    }
    else
    {
        x2 = (-b + sqrt(discrim)) / (2.0 * a);
        x3 = (-b - sqrt(discrim)) / (2.0 * a);
    }

    std::vector<double> results = {x1, x2, x3};

    return results;
}

double gx(double X, double a0, double a1, double a2)
{
    return X * X * X + a2 * X * X + a1 * X + a0;
}

double dgx(double X, double a1, double a2)
{
    return 3.0 * X * X + 2.0 * a2 * X + a1;
}

double d2gx(double X, double a2)
{
    return 6.0 * X + 2.0 * a2;
}

double minvalue(double num1, double num2, double num3)
{
    if ((num1 < num2) && (num1 < num3) && (fabs(num1) == num1) && (fabs(num2) == num2) && (fabs(num3) == num3))
    {
        return num1;
    }
    else if ((num2 < num1) && (num2 < num3) && (fabs(num1) == num1) && (fabs(num2) == num2) && (fabs(num3) == num3))
    {
        return num2;
    }
    else
    {
        return num3;
    }
}

double maxvalue(double num1, double num2, double num3)
{
    if ((num1 > num2) && (num1 > num3))
    {
        return num1;
    }
    else if ((num2 > num1) && (num2 > num3))
    {
        return num2;
    }
    else
    {
        return num3;
    }
}

void CheckValidPressure(double P)
{
    // TODO: improve out of range values in EoS
    // if (P < 0)
    // {
    //     std::cout << "Pressure values cannot be negative, P=" << P << std::endl;
    // }
}

/**
* @brief Calculates the temperature-dependent volume shift factor for a fluid.

* This function calculates the volume shift factor for a fluid at a given temperature
* using the Peng-Robinson equation of state. The volume shift factor is used to adjust
* the ideal gas law to account for deviations from ideal behavior in real gases.
*
* Source: Ahlers, J., & Gmehling, J. (2001). 
* Development of an universal group contribution equation of state.
* Fluid Phase Equilibria, 191(1-2), 177–188. doi:10.1016/s0378-3812(01)00626-4 
* 
* @param T Temperature in Kelvin.
* @param fluid Fluid object containing the critical properties of the fluid.
* @return The temperature-dependent volume shift factor.
*/
double GetTemperatureDependentVolumeShiftFactor(double T, Fluid fluid){
    double Tr = T/fluid.CriticalTemperature;
    double Cci = R*fluid.CriticalTemperature*(0.3074-fluid.CriticalCompressibility)/fluid.CriticalPressure; // [m3/mol]
    double Ci = 0.252*R*fluid.CriticalTemperature*(-0.4024 + 1.5448*fluid.CriticalCompressibility) / fluid.CriticalPressure;  // [m3/mol]
    double gamma = 12.67 - 107.21*fluid.CriticalCompressibility + 246.78*fluid.CriticalCompressibility*fluid.CriticalCompressibility; // [-]
    double eta = 26.966-74.458*fluid.CriticalCompressibility; // [-]
    double m = 0.37464+1.54226*fluid.AccentricFactor - 0.26992*fluid.AccentricFactor*fluid.AccentricFactor; // [-]
    double alpha = std::pow(1 + m*(1-std::sqrt(Tr)),2); // [-]

    return Ci - (0.35*Cci / (0.35 + std::pow(eta * std::fabs(alpha - Tr), gamma)));  // [m3/mol]
}

/** UNUSED
 * @brief Calculates the volume shift factor for a fluid.
 * 
 * This function calculates the volume shift factor for a fluid using the Peng-Robinson equation of state.
 * The volume shift factor is used to adjust the ideal gas law to account for deviations from ideal behavior
 * in real gases.
 * 
 * @param fluid Fluid object containing the critical properties of the fluid.
 * @return The fluid volume shift factor [m^3/mol].
 */
double GetFluidVolumeShiftFactor(Fluid fluid)
{
    return 0.40768 *
           R * fluid.CriticalTemperature * (0.29441 - fluid.CriticalCompressibility) / fluid.CriticalPressure;
           // [m3 Pa / K / mol * K / Pa] == [m3/mol]
}

double CalculatePurePenelouxVolumeTranslation(double vol, double T, Fluid fluid)
{
    return vol - GetTemperatureDependentVolumeShiftFactor(T, fluid); // [m3/mol] - [m3/mol]
}

double CalculateMixturePenelouxVolumeTranslation(double vol, double T, std::vector<double> molar_fractions, std::vector<Fluid> fluids)
{

    double MixtureVolumeShiftFactor = 0;
    for (std::size_t i = 0; i < molar_fractions.size(); i++)
    {
        MixtureVolumeShiftFactor += GetTemperatureDependentVolumeShiftFactor(T, fluids[i]) * molar_fractions[i];
    }

    return vol - MixtureVolumeShiftFactor;
}

Fluid CheckForFluidCriticalCompressibility(Fluid fluid){
    if (!fluid.CriticalCompressibility){
        std::cout << "Critical compressibility not provided for " << fluid.Name << ", Calculating from fluid's accentric factor" << std::endl;
        fluid.CriticalCompressibility  = 0.29506 - 0.08775* fluid.AccentricFactor;
    }

    return fluid;
}