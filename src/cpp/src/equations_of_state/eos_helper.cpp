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

double GetFluidVolumeShiftFactor(Fluid fluid)
{
    double CriticalCompressibility;
    if (fluid.CriticalCompressibility)
    {
        CriticalCompressibility = *fluid.CriticalCompressibility;
    }
    else
    {
        std::cout << "Critical compressibility not provided, calculating from critical properties" << std::endl;

        // double CriticalVolume = R * fluid.CriticalTemperature / fluid.CriticalPressure;

        // std::cout << CriticalVolume << std::endl;
        // std::cout << fluid.CriticalPressure << std::endl;
        // std::cout << fluid.Criti << std::endl;

        // CriticalCompressibility = (fluid.CriticalPressure * CriticalVolume) / (R * fluid.CriticalTemperature);
        CriticalCompressibility = pr77().get_mono_fluid_properties(fluid.CriticalPressure, fluid.CriticalTemperature, fluid).Z;
    }

    std::cout << R << std::endl;
    std::cout << CriticalCompressibility << std::endl;

    return 0.40768 *
           R * fluid.CriticalTemperature * (0.29441 - CriticalCompressibility) / fluid.CriticalPressure;
}

double CalculatePurePenelouxVolumeTranslation(double vol, Fluid fluid)
{
    return vol + GetFluidVolumeShiftFactor(fluid) / *fluid.MolecularWeight;
}

double CalculateMixturePenelouxVolumeTranslation(double vol, std::vector<double> molar_fractions, std::vector<Fluid> fluids)
{

    double MixtureMolarWeight = 0;
    double MixtureVolumeShiftFactor = 0;
    for (std::size_t i = 0; i < molar_fractions.size(); i++)
    {
        MixtureMolarWeight += *fluids[i].MolecularWeight * molar_fractions[i];
        MixtureVolumeShiftFactor += GetFluidVolumeShiftFactor(fluids[i]) * molar_fractions[i];
    }

    return vol + MixtureVolumeShiftFactor / MixtureMolarWeight;
}