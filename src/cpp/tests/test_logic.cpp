
#include <gtest/gtest.h>
#include <functional>
#include "../src/equations_of_state/eos.h"

#include <iostream>

typedef std::function<mono_eos(double, double)> call_eos;

std::function<mono_eos(double, double)> GetEquationOfStateInvoker(Fluid fluid)
{
    return [=](double P, double T)
    {
        return pr77().get_mono_fluid_properties(P, T, fluid);
    };
}

TEST(test_logic, DemonstrateEOSPointerLogic)
{

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;

    call_eos eos_pointer = GetEquationOfStateInvoker(co2);

    double fugacity = eos_pointer(1000000, 303).fug;

    EXPECT_NEAR(fugacity, 948684.4705, 1e-4);
}

TEST(test_logic, CalculateSRKMono)
{

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;

    double fugacity = srk().get_mono_fluid_properties(1000000, 303, co2).fug;

    EXPECT_NEAR(fugacity, 953627.1636, 1e-4);
}

TEST(test_logic, CalculatePR77PENELOUXMonoWithoutProvidingCriticalZ)
{

    Fluid co2;
    co2.CriticalPressure = 7.38e6;
    co2.CriticalTemperature = 304.19;
    co2.AccentricFactor = 0.22394;

    double density = pr77_peneloux().get_mono_fluid_properties(1000000, 303, co2).dens;

    EXPECT_NEAR(density, 415.80276, 1e-4);
}

TEST(test_logic, CalculatePR77PENELOUXMonoProvidingCriticalZ)
{

    Fluid co2;
    co2.CriticalPressure = 7.38e6;
    co2.CriticalTemperature = 304.19;
    co2.AccentricFactor = 0.22394;
    co2.CriticalCompressibility = 0.29;

    double density = pr77_peneloux().get_mono_fluid_properties(1000000, 303, co2).dens;


    EXPECT_NEAR(density, 418.524068, 1e-4);

}

TEST(test_logic, CalculateSRKMix)
{

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;

    Fluid ch4;
    ch4.CriticalPressure = 45.992e5;
    ch4.CriticalTemperature = 190.56;
    ch4.AccentricFactor = 0.01142;

    mix_eos props = srk().get_mixture_fluid_properties({0.5, 0.5}, 1000000, 303, {co2, ch4});

    double co2_fug = props.fug[0];
    double ch4_fug = props.fug[1];

    EXPECT_NEAR(co2_fug, 412439.4438, 1e-4);
    EXPECT_NEAR(ch4_fug, 460909.2332, 1e-4);
}


TEST(test_logic, CalculateSRKMixPeneloux)
{

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;

    Fluid ch4;
    ch4.CriticalPressure = 45.992e5;
    ch4.CriticalTemperature = 190.56;
    ch4.AccentricFactor = 0.01142;

    mix_eos props = srk_peneloux().get_mixture_fluid_properties({0.5, 0.5}, 1000000, 303, {co2, ch4});

    double co2_fug = props.fug[0];
    double ch4_fug = props.fug[1];

    EXPECT_NEAR(co2_fug, 412439.4438, 1e-4);
    EXPECT_NEAR(ch4_fug, 460909.2332, 1e-4);
}
