
#include <gtest/gtest.h>
#include <functional>
#include "../src/equations_of_state/pr77.h"

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

    EXPECT_NEAR(fugacity, 948684.4705, 1e-2);
}
