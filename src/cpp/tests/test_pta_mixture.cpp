#include <iostream>
#include <gtest/gtest.h>
#include <vector>
#include <cassert>
#include "../src/pta_mixture.h"

TEST(test_pta_mixture, testMixtureCalculation)
{
    auto pta_model = MixturePTA("DRA", "pr77", "excess", 155);

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;

    Fluid ch4;
    ch4.CriticalPressure = 45.992e5;
    ch4.CriticalTemperature = 190.56;
    ch4.AccentricFactor = 0.01142;

    std::vector<double> CO2_DRA_params = {7880.19, 0.29, 2.};
    std::vector<double> CH4_DRA_params = {5600, 0.36, 3.};

    std::vector<double> pressures = {10000.0, 37000.0, 52000.0, 77000.0, 101000.0, 504000.0, 966000.0, 1989000.0, 2692000.0, 3930000.0, 4912000.0, 5753000.0};

    std::vector<double> calculated = {0.0834321, 0.438571, 0.623465, 0.909272, 1.15948, 3.61979, 4.93091, 6.24615, 6.6598, 6.9871, 7.0454, 7.0007};

    std::vector<double> bulk_composition_fraction = {0.25, 0.75};

    std::vector<double>
        ns = pta_model.GetLoading(bulk_composition_fraction, pressures[7], 318.2, {CO2_DRA_params, CH4_DRA_params}, {co2, ch4});

    ASSERT_EQ(bulk_composition_fraction.size(), ns.size());
}

TEST(test_pta_mixture, testMixtureCalculationThreeComponents)
{
    auto pta_model = MixturePTA("DRA", "pr77", "excess", 155);

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;

    Fluid ch4;
    ch4.CriticalPressure = 45.992e5;
    ch4.CriticalTemperature = 190.56;
    ch4.AccentricFactor = 0.01142;

    Fluid n2;
    n2.CriticalPressure = 33.958e5;
    n2.CriticalTemperature = 126.192;
    n2.AccentricFactor = 0.0372;

    std::vector<double> CO2_DRA_params = {7880.19, 0.29, 2.};
    std::vector<double> CH4_DRA_params = {7000.19, 0.2, 3.};
    std::vector<double> N2_DRA_params = {5000.5, 0.57, 2.};

    std::vector<double> pressures = {10000.0, 37000.0, 52000.0, 77000.0, 101000.0, 504000.0, 966000.0, 1989000.0, 2692000.0, 3930000.0, 4912000.0, 5753000.0};

    std::vector<double> calculated = {0.0834321, 0.438571, 0.623465, 0.909272, 1.15948, 3.61979, 4.93091, 6.24615, 6.6598, 6.9871, 7.0454, 7.0007};

    std::vector<double> bulk_composition_fraction = {0.25, 0.25, 0.50};

    std::vector<double>
        ns = pta_model.GetLoading(bulk_composition_fraction, pressures[11], 318.2, {CO2_DRA_params, CH4_DRA_params, N2_DRA_params}, {co2, ch4, n2});

    ASSERT_EQ(bulk_composition_fraction.size(), ns.size());
}
