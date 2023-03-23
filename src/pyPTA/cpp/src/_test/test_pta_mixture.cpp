#include <iostream>
#include <gtest/gtest.h>
#include <vector>
#include <cassert>
#include "../pta_mixture.h"

TEST(test_pta_mixture, testMixtureCalculation)
{
    auto pta_model = MixturePTA("DRA", "pr77", "excess", 155);

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;
    co2.LennardJonnesDiameter = 3.941;

    std::vector<double> DRA_params = {7880.19, 0.29, 2.};

    std::vector<double> pressures = {10000.0, 37000.0, 52000.0, 77000.0, 101000.0, 504000.0, 966000.0, 1989000.0, 2692000.0, 3930000.0, 4912000.0, 5753000.0};

    std::vector<double> calculated = {0.0834321, 0.438571, 0.623465, 0.909272, 1.15948, 3.61979, 4.93091, 6.24615, 6.6598, 6.9871, 7.0454, 7.0007};

    std::vector<double>
        ns = pta_model.GetLoading({0.25, 0.25, 0.50}, pressures[7], 318.2, {DRA_params, DRA_params, DRA_params}, {co2, co2, co2});

    double diff;
    for (std::size_t i = 0; i < ns.size(); i++)
    {
        // std::cout << ns[i] << std::endl;
        std::cout << ns[i] << std::endl;
    }
}
