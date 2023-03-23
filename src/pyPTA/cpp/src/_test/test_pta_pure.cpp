#include <gtest/gtest.h>
#include <vector>
#include "../pta_pure.h"

TEST(test_pta_pure, TestGetMultiplePurePTA)
{
    PurePTA pta_model = PurePTA("DRA", "pr77", "excess", 155);

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;
    co2.LennardJonnesDiameter = 3.941;

    std::vector<double> CO2DRAParams = {7880.19, 0.29, 2.};

    std::vector<double> pressures = {10000.0, 37000.0, 52000.0, 77000.0, 101000.0, 504000.0, 966000.0, 1989000.0, 2692000.0, 3930000.0, 4912000.0, 5753000.0};

    std::vector<double> ExpectedLoadings = {0.0834321, 0.438571, 0.623465, 0.909272, 1.15948, 3.61979, 4.93091, 6.24615, 6.6598, 6.9871, 7.0454, 7.0007};

    std::vector<double>
        CalculatedLoadings = pta_model.GetMultipleLoadings(pressures, 318.2, CO2DRAParams, co2);

    // ASSERT_EQ(pressures.size(), ExpectedLoadings.size());
    for (std::size_t i = 0; i < CalculatedLoadings.size(); i++)
    {
        EXPECT_NEAR(ExpectedLoadings[i], CalculatedLoadings[i], 1e-2);
        // std::cout << CalculatedLoadings[i] << std::endl;
    }
};

TEST(test_pta_pure, TestSinglePTALoading)
{

    auto pta_model = PurePTA("DRA", "pr77", "excess", 100);

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;

    std::vector<double> DRA_params = {7880.19, 0.29, 2.};

    double adsorbed_loading = pta_model.GetLoading(1.e6, 303., DRA_params, co2);

    ASSERT_EQ(adsorbed_loading, 5.8862833281894718);
};

TEST(test_pta_pure, TestGetMultiplePureDeviation)
{
    PurePTA pta_model = PurePTA("DRA", "pr77", "excess", 155);

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;
    co2.LennardJonnesDiameter = 3.941;

    std::vector<double> CO2DRAParams = {7880.19, 0.29, 2.};

    std::vector<double> Pressures = {10000.0, 37000.0, 52000.0, 77000.0, 101000.0, 504000.0, 966000.0, 1989000.0, 2692000.0, 3930000.0, 4912000.0, 5753000.0};

    std::vector<double> ExperimentalLoadings = {0.0834321, 0.438571, 0.623465, 0.909272, 1.15948, 3.61979, 4.93091, 6.24615, 6.6598, 6.9871, 7.0454, 7.0007};

    double CalculatedDeviation = pta_model.GetDeviationRange("absolute", Pressures, ExperimentalLoadings, 318.2, CO2DRAParams, co2);

    ASSERT_EQ(Pressures.size(), ExperimentalLoadings.size());

    EXPECT_EQ(CalculatedDeviation, 0.0);
    // for (std::size_t i = 0; i < CalculatedDeviation.size(); i++)
    // {
    //     // EXPECT_EQ(ExpectedLoadings[i], CalculatedDeviation[i]);
    //     // std::cout << CalculatedDeviation[i] << std::endl;
    // }
};