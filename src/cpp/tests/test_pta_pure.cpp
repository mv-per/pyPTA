#include <gtest/gtest.h>
#include <vector>
#include "../src/pta_pure.h"

TEST(test_pta_pure, TestSinglePTALoading)
{

    std::string potential = DRA_POTENTIAL;
    auto pta_model = PurePTA(potential, "pr77", "excess", 100);

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;

    std::vector<double> DRA_params = {7880.19, 0.29, 2.};

    double adsorbed_loading = pta_model.GetLoading(1.e6, 303., DRA_params, co2);

    EXPECT_NEAR(adsorbed_loading, 5.8862, 1e-4);
};

TEST(test_pta_pure, TestSinglePTALoading2)
{
    std::string potential = DRA_POTENTIAL;
    auto pta_model = PurePTA(potential, "pr77", "excess", 100);

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;

    std::vector<double> DRA_params = {7880.19, 0.29, 2.};

    double adsorbed_loading = pta_model.GetLoading(10.e6, 323., DRA_params, co2);

    EXPECT_NEAR(adsorbed_loading, 5.62658, 1e-4);
};

TEST(test_pta_pure, TestGetMultiplePurePTA)
{
    std::string potential = DRA_POTENTIAL;
    PurePTA pta_model = PurePTA(potential, "pr77", "excess", 155);

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;
    co2.LennardJonnesDiameter = 3.941;

    std::vector<double> CO2DRAParams = {7880.19, 0.29, 2.};

    std::vector<double> pressures = {10000.0, 37000.0, 52000.0, 77000.0, 101000.0, 504000.0, 966000.0, 1989000.0, 2692000.0, 3930000.0, 4912000.0, 5753000.0};

    std::vector<double> ExpectedLoadings = {0.0834321, 0.438571, 0.623465, 0.909272, 1.15948, 3.61979, 4.93091, 6.24615, 6.6598, 6.9871, 7.0454, 7.0007};

    EXPECT_NEAR(0.0001, 0.0001, 0.2);
    std::vector<double>
        CalculatedLoadings = pta_model.GetLoadings(pressures, 318.2, CO2DRAParams, co2);

    ASSERT_EQ(pressures.size(), ExpectedLoadings.size());
    for (std::size_t i = 0; i < CalculatedLoadings.size(); i++)
    {
        EXPECT_NEAR(ExpectedLoadings[i], CalculatedLoadings[i], 1e-2);
        ;
    }
};

TEST(test_pta_pure, TestGetMultiplePureDeviationDRA)
{
    std::string potential = DRA_POTENTIAL;
    PurePTA pta_model = PurePTA(potential, "pr77", "excess", 155);

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;
    co2.LennardJonnesDiameter = 3.941;

    std::vector<double> CO2DRAParams = {7880.19, 0.29, 2.};

    std::vector<double> pressures = {10000.0, 37000.0, 52000.0, 77000.0, 101000.0, 504000.0, 966000.0, 1989000.0, 2692000.0, 3930000.0, 4912000.0, 5753000.0};

    std::vector<double> experimental = {0.0834321, 0.438571, 0.623465, 0.909272, 1.15948, 3.61979, 4.93091, 6.24615, 6.6598, 6.9871, 7.0454, 7.0007};

    double CalculatedDeviation = pta_model.GetDeviationRange("absolute", experimental, pressures, 318.2, CO2DRAParams, co2);

    ASSERT_EQ(pressures.size(), experimental.size());

    EXPECT_NEAR(CalculatedDeviation, 0.0, 1e-4);
};

TEST(test_pta_pure, TestGetMultiplePureDeviationLEE)
{
    std::string potential = LEE_POTENTIAL;
    PurePTA pta_model = PurePTA(potential, "pr77", "excess", 155);

    Adsorbent Z01x = Adsorbent("Z01x", 3.35, 0.382);

    pta_model.SetAdsorbent(Z01x);

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;
    co2.LennardJonnesDiameter = 3.941;

    std::vector<double> CO2LEEParams = {125.63, 12.26, 765.70}; // eps/K, L, A

    std::vector<double> pressures = {10000.0, 37000.0, 52000.0, 77000.0, 101000.0, 504000.0, 966000.0, 1989000.0, 2692000.0, 3930000.0, 4912000.0, 5753000.0};

    std::vector<double> experimental = {0.0834321, 0.438571, 0.623465, 0.909272, 1.15948, 3.61979, 4.93091, 6.24615, 6.6598, 6.9871, 7.0454, 7.0007};

    double CalculatedDeviation = pta_model.GetDeviationRange("absolute", experimental, pressures, 318.2, CO2LEEParams, co2);

    ASSERT_EQ(pressures.size(), experimental.size());

    EXPECT_NEAR(CalculatedDeviation, 2.5183, 1e-4);
};

TEST(test_pta_pure, TestGetMultiplePureDeviationSTEELE)
{
    std::string potential = STEELE_POTENTIAL;
    PurePTA pta_model = PurePTA(potential, "pr77", "excess", 155);

    Adsorbent adsorbent = Adsorbent("Z01x", 3.35, 0.382);

    pta_model.SetAdsorbent(adsorbent);

    Fluid fluid;
    fluid.CriticalPressure = 73.773e5;
    fluid.CriticalTemperature = 304.13;
    fluid.AccentricFactor = 0.22394;
    fluid.LennardJonnesDiameter = 3.941;

    std::vector<double> params = {109.32, 13.34, 611.88}; // eps/K, L, A

    std::vector<double> pressures = {10000.0, 37000.0, 52000.0, 77000.0, 101000.0, 504000.0, 966000.0, 1989000.0, 2692000.0, 3930000.0, 4912000.0, 5753000.0};

    std::vector<double> experimental = {0.0834321, 0.438571, 0.623465, 0.909272, 1.15948, 3.61979, 4.93091, 6.24615, 6.6598, 6.9871, 7.0454, 7.0007};

    double CalculatedDeviation = pta_model.GetDeviationRange("absolute", experimental, pressures, 318.2, params, fluid);

    ASSERT_EQ(pressures.size(), experimental.size());

    EXPECT_NEAR(CalculatedDeviation, 41.7384, 1e-4);
};

// equations of state

TEST(test_pta_pure, TestGetMultiplePureDeviationDRAWithSRK)
{
    std::string potential = DRA_POTENTIAL;
    PurePTA pta_model = PurePTA(potential, "srk", "excess", 155);

    Fluid co2;
    co2.CriticalPressure = 73.773e5;
    co2.CriticalTemperature = 304.13;
    co2.AccentricFactor = 0.22394;
    co2.LennardJonnesDiameter = 3.941;

    std::vector<double> CO2DRAParams = {7880.19, 0.29, 2.};

    std::vector<double> pressures = {10000.0, 37000.0, 52000.0, 77000.0, 101000.0, 504000.0, 966000.0, 1989000.0, 2692000.0, 3930000.0, 4912000.0, 5753000.0};

    std::vector<double> expected = {0.076485, 0.38631, 0.546913, 0.795782, 1.013569, 3.157997, 4.30760, 5.46686, 5.83229, 6.119404, 6.16664, 6.12293};

    std::vector<double> calculated = pta_model.GetLoadings(pressures, 318.2, CO2DRAParams, co2);

    ASSERT_EQ(pressures.size(), calculated.size());

    for (std::size_t i = 0; i < pressures.size(); i++)
    {
        EXPECT_NEAR(calculated[i], expected[i], 1e-4);
    }
};