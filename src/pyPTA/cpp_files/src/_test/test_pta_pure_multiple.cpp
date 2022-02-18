#include <iostream>
#include <vector>
#include <cassert>
#include "../pta_pure.h"

int main(void)
{
    auto pta_model = PurePTA("DRA", "pr77", "excess", 100);

    Fluid co2;
    co2.Pc = 73.773e5;
    co2.Tc = 304.13;
    co2.w = 0.22394;

    std::vector<double> DRA_params = {15000., 1., 2.};

    std::vector<double> pressures = {1.e6,
                                     1.e6,
                                     1.e6};

    std::vector<double> ns = pta_model.get_multiple_loadings(pressures, 303., DRA_params, co2);

    for (std::size_t i = 0; i < ns.size(); i++)
    {
        assert(ns[i] == 27.632277562030094);
        std::cout << "Expression valid....\n"
                  << std::endl;
    }

    double abs_deviation = pta_model.get_absolute_deviation(ns, pressures, 303., DRA_params, co2);

    assert(abs_deviation == 0.);
    std::cout << "Expression valid....\n"
              << std::endl;

    return 0;
}