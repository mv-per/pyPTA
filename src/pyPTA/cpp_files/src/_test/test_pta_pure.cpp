#include <iostream>
#include <vector>
#include <cassert>
#include "../pta_pure.h"

int main(void)
{

    auto pta_model = PurePTA("DRA", "pr77", "excess", 100);

    Fluid co2;
    co2.critical_pressure = 73.773e5;
    co2.critical_temperature = 304.13;
    co2.acentric_factor = 0.22394;

    std::vector<double> DRA_params = {15000., 1., 2.};

    std::vector<double> n2 = pta_model.get_multiple_loadings({1.e6, 1.e6}, 303., DRA_params, co2, 100);

    double n = pta_model.get_loading(1.e6, 303., DRA_params, co2, 100);
    // double n3 = pta_model.get_loading(1.e6, 303., DRA_params, co2, 100);
    assert(n == 27.632277562030094);
    std::cout << "Expression valid....\n";

    assert(n2[1] == 27.632277562030094);
    std::cout << "Expression valid....\n";

    return 0;
}