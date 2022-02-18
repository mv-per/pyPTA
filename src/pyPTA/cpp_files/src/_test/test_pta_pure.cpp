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

    double n = pta_model.get_loading(1.e6, 303., DRA_params, co2);
    assert(n == 27.632277562030094);
    cout << "Expression valid....\n";

    std::vector<double> n2 = pta_model.get_multiple_loadings({1.e6, 1.e6}, 303., DRA_params, co2);
    assert(n2[1] == 27.632277562030094);
    cout << "Expression valid....\n";

    return 0;
}