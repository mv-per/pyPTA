"""
`cpta`, or `cPTA` is the builtin package with the C++ methods to run the PTA package
"""

from cPTA import Fluid
from cPTA import PurePTA


def test_cpta_pure_ncalc():

    co2: Fluid = Fluid()
    co2.critical_pressure = 73.773e5
    co2.critical_temperature = 304.13
    co2.acentric_factor = 0.22394
    DRA_params = [15000.0, 1.0, 2.0]

    pta_model: PurePTA = PurePTA("DRA", "pr77", "excess", 100)

    n = pta_model.get_loading(1.0e6, 303.0, DRA_params, co2)

    assert n == 27.632277562030094

    ns = pta_model.get_multiple_loadings([1.0e6], 303.0, DRA_params, co2)

    assert ns == [n]
