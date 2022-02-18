"""
`cpta`, or `cPTA` is the builtin package with the C++ methods to run the PTA package
"""

import cpta


def test_cpta_pure_ncalc():

    eos_params = [73.773e5, 304.13, 0.22394]
    DRA_params = [15000.0, 1.0, 2.0]

    pta_model = cpta.PurePTA("DRA", "pr77", "excess", 100)

    n = pta_model.get_loading(1.0e6, 303.0, DRA_params, eos_params)

    assert n == 27.632277562030094

    ns = pta_model.get_multiple_loadings([1.0e6], 303.0, DRA_params, eos_params)

    assert ns == [n]
