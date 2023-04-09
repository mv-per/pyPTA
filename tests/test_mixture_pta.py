from pyPTA import MixturePTA, Fluid, DRA_POTENTIAL, LEE_POTENTIAL, STEELE_POTENTIAL, Adsorbent
import pytest
from typing import List

@pytest.fixture
def setup_fluids() -> List[Fluid]:
    fluids =  [Fluid('CO2', 73.773e5, 304.13, 0.22394), Fluid('CH4', 45.99e5, 190.56, 0.01142), Fluid('N2', 3395800.0, 126.192, 0.0372)]

    fluids[0].lennard_jonnes_diameter = 3.941
    fluids[1].lennard_jonnes_diameter = 3.798
    fluids[2].lennard_jonnes_diameter = 3.758

    return fluids

@pytest.mark.parametrize(
    'potential, CO2_params, CH4_params, expected',
    [
        (DRA_POTENTIAL, [7880.19, 0.29, 2.], [5600, 0.36, 3.], [2.50608, 0.87942]),
        (STEELE_POTENTIAL, [109.32, 13.34, 611.8], [109.32, 13.34, 611.8], [6.56534, 0.83682]),
        (LEE_POTENTIAL, [125.63, 12.26, 765.70], [112.63, 10.26, 720.70], [1.74713, 1.35373]),
    ]
)
def test_get_mixture_loading_two_components(setup_fluids:List[Fluid], potential:str, CO2_params:List[float], CH4_params:List[float], expected:List[float])->None:

    bulk_composition = [.25, .75]
    mixture_params = [CO2_params, CH4_params]
    fluids = setup_fluids[:2]
    adsorbent = Adsorbent("Z01x", 3.35, 0.382)

    model = MixturePTA(potential, 'pr77', 'excess', 555)
    model.set_adsorbent(adsorbent)

    calculated = model.get_loading(bulk_composition, 1e6, 310.2, mixture_params, fluids)

    calculated = [round(val, 5) for val in calculated]

    assert len(calculated) == len(bulk_composition)
    assert calculated == expected


@pytest.mark.parametrize(
    'potential, CO2_params, CH4_params, expected',
    [
        (DRA_POTENTIAL, [7880.19, 0.29, 2.], [5600, 0.36, 3.], [2.90177, 1.45198]),
        (STEELE_POTENTIAL, [109.32, 18, 611.8], [92.32, 16, 456.8], [4.80308, 2.51869]),
        (LEE_POTENTIAL, [125.63, 12.26, 611.70], [92.32, 13.34, 456.8], [2.69445, 0.77108]),
    ]
)
def test_get_mixture_loading_two_components_srk(setup_fluids:List[Fluid], potential:str, CO2_params:List[float], CH4_params:List[float], expected:List[float])->None:

    bulk_composition = [0.5, 0.5]
    mixture_params = [CO2_params, CH4_params]
    fluids = setup_fluids[:2]
    adsorbent = Adsorbent("Z01x", 3.35, 0.382)

    model = MixturePTA(potential, 'srk', 'excess', 555)
    model.set_adsorbent(adsorbent)

    calculated = model.get_loading(bulk_composition, 1e6, 310.2, mixture_params, fluids)

    calculated = [round(val, 5) for val in calculated]

    assert len(calculated) == len(bulk_composition)
    assert calculated == expected


@pytest.mark.parametrize(
    'potential, CO2_params, CH4_params, expected',
    [
        (STEELE_POTENTIAL, [125.63, 12.26, 765.70], [112.63, 12.26, 720.70], [3.85778, 0.46537]),
        (LEE_POTENTIAL, [125.63, 12.26, 765.70], [112.63, 10.26, 720.70], [3.85778, 0.46537]),
    ]
)
def test_raise_error_get_mixture_loading_Lj_without_adsorbent(setup_fluids:List[Fluid], potential:str, CO2_params:List[float], CH4_params:List[float], expected:List[float])->None:

    bulk_composition = [0.5, 0.5]
    mixture_params = [CO2_params, CH4_params]
    fluids = setup_fluids[:2]

    model = MixturePTA(potential, 'pr77', 'excess', 555)

    expected_error_message = """Adsorbent properties are needed for LJ-based potentials and is not defined."""
    with pytest.raises(ValueError) as e:
        model.get_loading(bulk_composition, 1e6, 310.2, mixture_params, fluids)

    assert expected_error_message in str(e)

    

@pytest.mark.parametrize(
    'potential, CO2_params, CH4_params, N2_params, expected',
    [(DRA_POTENTIAL, [7880.19, 0.29, 2.], [5600, 0.36, 3.], [5000.5, 0.57, 2.] , [2.52491, 0.28772, 0.38486])]
)
def test_get_mixture_loading_three_components_pr(setup_fluids:List[Fluid], potential:str, CO2_params:List[float], CH4_params:List[float], N2_params:List[float], expected:List[float])->None:

    bulk_composition = [0.25, 0.25, 0.5]
    mixture_params = [CO2_params, CH4_params, N2_params]
    fluids = setup_fluids[:]

    model = MixturePTA(potential, 'pr77', 'excess', 555)

    calculated = model.get_loading(bulk_composition, 1e6, 310.2, mixture_params, fluids)

    calculated = [round(val, 5) for val in calculated]

    assert len(calculated) == len(bulk_composition)
    assert calculated == expected
    
@pytest.mark.parametrize(
    'potential, CO2_params, CH4_params, N2_params, expected',
    [(DRA_POTENTIAL, [7880.19, 0.29, 2.], [5600, 0.36, 3.], [5000.5, 0.57, 2.] , [1.2, 0.3971, 1.18704])]
)
def test_get_mixture_loading_three_components_srk(setup_fluids:List[Fluid], potential:str, CO2_params:List[float], CH4_params:List[float], N2_params:List[float], expected:List[float])->None:

    bulk_composition = [0.25, 0.25, 0.5]
    mixture_params = [CO2_params, CH4_params, N2_params]
    fluids = setup_fluids[:]

    model = MixturePTA(potential, 'srk', 'excess', 555)

    calculated = model.get_loading(bulk_composition, 1e6, 310.2, mixture_params, fluids)

    calculated = [round(val, 5) for val in calculated]

    assert len(calculated) == len(bulk_composition)
    assert calculated == expected
    