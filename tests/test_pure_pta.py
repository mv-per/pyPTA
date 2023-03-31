from pyPTA import PurePTA, Fluid, Adsorbent, DRA_POTENTIAL, STEELE_POTENTIAL, LEE_POTENTIAL
import pytest

@pytest.fixture
def setup_fluid() -> Fluid:
    return Fluid('CO2', 73.773e5, 304.13, 0.22394)

@pytest.mark.parametrize(
    'adsorption_potential', (DRA_POTENTIAL, LEE_POTENTIAL, STEELE_POTENTIAL)
)
def test_create_pure_pta(adsorption_potential:str)->None:

    solver = PurePTA(adsorption_potential, 'pr77', 'excess', 555)

    assert solver is not None

@pytest.mark.parametrize(
    'adsorption_potential, parameters, expected_calculated_loading',
    [
        (DRA_POTENTIAL, [7880.19, 0.29, 2.], 5.7589324),
        (LEE_POTENTIAL, [125.63, 12.26, 765.70], 0.2207866),
        (STEELE_POTENTIAL, [109.32, 13.34, 611.88], 0.06526460),
    ]
)
def test_pure_pta_calculate_loading(setup_fluid:Fluid, adsorption_potential:str, parameters:list[float], expected_calculated_loading:float)->None:

    pure_pta = PurePTA(adsorption_potential, 'pr77', 'excess', 555)

    adsorbent = Adsorbent("Z01x", 3.35, 0.382)
    pure_pta.set_adsorbent(adsorbent)

    loading = pure_pta.get_loading(1e6, 305, parameters, setup_fluid)

    assert loading == pytest.approx(expected_calculated_loading)

@pytest.mark.parametrize(
    'adsorption_potential, parameters',
    [
        (LEE_POTENTIAL, [125.63, 12.26, 765.70]),
        (STEELE_POTENTIAL, [109.32, 13.34, 611.88]),
    ]
)
def test_pure_pta_calculate_loading_error_without_adsorbent(setup_fluid:Fluid, adsorption_potential:str, parameters:list[float])->None:

    pure_pta = PurePTA(adsorption_potential, 'pr77', 'excess', 555)

    adsorbent = Adsorbent("Z01x", 3.35, 0.382)

    with pytest.raises(Exception) as exception:
        pure_pta.get_loading(1e6, 305, parameters, setup_fluid)

    print(exception)

def test_pure_pta_calculate_loadings(setup_fluid:Fluid)->None:
    import numpy as np

    pure_pta = PurePTA(DRA_POTENTIAL, 'pr77', 'excess', 555)

    DRA_PARAMETERS = [7880.19, 0.29, 2.]

    pressures = np.arange(1e6, 11e6, 1e6)

    loadings = pure_pta.get_loadings(pressures, 305, DRA_PARAMETERS, setup_fluid)

    assert len(loadings) == len(pressures)

    assert loadings == [5.758932400768459,
            6.895706537068475,
            7.304286537152829,
            7.42048187387874,
            7.368745932028131,
            7.173182017172223,
            6.720995771700607,
            4.603105939555204,
            4.11310650481496,
            3.8378197855658422]

@pytest.mark.parametrize(
    'type_of_deviation, expected',
    [
        ('absolute_relative', 0.0422067131765064),
        ('absolute', 0.023953988843656404),
        ('relative', -0.011296498267532325),
    ]
)
def test_pure_pta_calculate_deviation_of_range(setup_fluid:Fluid, type_of_deviation:str, expected:float)->None:
    import numpy as np

    pure_pta = PurePTA(DRA_POTENTIAL, 'pr77', 'excess', 555)

    DRA_PARAMETERS = [7880.19, 0.29, 2.]

    pressures = np.arange(1e6, 11e6, 1e6)

    experimental_loadings = [5.758932400768459,
            6.895706537068475,
            7.304286537152829,
            7.42048187387874,
            7.368745932028131,
            7.173182017172223,
            6.720995771700607,
            4.603105939555204,
            4.11310650481496,
            3.8378197855658422]

    experimental_loadings = [round(loading, 2) for loading in experimental_loadings]

    deviation = pure_pta.get_deviation_range(type_of_deviation, experimental_loadings, pressures, 305, DRA_PARAMETERS, setup_fluid)

    assert deviation == expected
