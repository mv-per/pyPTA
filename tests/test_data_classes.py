from pyPTA import Fluid, Adsorbent, PurePTA

def test_create_fluid() -> None:

    fluid = Fluid()
    assert fluid is not None

def test_create_fluid_with_kwargs() -> None:

    fluid = Fluid(name='CO2', critical_pressure=100.0, critical_temperature=200, accentric_factor=222)

    assert fluid is not None

def test_create_adsorbent() -> None:

    adsorbent = Adsorbent()

    assert adsorbent is not None

def test_create_adsorbent_with_args() -> None:

    adsorbent = Adsorbent('Zeolite', 0.01, 102.4)

    assert adsorbent is not None
    assert adsorbent.name == 'Zeolite'
    assert adsorbent.diameter == 0.01

def test_create_adsorbent_with_kwargs() -> None:

    adsorbent = Adsorbent(name='Zeolite', diameter=0.01, atomic_density=102.4)

    assert adsorbent is not None
    assert adsorbent.name == 'Zeolite'
    assert adsorbent.diameter == 0.01
