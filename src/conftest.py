import copy

import pytest
from pytest_mock_resources import Rows, create_sqlite_fixture

from database.models import Base, Fluid

co2 = Fluid()
co2.name = "Carbon Dioxide"
co2.formula = "CO2"
co2.molar_mass = 44.0098               # Molar mass, g / mol
co2.critical_pressure = 73.773e5              # Critical Pressure, Pa
co2.critical_temperature = 304.13                # Critical Temperature, K
co2.critical_compressibility = 0.2746                # Critical compressibility
co2.accentric_factor = 0.22394                # Accentric factor
co2.lj_diameter = 3.941           # Lennard-Jonnes diameter

ch4 = copy.copy(co2)
ch4.name = "Methane"
ch4.formula = "CH4"


n2 = copy.copy(co2)
n2.name = "Nitrogen"
n2.formula = "N2"

rows = Rows(
    co2, ch4, n2
)

@pytest.fixture
def sample_fluid():
    return co2

@pytest.fixture
def sample_fluid_2():
    return ch4


db_session = create_sqlite_fixture(Base, rows,session=True)
empty_db_session = create_sqlite_fixture(Base,session=True)
