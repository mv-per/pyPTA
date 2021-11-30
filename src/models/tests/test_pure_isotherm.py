

# import pytest

from models.pure_isotherm import PureIsotherm
from database.models import PureIsothermModel

def test_loading_txt(datadir):
    path = str(datadir / 'test_data.txt')

    test_isotherm = PureIsotherm('test_isotherm', 'CO2')
    test_isotherm.load_data_from_txt(path, 0, 2)

    loaded_isotherm = test_isotherm.get_isotherm()
    assert len(loaded_isotherm.experimental_data.loading) > 1
    assert 'test_isotherm' == loaded_isotherm.name

def test_loading_csv(datadir):
    path = str(datadir / 'test_data.csv')

    test_isotherm = PureIsotherm('test_isotherm', 'CO2')
    test_isotherm.load_data_from_csv(path, 0, 2)

    loaded_isotherm = test_isotherm.get_isotherm()
    assert len(loaded_isotherm.experimental_data.loading) > 1
    assert 'test_isotherm' == loaded_isotherm.name


def test_loading_xlsx(datadir):
    path = str(datadir / 'test_data.xlsx')

    test_isotherm = PureIsotherm('test_isotherm', 'CO2')
    test_isotherm.load_data_from_excel(path, 0, 2)

    loaded_isotherm = test_isotherm.get_isotherm()
    assert len(loaded_isotherm.experimental_data.loading) > 1
    assert 'test_isotherm' == loaded_isotherm.name
