

# import pytest

from pypta.models.mono_component import MonocomponentBaseModel


def test_loading_txt(datadir):
    path = str(datadir / 'test_data.txt')

    test_isotherm = MonocomponentBaseModel('test_isotherm')
    test_isotherm.load_data_from_txt(path, 0, 2)

    loaded_isotherm = test_isotherm.get_isotherm()
    assert len(loaded_isotherm) == 2
    assert 'test_isotherm' == test_isotherm.isotherm_name


