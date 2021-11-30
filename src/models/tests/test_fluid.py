
import copy
import os

import pytest
from database.models import Fluid
from models.fluid import FluidData


def test_mocked_db(db_session):
    assert len(db_session.query(Fluid).all()) == 3

def test_fluid_data_get_all(db_session, mocker):

    mocker.patch('models.fluid_storage.FluidData._get_db_session', return_value=db_session)

    fluid_storage = FluidData()

    fluids = fluid_storage.get_all()

    assert len(fluids) == 3

def test_fluid_data_add(empty_db_session, mocker, sample_fluid):
    mocker.patch('models.fluid_storage.FluidData._get_db_session', return_value=empty_db_session)

    fluid_storage = FluidData()

    fluid_storage.add_fluid(sample_fluid)

    fluids = fluid_storage.get_all()

    assert len(fluids) == 1


def test_fluid_data_update_fail(empty_db_session, mocker, sample_fluid):
    mocker.patch('models.fluid_storage.FluidData._get_db_session', return_value=empty_db_session)

    fluid_to_update = 'co2'
    with pytest.raises(Exception, match=f"Fluid {fluid_to_update} not found in database"):
        fluid_storage = FluidData()
        fluid_storage.update_fluid(fluid_to_update, sample_fluid)

def test_fluid_data_update(empty_db_session, mocker, sample_fluid, sample_fluid_2):
    mocker.patch('models.fluid_storage.FluidData._get_db_session', return_value=empty_db_session)

    fluid_storage = FluidData()

    fluid_storage.add_fluid(sample_fluid)

    fluid_storage.update_fluid(sample_fluid.name, sample_fluid_2)

    fluid = fluid_storage.get_fluid(sample_fluid_2.name)

    assert fluid.name == sample_fluid_2.name
    assert fluid.formula == sample_fluid_2.formula

def test_fluid_data_load_fluids_from_csv(empty_db_session, mocker, datadir):
    file_path = os.path.join(datadir, 'fluids.csv')
    mocker.patch('models.fluid_storage.FluidData._get_db_session', return_value=empty_db_session)
    fluid_storage = FluidData()

    fluid_storage.load_fluids_from_csv(file_path)

    fluids = fluid_storage.get_all()

    assert len(fluids) == 5

    fluid_storage.load_fluids_from_csv(file_path)

def test_fluid_data_load_fluids_from_excel(empty_db_session, mocker, datadir):
    file_path = os.path.join(datadir, 'fluidsxlsx.xlsx')
    mocker.patch('models.fluid_storage.FluidData._get_db_session', return_value=empty_db_session)
    fluid_storage = FluidData()

    fluid_storage.load_fluids_from_excel(file_path)

    fluids = fluid_storage.get_all()

    assert len(fluids) == 5

def test_fluid_data_get_fluids(empty_db_session, mocker, datadir):
    file_path = os.path.join(datadir, 'fluidsxlsx.xlsx')
    mocker.patch('models.fluid_storage.FluidData._get_db_session', return_value=empty_db_session)
    fluid_storage = FluidData()

    fluid_storage.load_fluids_from_excel(file_path)

    fluids = fluid_storage.get_fluids(['CO2', 'Methane', 'Nitrogen'])

    assert len(fluids) == 3
    