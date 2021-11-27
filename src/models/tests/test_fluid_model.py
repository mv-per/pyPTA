from pypta.models.fluid_model import FluidModel

TEST_CO2 = {
            'name':'CO2',
            'critical_pressure': 20,
            'critical_temperature': 304,
            'critical_compressibility':0.9,
            'accentric_factor':0.222,
        }


def test_fluid_model():

    test_fluid = FluidModel(**TEST_CO2)

    assert TEST_CO2 == test_fluid.fluid_dict
    assert test_fluid.name == 'CO2'





    # def mock_load_db(self):
    #     return 'xyz'

    # mocker.patch(FluidModel._load_db, mock_load_db)
