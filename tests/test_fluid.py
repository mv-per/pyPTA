# import copy

# from pyPTA.fluid_models import Fluid
# from pyPTA.fluid_models import FluidData


# def test_fluid():
#     """Tests the fluid class"""

#     co2_dict = {
#         "name": "Carbon Dioxide",
#         "critical_pressure": 73.773e5,  # Critical Pressure, Pa
#         "critical_temperature": 304.13,  # Critical Temperature, K
#         "accentric_factor": 0.22394,  # Accentric factor
#     }

#     co2 = Fluid(**co2_dict)

#     co2_2 = copy.copy(co2)
#     assert co2 == co2_2
#     assert co2 is not co2_2


# def test_fluid_data():
#     """Tests for the FluidData"""

#     co2_dict = {
#         "name": "Carbon Dioxide",
#         "critical_pressure": 73.773e5,  # Critical Pressure, Pa
#         "critical_temperature": 304.13,  # Critical Temperature, K
#         "accentric_factor": 0.22394,  # Accentric factor
#     }

#     ch4_dict = {
#         "name": "Methane",
#         "critical_pressure": 73.773e5,  # Critical Pressure, Pa
#         "critical_temperature": 304.13,  # Critical Temperature, K
#         "accentric_factor": 0.22394,  # Accentric factor
#     }

#     co2 = Fluid(**co2_dict)
#     ch4 = Fluid(**ch4_dict)

#     fluid_data = FluidData()
#     fluid_data.add_fluid(co2)
#     fluid_data.add_fluid(ch4)
#     assert len(fluid_data.fluids) == 2

#     fluid_data_2 = FluidData()
#     fluid_data_2.add_fluids(co2, ch4)
#     assert len(fluid_data.fluids) == 2

#     assert fluid_data.fluids == fluid_data_2.fluids
