# import pytest
# from attrs_strict._error import AttributeTypeError

# from pyPTA.fluid_models import Fluid
# from pyPTA.isotherm_models import Isotherm


# def test_cpta():
#     a = 1

#     assert a == 2


# def test_pure_isotherm_success():
#     """Test Success of adding isotherm"""

#     co2_dict = {
#         "name": "Carbon Dioxide",
#         "critical_pressure": 73.773e5,  # Critical Pressure, Pa
#         "critical_temperature": 304.13,  # Critical Temperature, K
#         "accentric_factor": 0.22394,  # Accentric factor
#     }

#     co2 = Fluid(**co2_dict)
#     pressure = [1.0, 2.0, 3.0]
#     loading = [1.0, 2.0, 3.0]

#     co2_isotherm = Isotherm(co2, pressure, loading)

#     assert co2_isotherm


# def test_pure_isotherm_fail_type():
#     """Load list of integers instead of list o floats"""

#     co2_dict = {
#         "name": "Carbon Dioxide",
#         "critical_pressure": 73.773e5,  # Critical Pressure, Pa
#         "critical_temperature": 304.13,  # Critical Temperature, K
#         "accentric_factor": 0.22394,  # Accentric factor
#     }

#     co2 = Fluid(**co2_dict)
#     pressure = [1, 2, 3]
#     loading = [1, 2, 3]

#     with pytest.raises(AttributeTypeError):
#         Isotherm(co2, pressure, loading)


# def test_pure_isotherm_fail_len_data():
#     """Fails to load data with same size"""

#     co2_dict = {
#         "name": "Carbon Dioxide",
#         "critical_pressure": 73.773e5,  # Critical Pressure, Pa
#         "critical_temperature": 304.13,  # Critical Temperature, K
#         "accentric_factor": 0.22394,  # Accentric factor
#     }

#     co2 = Fluid(**co2_dict)
#     pressure = [1.0, 2.0]
#     loading = [1.0, 2.0, 3.0]

#     with pytest.raises(Exception, match=r"DIFFERENT LENGHT: .*"):
#         Isotherm(co2, pressure, loading)


# def test_pure_isotherm_fail_fluid():
#     """Fails to load proper fluid"""

#     co2_dict = {
#         "name": "Carbon Dioxide",
#         "critical_pressure": 73.773e5,  # Critical Pressure, Pa
#         "critical_temperature": 304.13,  # Critical Temperature, K
#         "accentric_factor": 0.22394,  # Accentric factor
#     }

#     pressure = [1.0, 2.0, 3.0]
#     loading = [1.0, 2.0, 3.0]

#     with pytest.raises(TypeError):
#         Isotherm(co2_dict, pressure, loading)


# # def test_loading_txt(datadir):
# #     path = str(datadir / 'test_data.txt')

# #     test_isotherm = PureIsotherm('test_isotherm', 'CO2')
# #     test_isotherm.load_data_from_txt(path, 0, 2)

# #     loaded_isotherm = test_isotherm.get_isotherm()
# #     assert len(loaded_isotherm.experimental_data.loading) > 1
# #     assert 'test_isotherm' == loaded_isotherm.name

# # def test_loading_csv(datadir):
# #     path = str(datadir / 'test_data.csv')

# #     test_isotherm = PureIsotherm('test_isotherm', 'CO2')
# #     test_isotherm.load_data_from_csv(path, 0, 2)

# #     loaded_isotherm = test_isotherm.get_isotherm()
# #     assert len(loaded_isotherm.experimental_data.loading) > 1
# #     assert 'test_isotherm' == loaded_isotherm.name


# # def test_loading_xlsx(datadir):
# #     path = str(datadir / 'test_data.xlsx')

# #     test_isotherm = PureIsotherm('test_isotherm', 'CO2')
# #     test_isotherm.load_data_from_excel(path, 0, 2)

# #     loaded_isotherm = test_isotherm.get_isotherm()
# #     assert len(loaded_isotherm.experimental_data.loading) > 1
# #     assert 'test_isotherm' == loaded_isotherm.name


# # class PureIsotherm:
# #     def __init__(
# #         self, isotherm_name: str, fluid: str, pressure_unit="Pa", loading_unit="mmol/g"
# #     ):
# #         """[summary]

# #         Args:
# #             isotherm_name (str): [description]
# #         """
# #         self.IsothermModel = PureIsothermModel()
# #         self.IsothermModel.name = isotherm_name
# #         self.IsothermModel.fluid_name = fluid
# #         self.loading_unit = loading_unit
# #         self.pressure_unit = pressure_unit
# #         self.pressure = []
# #         self.loading = []

# #     def get_isotherm(self):
# #         """get isotherm

# #         Returns:
# #             PureIsothermModel
# #         """
# #         return self.IsothermModel

# #     def load_data_from_txt(
# #         self,
# #         txtFilePath: str,
# #         pressure_col,
# #         loading_col,
# #         pressure_unit=None,
# #         loading_unit=None,
# #     ):
# #         r"""Load monocomponent data from TXT file

# #         Args:
# #             txtFilePath (str): Path to TXT file
# #             pressure_col (int): Column with pressure data (starts at 0)
# #             loading_col (int): Column with loading data (starts at 0)
# #             row_index (int, optional): Row which data starts, Defaults to 0
# #             pressure_unit (str, optional): Unit of the pressure, Defaults to None
# #             loading_unit (str, optional): Unit of the loading, Defaults to None
# #         """
# #         data = pd.read_csv(txtFilePath, sep=" ", header=None).to_dict()
# #         self.pressure = data[pressure_col]
# #         self.loading = data[loading_col]

# #         self._fix_units(pressure_unit, loading_unit)
# #         self._check_data_consistency()

# #     def _add_data_to_Model(self):
# #         self.IsothermModel.experimental_data.experimental_data = PureExperimentalData()
# #         self.IsothermModel.experimental_data.pressure = self.pressure
# #         self.IsothermModel.experimental_data.loading = self.loading

# #     def load_data_from_csv(
# #         self,
# #         csvFilePath: str,
# #         pressure_col,
# #         loading_col,
# #         pressure_unit=None,
# #         loading_unit=None,
# #     ):
# #         """Load monocomponent data from CSV file

# #         Args:
# #             csvFilePath (str): [description]
# #         """
# #         data = pd.read_csv(csvFilePath, sep=",", header=None).to_dict()
# #         self.pressure = data[pressure_col]
# #         self.loading = data[loading_col]
# #         self._fix_units(pressure_unit, loading_unit)
# #         self._check_data_consistency()

# #     def load_data_from_excel(
# #         self,
# #         excelFilePath: str,
# #         pressure_col,
# #         loading_col,
# #         pressure_unit=None,
# #         loading_unit=None,
# #         **kwargs,
# #     ):
# #         """[summary]

# #         Args:
# #             excelFilePath (str): Path to excel file
# #             kwargs:
# #                 pressure_range = [ini_row, ini_col, last_row, last_col]
# #                 loading_range = [ini_row, ini_col, last_row, last_col]
# #         """
# #         data = pd.read_excel(excelFilePath, header=None).to_dict()
# #         self.pressure = data[pressure_col]
# #         self.loading = data[loading_col]
# #         self._fix_units(pressure_unit, loading_unit)
# #         self._check_data_consistency()
# #         pass

# #     def _update_loading_unit(self, loading_unit):
# #         try:
# #             self.loading = [
# #                 loading * PRESSURE_UNIT_DICT[loading_unit] for loading in self.loading
# #             ]
# #         except:
# #             raise Exception(
# #                 f"Loading unit not found\nAvailable units: {LOADING_UNIT_DICT.items()}"
# #             )

# #     def _update_pressure_unit(self, pressure_unit):
# #         try:
# #             self.pressure = [
# #                 pressure * PRESSURE_UNIT_DICT[pressure_unit]
# #                 for pressure in self.pressure
# #             ]
# #         except:
# #             raise Exception(
# #                 f"Pressure unit not found\nAvailable units: {PRESSURE_UNIT_DICT.items()}"
# #             )

# #     def _fix_units(self, pressure_unit, loading_unit):
# #         if isinstance(pressure_unit, str):
# #             self._update_pressure_unit(pressure_unit)
# #         if isinstance(loading_unit, str):
# #             self._update_loading_unit(loading_unit)

# #     def _check_data_consistency(self):
# #         if len(self.pressure) != len(self.loading):
# #             raise Exception("Pass list of pure component data with same number of data")
# #         else:
# #             self._add_data_to_Model()
