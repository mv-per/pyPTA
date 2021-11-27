from typing import List

import numpy as np

PRESSURE_UNIT_DICT = {
    'Pa':1.0,
    'bar':1.0e5,
    'MPa':1.0e6,
    'kPa':1.0e3,
    'psi':6895.0,
    'atm':101330.0
}

LOADING_UNIT_DICT = {
    'mmol/g':1.0,
    'mol/kg':1.0,
    'mol/g': 1.0e3,
}

class MonocomponentDataSerializer:

    def __init__(self, isotherm_name:str, pressure_unit="Pa", loading_unit="mmol/g"):
        """[summary]

        Args:
            isotherm_name (str): [description]
        """
        self.isotherm_name = isotherm_name
        self.loading_unit = loading_unit
        self.pressure_unit = pressure_unit
        self.pressure = []
        self.loading = []

    def get_isotherm(self):
        """ get isotherm data

        Returns:
            list[pressure, loading]
        """
        return [self.pressure, self.loading]

    def load_data_from_txt(self, txtFilePath:str, pressure_col, loading_col, row_index=0, pressure_unit=None, loading_unit=None):
        r""" Load monocomponent data from TXT file

        Args:
            txtFilePath (str): Path to TXT file
            pressure_col (int): Column with pressure data (starts at 0)
            loading_col (int): Column with loading data (starts at 0)
            row_index (int, optional): Row which data starts, Defaults to 0
            pressure_unit (str, optional): Unit of the pressure, Defaults to None
            loading_unit (str, optional): Unit of the loading, Defaults to None
        """
        self.pressure = np.loadtxt(txtFilePath)[row_index:, pressure_col]
        self.loading = np.loadtxt(txtFilePath)[row_index:, loading_col]

        self._fix_units(pressure_unit, loading_unit)
        self._check_data_consistency()

    def load_data_from_csv(self, csvFilePath:str, pressure_unit=None, loading_unit=None):
        """ Load monocomponent data from CSV file

        Args:
            csvFilePath (str): [description]
        """

        self._fix_units(pressure_unit, loading_unit)
        self._check_data_consistency()
        pass

    def load_data_from_excel(self, excelFilePath:str, pressure_unit=None, loading_unit=None, **kwargs):
        """[summary]

        Args:
            excelFilePath (str): Path to excel file
            kwargs:
                pressure_range = [ini_row, ini_col, last_row, last_col]
                loading_range = [ini_row, ini_col, last_row, last_col]
        """
        if pressure_unit or loading_unit:
            self._fix_units(pressure_unit, loading_unit)
        self._check_data_consistency()
        pass

    def _update_loading_unit(self, loading_unit):
        try:
            self.loading = [loading*PRESSURE_UNIT_DICT[loading_unit] for loading in self.loading]       
        except:
            raise Exception(f"Loading unit not found\nAvailable units: {LOADING_UNIT_DICT.items()}")

    def _update_pressure_unit(self, pressure_unit):
        try:
            self.pressure = [pressure*PRESSURE_UNIT_DICT[pressure_unit] for pressure in self.pressure]       
        except:
            raise Exception(f"Pressure unit not found\nAvailable units: {PRESSURE_UNIT_DICT.items()}")

    def _fix_units(self, pressure_unit, loading_unit):
        if isinstance(pressure_unit, str):
            self._update_pressure_unit(pressure_unit)
        if isinstance(loading_unit, str):
            self._update_loading_unit(loading_unit)

    def _check_data_consistency(self):
        if len(self.pressure) != len(self.loading):
            raise Exception("Pass list of pure component data with same number of data")
        
