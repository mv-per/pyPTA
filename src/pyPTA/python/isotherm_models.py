from typing import List

import attr
from attr import validators
from attrs_strict import type_validator

from pyPTA.fluid_models import Fluid

PRESSURE_UNIT_DICT = {
    "Pa": 1.0,
    "bar": 1.0e5,
    "MPa": 1.0e6,
    "kPa": 1.0e3,
    "psi": 6895.0,
    "atm": 101330.0,
}

LOADING_UNIT_DICT = {
    "mmol/g": 1.0,
    "mol/kg": 1.0,
    "mol/g": 1.0e3,
}


@attr.s(auto_attribs=True)
class Isotherm:
    """
    A Monocomponent Isotherm

    :param: fluid
    :type: :class:`Fluid`
    :param: A list containing N float elements of pressure
    :type: list[float]
    :param: A list containing N float elements of loading amount
    :type: list[float]
    """

    fluid: Fluid = attr.ib(validator=validators.instance_of(Fluid))
    pressure: List[float] = attr.ib(validator=type_validator())
    loading: List[float] = attr.ib(validator=type_validator())

    def __attrs_post_init__(self):
        """Check if pressure and loading have the same length"""
        if not len(self.pressure) == len(self.loading):
            raise Exception(
                f"""DIFFERENT LENGHT: pressure (len={len(self.pressure)})
                 and loading (len={len(self.loading)}) must have the same
                 length of data"""
            )
