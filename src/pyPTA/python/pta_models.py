from typing import List
from typing import Union

import attr
from attrs_strict import type_validator

from pyPTA.isotherm_models import Isotherm

POTENTIALS = {
    "DRA": [15000.0, 0.7, 2.0],
    # "STEELE": [15000., 0.7, 2.],
    # "LEE": [15000., 0.7, 2.],
}


attr.s(auto_attribs=True, frozen=True)


class MonocomponentModel:
    """_summary_

    :raises AttributeError: _description_
    """

    energy_potential: str = attr.ib(
        validator=type_validator(), on_setattr=attr.setters.validate
    )
    isotherms: List[Isotherm] = attr.ib(validator=type_validator())

    energy_parameters: List[float] = attr.ib(init=False)
    distance_parameters: Union[List[float], float] = attr.ib(init=False)
    heterogeneity_parameters: Union[List[float], float] = attr.ib(init=False)
    initial_conditions: List[float] = attr.ib(init=False)
    single_distance_parameter: bool = attr.ib(init=False, default=True)
    single_heterogeneity_parameter: bool = attr.ib(init=False, default=True)

    def __attrs_post_init__(self):
        #
        num_of_components = len(self.isotherms)
        # Validate model
        if self.energy_potential not in POTENTIALS:
            raise AttributeError(
                f"{self.energy_potential} not found. Available potentials: {POTENTIALS.keys()}"
            )

        self.energy_parameters = [
            POTENTIALS[self.energy_potential][0]
        ] * num_of_components

        if self.single_distance_parameter:
            self.distance_parameters = [POTENTIALS[self.energy_potential][1]] * 1
        else:
            self.distance_parameters = [
                POTENTIALS[self.energy_potential][1]
            ] * num_of_components

        if self.single_heterogeneity_parameter:
            self.heterogeneity_parameters = [POTENTIALS[self.energy_potential][2]] * 1
        else:
            self.heterogeneity_parameters = [
                POTENTIALS[self.energy_potential][2]
            ] * num_of_components
