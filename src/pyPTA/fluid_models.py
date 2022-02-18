from typing import List
from typing import Optional

import attr
from attrs_strict import type_validator


@attr.s(auto_attribs=True)
class IdealCpParameters:
    """
    Ideal Cp Parameters

    :ivar float a
    :ivar float b
    :ivar float c
    :ivar float d
    :ivar float e
    """

    a: float = attr.ib(factory=float, validator=type_validator())
    b: float = attr.ib(factory=float, validator=type_validator())
    c: float = attr.ib(factory=float, validator=type_validator())
    d: float = attr.ib(factory=float, validator=type_validator())
    e: float = attr.ib(factory=float, validator=type_validator())


@attr.s(auto_attribs=True)
class Fluid:
    """
    Base fluid class

    :ivar str name:
        Name of the fluid. e.g: "Carbon Dioxide"
    :ivar float critical_temperature:
        Critical Temperature of the fluid. Unit = Kelvin (K)
    :ivar float critical_pressure:
        Critical Pressure of the fluid. Unit = Pascal (Pa)
    :ivar float accentric_factor:
        Accentric factor of the fluid.
    :ivar Optional[float] lj_diameter:
        Lennard-Jonnes diameter of the fluid. Unit= nanometer (nm)
    :ivar Optional[float] lj_energy:
        Lennard-Jonnes diameter of the fluid. Unit = kJ/mol
    """

    name: str
    critical_temperature: float
    critical_pressure: float
    accentric_factor: float
    ideal_cp_parameters: Optional[IdealCpParameters] = attr.ib(
        factory=IdealCpParameters, validator=type_validator()
    )
    lj_diameter: Optional[float] = attr.ib(factory=float, validator=type_validator())
    lj_energy: Optional[float] = attr.ib(factory=float, validator=type_validator())


@attr.s(auto_attribs=True)
class FluidData:
    """_summary_

    :raises Exception: _description_
    :raises Exception: _description_
    :raises Exception: _description_
    :return: _description_
    :rtype: _type_
    """

    fluids: List[Fluid] = attr.ib(
        init=False, validator=type_validator(), default=attr.Factory(list)
    )

    def add_fluid(self, fluid: Fluid, position: Optional[int] = None):
        """
        Add a fluid to the list

        :param fluid: Fluid instance to load in the fluids list
        :type fluid: Fluid
        :param position: position of the fluid in the list
                 (must be the same as the ExperimentalData),
                 defaults to None
        :type position: Optional[int], optional
        """
        if position is not None:
            self.fluids.insert(position, fluid)
        else:
            self.fluids.append(fluid)

    def add_fluids(self, *args):
        """
        Refreshes a list of fluids and add multiple fluids in the input order

        :param `Fluid` args:
            Fluids to add to list

        """
        self.fluids = []
        for fluid in args:
            self.fluids.append(fluid)

    def __repr__(self):
        return f"{[fluid.name for fluid in self.fluids]}"

    def update_fluid(self, position: int, fluid: Fluid):
        """
        Update a specific fluid

        :param position: Position of the fluid to be updated
        :type position: int
        :param fluid: Values to write in that position
        :type fluid: Fluid
        """
        if self.fluids[position]:
            self.fluids[position] = fluid
        else:
            raise Exception(f"Fluid position {position} not found in list")
