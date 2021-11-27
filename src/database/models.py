from typing import List

from sqlalchemy import Column, Float, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship

from database.db_setup import engine

Base = declarative_base()


class IdealCpParameters(Base):
    __tablename__ = 'ideal_cp_parameters'

    id = Column(Integer, primary_key=True)
    a = Column(Float, nullable=True)
    b = Column(Float, nullable=True)
    c = Column(Float, nullable=True)
    d = Column(Float, nullable=True)
    e = Column(Float, nullable=True)
    fluid_id = Column(Integer, ForeignKey('fluid.id'))
    fluid = relationship("Fluid", back_populates="ideal_cp_parameters")

class Fluid(Base):
    __tablename__ = 'fluid'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    formula = Column(String, nullable=True)
    critical_temperature = Column(Float)
    critical_pressure = Column(Float)
    critical_volume = Column(Float, nullable=True)
    critical_compressibility = Column(Float, nullable=True)
    accentric_factor = Column(Float)
    lj_diameter = Column(Float, nullable=True)
    lj_energy = Column(Float, nullable=True)
    molar_mass = Column(Float)
    ideal_cp_parameters = relationship("IdealCpParameters", back_populates="fluid", cascade="all, delete-orphan")


class ExperimentalData(Base):
    __tablename__ = 'experimental_data'

    id = Column(Integer, primary_key=True)
    _pressure = Column(String, default='')
    _loading = Column(String, default='')

    mono_isotherm_id = Column(Integer, ForeignKey('mono_isotherm.id'))
    mono_isotherm = relationship("MonoIsotherm", back_populates="experimental_data")

    @property
    def pressure(self):
        return [float(x) for x in self._pressure.split(';')]
    @pressure.setter
    def pressure(self, pressure_list: List[float]):
        self._pressure = [str(x)+';' for x in pressure_list]

    
    @property
    def loading(self):
        return [float(x) for x in self._loading.split(';')]
    @loading.setter
    def loading(self, loading_list: List[float]):
        self._loading = [str(x)+';' for x in loading_list]

class MonoIsotherm(Base):
    __tablename__ = 'mono_isotherm'
    
    id = Column(Integer, primary_key=True)
    name = Column(String)
    fluid_name = Column(String)
    experimental_data = relationship("ExperimentalData", back_populates="mono_isotherm", cascade="all, delete-orphan")



# Create tables if doesn't exists
Base.metadata.create_all(engine)
