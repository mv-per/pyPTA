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
    owner_fluid_id = Column(Integer, ForeignKey('fluid.id'))
    owner_fluid = relationship("Fluid", back_populates="ideal_cp_parameters")

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
    ideal_cp_parameters = relationship("IdealCpParameters", back_populates="owner_fluid", cascade="all, delete-orphan")


class PureExperimentalData(Base):
    __tablename__ = 'experimental_data'

    id = Column(Integer, primary_key=True)
    _pressure = Column(String, default='')
    _loading = Column(String, default='')

    owner_pure_isotherm_id = Column(Integer, ForeignKey('pure_isotherm.id'))
    owner_pure_isotherm = relationship("PureIsothermModel", back_populates="experimental_data")

    @property
    def pressure(self):
        return [float(x) for x in self._pressure.split(';')]
    @pressure.setter
    def pressure(self, pressure_list):
        self._pressure = ';'.join(map(str, pressure_list))
 
    @property
    def loading(self):
        return [float(x) for x in self._loading.split(';')]
    @loading.setter
    def loading(self, loading_list):   
        self._loading = ';'.join(map(str, loading_list))

class PureIsothermModel(Base):
    __tablename__ = 'pure_isotherm'
    
    id = Column(Integer, primary_key=True)
    name = Column(String)
    fluid_name = Column(String)
    experimental_data = relationship("PureExperimentalData", back_populates="owner_pure_isotherm", cascade="all, delete-orphan")



# Create tables if doesn't exists
Base.metadata.create_all(engine)
