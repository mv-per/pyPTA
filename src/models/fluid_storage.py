
# db.add(co2)
# db.commit()
from collections import defaultdict

import pandas as pd
from database.db_setup import Session
from database.models import Fluid
from sqlalchemy import or_

# co2 = Fluid()
# co2.name = "Carbon Dioxide"
# co2.formula = "CO2"
# co2.molar_mass = 44.0098               # Molar mass, g / mol
# co2.critical_pressure = 73.773e5              # Critical Pressure, Pa
# co2.critical_temperature = 304.13                # Critical Temperature, K
# co2.critical_compressibility = 0.2746                # Critical compressibility
# co2.accentric_factor = 0.22394                # Accentric factor
# co2.lj_diameter = 3.941           # Lennard-Jonnes diameter



class FluidData:

    def __init__(self):
        self.mixture = None
        self.data = defaultdict(list)
        self.db = Session()

    def getProperties(self, mixture):
        self.mixture = mixture
        self._load_fluid_data_from_db()
        return self.data

    def _store_to_list(self, data):
        for key, value in data.items():
            self.data[key].append(value)

    def _load_fluid_data_from_db(self):
        for component in self.mixture:
            data = self.db.query(Fluid).filter(or_(Fluid.name == component, Fluid.formula == component)).first()
            if data:
                self._store_to_list(data.__dict__)
            else:
                raise Exception(f"Component {component} not found in Database")

    def updateProperty(self, component, data):
        self.db.query(Fluid).filter(or_(Fluid.name == component, Fluid.formula == component)).update(**data)
        self.db.commit()
        self.db.flush()
 
    def loadPropertiesFromCSV(self, csvPath):
       #get all current fluids 
       self.db.query(Fluid).all()

       # load csv


       #loop through CSV data



storage = FluidData()
fluid_properties = storage.getProperties(['CO2', 'Carbon Dioxide'])
# storage.getProperties(['CO2'])

print("a")
