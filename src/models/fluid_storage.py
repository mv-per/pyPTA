
# db.add(co2)
# db.commit()
from collections import defaultdict

import pandas as pd
from database.db_setup import Session
from database.models import Fluid
from sqlalchemy import or_


class FluidData:

    def __init__(self):
        self.mixture = None
        self.fluids = []
        self.db_session = self._get_db_session()

    def _get_db_session(self):
        return Session()

    def get_all(self):
        return self.db_session.query(Fluid).all()

    def getProperties(self, mixture):
        self.mixture = mixture
        self._load_fluid_data_from_db()
        return self.data

    def _load_fluid_data_from_db(self):
        for fluid in self.mixture:
            data = self.db_session.query(Fluid).filter(or_(Fluid.name == fluid, Fluid.formula == fluid)).first()
            if data:
                self.fluids.append(data)
            else:
                raise Exception(f"Fluid {fluid} not found in Database")

    def get_fluid(self, fluid):
        _fluid = self.db_session.query(Fluid).filter(or_(Fluid.name == fluid, Fluid.formula == fluid)).first()
        if _fluid:
            return _fluid
        else:
            raise Exception(f"Fluid {fluid} not found in database")

    def update_fluid(self, fluid, data):
        _fluid = self.db_session.query(Fluid).filter(or_(Fluid.name == fluid, Fluid.formula == fluid)).first()
        if _fluid:
            for key, value in data.__dict__.items():
                try:
                    setattr(_fluid, key, value)
                except:
                    pass
            self.db_session.commit()
            self.db_session.flush()
        else:
            raise Exception(f"Fluid {fluid} not found in database")

    def add_fluid(self, data:Fluid):
        self.db_session.add(data)
        self.db_session.commit()

 
    def load_fluids_from_dict(self, data_dict):
        for i in range(len(data_dict['name'])):
            _fluid_name =data_dict['name'][i]
            _fluid_formula =data_dict['formula'][i]
            exists = self.db_session.query(Fluid).filter(or_(Fluid.name == _fluid_name, Fluid.formula == _fluid_formula)).first()

            fluid = Fluid()
            for key, values in data_dict.items():
                setattr(fluid, key, values[i])

            if exists:
                print(f"Fluid {_fluid_name} already on database, updating...")
                self.update_fluid(_fluid_name, fluid)
            else:
                self.add_fluid(fluid)
            del fluid


    def load_fluids_from_csv(self, csvPath):
        # load csv to dict
        data_dict = pd.read_csv(csvPath).to_dict()
        
        # store fluids
        self.load_fluids_from_dict(data_dict)

    
    def load_fluids_from_excel(self, excelPath):
        # load excel to dict
        data_dict = pd.read_excel(excelPath, engine='openpyxl').to_dict()
        
        # store fluids
        self.load_fluids_from_dict(data_dict)



       #loop through CSV data



