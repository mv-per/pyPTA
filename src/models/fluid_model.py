

class FluidModel:

    def __init__(self, **kwargs):
        self.fluid_dict = {
            'name':None,
            'critical_pressure':None,
            'critical_temperature':None,
            'critical_compressibility':None,
            'accentric_factor':None,
        }
        self.load_data(**kwargs)
    
    def load_data(self, **kwargs):
        for key,value in kwargs.items():
            try:
                self.fluid_dict[key] = value
                setattr(self, key, value)
            except:
                raise Exception(f"{key} not found in fluid data")

    def get_dict_data(self):
        return self.fluid_dict
        