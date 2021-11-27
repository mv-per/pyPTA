


class Unit:
    # TODO
    def __init__(self, value, unit, description=None):

        self.value = value
        self.unit = unit

        if description:
            self.description = description
        else:
            self.description = '-'

    