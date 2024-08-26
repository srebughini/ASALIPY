from abc import ABC
from types import SimpleNamespace


class BasicSupporter(ABC):
    def __init__(self, cls):
        """
        Abstract class to create support classes (e.g. plotters, savers)
        :param cls: Object of class to be plotted
        """
        self.cls = cls

        self._mass_fraction = None
        self._mole_fraction = None
        self._coverage = None
        self._temperature = None

        self._mass_fraction_wall = None
        self._mole_fraction_wall = None
        self._temperature_wall = None

        self._x = None
        self._length = None

        self._variables_name = SimpleNamespace(time="Time (s)",
                                               length="Length (m)",
                                               gas_mole_fraction="Gas mole fraction (-)",
                                               solid_mole_fraction="Solid mole fraction (-)",
                                               gas_mass_fraction="Gas mass fraction (-)",
                                               solid_mass_fraction="Solid mass fraction (-)",
                                               coverage="Site fraction (-)",
                                               gas_temperature="Gas temperature (K)",
                                               solid_temperature="Solid temperature (K)")

    def get_mass_fraction(self):
        """
        Mass fraction getter
        :return: Mass fraction
        """
        return self._mass_fraction

    def set_mass_fraction(self, value):
        """
        Mass fraction setter
        :param value: Mass fraction
        :return:
        """
        self._mass_fraction = value

    # Creating a property object for Mass fraction
    mass_fraction = property(get_mass_fraction, set_mass_fraction)

    def get_mass_fraction_wall(self):
        """
        Solid mass fraction getter
        :return: Mass fraction
        """
        return self._mass_fraction_wall

    def set_mass_fraction_wall(self, value):
        """
        Solid mass fraction setter
        :param value: Mass fraction
        :return:
        """
        self._mass_fraction_wall = value

    # Creating a property object for Solid mass fraction
    mass_fraction_wall = property(get_mass_fraction_wall, set_mass_fraction_wall)

    def get_mole_fraction(self):
        """
        Mole fraction getter
        :return: Mole fraction
        """
        return self._mole_fraction

    def set_mole_fraction(self, value):
        """
        Mole fraction setter
        :param value: Mole fraction
        :return:
        """
        self._mole_fraction = value

    # Creating a property object for Mole fraction
    mole_fraction = property(get_mole_fraction, set_mole_fraction)

    def get_mole_fraction_wall(self):
        """
        Solid mole fraction getter
        :return: Mole fraction
        """
        return self._mole_fraction_wall

    def set_mole_fraction_wall(self, value):
        """
        Solid mole fraction setter
        :param value: Mole fraction
        :return:
        """
        self._mole_fraction_wall = value

    # Creating a property object for Solid mole fraction
    mole_fraction_wall = property(get_mole_fraction_wall, set_mole_fraction_wall)

    def get_coverage(self):
        """
        Site fraction getter
        :return: Site fraction
        """
        return self._coverage

    def set_coverage(self, value):
        """
        Site fraction setter
        :param value: Site fraction
        :return:
        """
        self._coverage = value

    # Creating a property object for Site fraction
    coverage = property(get_coverage, set_coverage)

    def get_temperature(self):
        """
        Temperature getter
        :return: Temperature
        """
        return self._temperature

    def set_temperature(self, value):
        """
        Temperature setter
        :param value: Temperature
        :return:
        """
        self._temperature = value

    # Creating a property object for Temperature
    temperature = property(get_temperature, set_temperature)

    def get_temperature_wall(self):
        """
        Solid temperature getter
        :return: Temperature
        """
        return self._temperature_wall

    def set_temperature_wall(self, value):
        """
        Solid temperature setter
        :param value: Temperature
        :return:
        """
        self._temperature_wall = value

    # Creating a property object for Solid temperature
    temperature_wall = property(get_temperature_wall, set_temperature_wall)

    def get_x(self):
        """
        Independent variable getter
        :return: Independent variable
        """
        return self._x

    def set_x(self, value):
        """
        Independent variable setter
        :param value: Independent variable
        :return:
        """
        self._x = value

    # Creating a property object for Independent variable
    x = property(get_x, set_x)

    def get_length(self):
        """
        Length getter
        :return: Length
        """
        return self._length

    def set_length(self, value):
        """
        Length setter
        :param value: Length
        :return:
        """
        self._length = value

    # Creating a property object for Length
    length = property(get_length, set_length)
