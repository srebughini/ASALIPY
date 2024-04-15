from abc import ABC, abstractmethod

from asali.utils.unit_converter import UnitConverter


class BasicReactorShape(ABC):
    def __init__(self):
        """
        Class that describe reactor shape
        """
        self.uc = UnitConverter()

        self._nusselt = 0.
        self._sherwood = 0.
        self._section_area = 0.
        self._reactor_model = 0.
        self._void_fraction = 0.
        self._tube_diameter = 0.
        self._specific_area = 0.
        self._length = None

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

    @property
    def nusselt(self):
        """
        Nusselt number
        :return: Nusselt number
        """
        return self._nusselt

    @property
    def sherwood(self):
        """
        Sherwood number
        :return: Sherwood number
        """
        return self._sherwood

    @property
    def section_area(self):
        """
        Section area
        :return: Section area
        """
        return self._section_area

    @property
    def reactor_model(self):
        """
        Reactor model
        :return: Reactor model
        """
        return self._reactor_model

    @property
    def void_fraction(self):
        """
        Void fraction
        :return: Void fraction
        """
        return self._void_fraction

    @property
    def tube_diameter(self):
        """
        Tube diameter
        :return: Tube diameter
        """
        return self._tube_diameter

    @property
    def specific_area(self):
        """
        Specific area
        :return: Specific area
        """
        return self._specific_area

    @abstractmethod
    def estimate_mass_transfer_coefficient(self, mass_flow_rate, viscosity, density, diffusivity):
        """
        Estimate mass transfer coefficient based on ReactorModel
        :param mass_flow_rate: Mass flow rate in [kg/s]
        :param viscosity: Gas viscosity in [Pas]
        :param density: Gas density in [kg/m3]
        :param diffusivity: Gas mixture diffusivity in [m2/s]
        :return: Mass transfer coefficient
        """
        pass

    @abstractmethod
    def estimate_heat_transfer_coefficient(self, mass_flow_rate, viscosity, conductivity, specific_heat):
        """
        Estimate heat transfer coefficient based on ReactorModel
        :param mass_flow_rate: Mass flow rate in [kg/s]
        :param viscosity: Gas viscosity in [Pas]
        :param conductivity: Gas thermal conductivity in [W/m/K]
        :param specific_heat: Gas specific heat in [J/kg/K]
        :return: Heat transfer coefficient
        """
        pass
