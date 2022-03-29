from asali.reactors.shapes.basic import BasicReactorShape
from asali.utils.input_parser import ReactorSection, ReactorModel

import numpy as np


class PackedBedReactorShape(BasicReactorShape):
    def __init__(self):
        """
        Class that describe packed bed reactor shape
        """
        super().__init__()
        self._particle_diameter = 0.

    @property
    def particle_diameter(self):
        """
        Particle diameter
        :return: Particle diameter
        """
        return self._particle_diameter

    def set_properties(self, particle_diameter, particle_diameter_ud, tube_diameter, tube_diameter_ud, void_fraction):
        """
        Set packed bed reactor shape properties
        :param particle_diameter: Particle diameter
        :param particle_diameter_ud: Particle diameter unit dimensions
        :param tube_diameter: Tube diameter
        :param tube_diameter_ud: Tube diameter unit dimensions
        :param void_fraction: Void fraction
        :return: List of estimated reactor properties
        """
        self._tube_diameter = self.uc.convert_to_meter(tube_diameter, tube_diameter_ud)
        self._void_fraction = void_fraction
        self._particle_diameter = self.uc.convert_to_meter(particle_diameter, particle_diameter_ud)
        self._specific_area = 6. * (1. - self.void_fraction) / self.particle_diameter
        self._section_area = np.pi * 0.25 * np.square(self.tube_diameter)

        self._reactor_model = ReactorModel.PACKEDBED

        return [self.section_area,
                self.void_fraction,
                self.tube_diameter,
                self.reactor_model,
                self.specific_area,
                self.particle_diameter]

    def estimate_mass_transfer_coefficient(self, mass_flow_rate, viscosity, density, diffusivity):
        """
        Estimate mass transfer coefficient based on ReactorModel
        :param mass_flow_rate: Mass flow rate in [kg/s]
        :param viscosity: Gas viscosity in [Pas]
        :param density: Gas density in [kg/m3]
        :param diffusivity: Gas mixture diffusivity in [m2/s]
        :return: Mass transfer coefficient
        """
        reynolds_real = mass_flow_rate * self.particle_diameter / (viscosity * 0.25 * np.square(self.tube_diameter))
        reynolds = reynolds_real / ((1. - self.void_fraction) * 6)

        laminar_flow = np.logical_and(reynolds < 50, reynolds > 0)
        turbulent_flow = reynolds >= 50

        j_mass = np.zeros_like(reynolds_real)
        j_mass[laminar_flow] = 0.91 / (np.power(reynolds[laminar_flow], 0.51))
        j_mass[turbulent_flow] = 0.61 / (np.power(reynolds[turbulent_flow], 0.41))

        Sc = ((viscosity / density) / diffusivity.T).T

        return ((j_mass * reynolds_real) * (np.power(Sc, 1. / 3.) * diffusivity).T).T / self.particle_diameter

    def estimate_heat_transfer_coefficient(self, mass_flow_rate, viscosity, conductivity, specific_heat):
        """
        Estimate heat transfer coefficient based on ReactorModel
        :param mass_flow_rate: Mass flow rate in [kg/s]
        :param viscosity: Gas viscosity in [Pas]
        :param conductivity: Gas thermal conductivity in [W/m/K]
        :param specific_heat: Gas specific heat in [J/kg/K]
        :return: Heat transfer coefficient
        """
        reynolds_real = mass_flow_rate * self.particle_diameter / (
                viscosity * 0.25 * np.square(self.tube_diameter))
        reynolds = reynolds_real / ((1. - self.void_fraction) * 6)

        laminar_flow = np.logical_and(reynolds < 50, reynolds > 0)
        turbulent_flow = reynolds >= 50

        j_heat = np.zeros_like(reynolds_real)
        j_heat[laminar_flow] = 0.91 / (np.power(reynolds[laminar_flow], 0.51))
        j_heat[turbulent_flow] = 0.61 / (np.power(reynolds[turbulent_flow], 0.41))

        Pr = (specific_heat * viscosity) / conductivity

        return j_heat * reynolds_real * np.power(Pr, 1. / 3.) * conductivity / self.particle_diameter
