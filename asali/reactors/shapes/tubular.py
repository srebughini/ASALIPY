from asali.reactors.shapes.basic import BasicReactorShape
from asali.utils.input_parser import ReactorSection, ReactorModel

import numpy as np


class TubularReactorShape(BasicReactorShape):
    def __init__(self):
        """
        Class that describe tubular reactor shape
        """
        super().__init__()

    def set_properties(self, reactor_section, tube_diameter, tube_diameter_ud, wall_thickness, wall_thickness_ud):
        """
        Set tubular reactor shape properties
        :param reactor_section: Reactor section shape
        :param tube_diameter: Tube diameter
        :param tube_diameter_ud: Tube diameter unit dimensions
        :param wall_thickness: Wall thickness
        :param wall_thickness_ud: Wall thickness unit dimensions
        :return: List of estimated reactor properties
        """
        Dt = self.uc.convert_to_meter(tube_diameter, tube_diameter_ud)
        tw = self.uc.convert_to_meter(wall_thickness, wall_thickness_ud)
        self._tube_diameter = Dt - 2. * tw
        self._specific_area = 4. / self._tube_diameter

        if reactor_section == ReactorSection.CIRCLE:
            self._nusselt = 3.659
            self._sherwood = 3.659
        elif reactor_section == ReactorSection.SQUARE:
            self._nusselt = 2.977
            self._sherwood = 2.977
        elif reactor_section == ReactorSection.TRIANGLE:
            self._nusselt = 2.494
            self._sherwood = 2.494

        self._void_fraction = np.square(self.tube_diameter / Dt)
        self._section_area = np.pi * 0.25 * np.square(self.tube_diameter)
        self._reactor_model = ReactorModel.TUBULAR

        return [self.nusselt,
                self.sherwood,
                self.section_area,
                self.reactor_model,
                self.void_fraction,
                self.tube_diameter,
                self.specific_area]

    def estimate_mass_transfer_coefficient(self, mass_flow_rate, viscosity, density, diffusivity):
        """
        Estimate mass transfer coefficient based on ReactorModel
        :param mass_flow_rate: Mass flow rate in [kg/s]
        :param viscosity: Gas viscosity in [Pas]
        :param density: Gas density in [kg/m3]
        :param diffusivity: Gas mixture diffusivity in [m2/s]
        :return: Mass transfer coefficient
        """
        reynolds_real = mass_flow_rate * self.tube_diameter / (
                viscosity * self.void_fraction * 0.25 * np.square(self.tube_diameter))
        Sc = ((viscosity / density) / diffusivity.T).T
        z_star = np.fabs(np.maximum(1e-06, self._length)) / (self.tube_diameter * reynolds_real)
        z_star = (z_star / Sc.T).T
        return ((self.sherwood + 6.874 * np.power(1000 * z_star, -0.488) * np.exp(
            -57.2 * z_star)) * diffusivity) / self.tube_diameter

    def estimate_heat_transfer_coefficient(self, mass_flow_rate, viscosity, conductivity, specific_heat):
        """
        Estimate heat transfer coefficient based on ReactorModel
        :param mass_flow_rate: Mass flow rate in [kg/s]
        :param viscosity: Gas viscosity in [Pas]
        :param conductivity: Gas thermal conductivity in [W/m/K]
        :param specific_heat: Gas specific heat in [J/kg/K]
        :return: Heat transfer coefficient
        """
        reynolds_real = mass_flow_rate * self.tube_diameter / (
                viscosity * self.void_fraction * 0.25 * np.square(self.tube_diameter))
        Pr = specific_heat * viscosity / conductivity
        z_star = np.fabs(np.maximum(1e-06, self.length)) / (self.tube_diameter * reynolds_real * Pr)
        return ((self.nusselt + 8.827 * np.power(1000 * z_star, -0.545) * np.exp(
            -48.2 * z_star)) * conductivity) / self.tube_diameter
