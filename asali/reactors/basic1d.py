from abc import ABC

from asali.reactors.basic import BasicReactor

import numpy as np

from asali.utils.input_parser import InputParser, ReactorType


class Basic1DReactor(BasicReactor, ABC):
    def __init__(self, cantera_input_file, gas_phase_name, surface_phase_name):
        """
        Abstract class representing 1D reactor models
        :param cantera_input_file: Cantera file path
        :param gas_phase_name: Cantera gas phase name
        :param surface_phase_name: Cantera interface phase name
        """
        super().__init__(cantera_input_file=cantera_input_file, gas_phase_name=gas_phase_name,
                         surface_phase_name=surface_phase_name)

        self.is_mass_flow_rate = True
        self.gas_diffusion = False

        self.area = 0.
        self.inlet_mass_flow_rate = 0.
        self.inlet_volumetric_flow_rate = 0.
        self.inlet_temperature = 0.
        self.inert_specie_index = 0
        self.inert_coverage_index = 0

        self.inlet_mass_fraction = None
        self.inlet_mole_fraction = None
        self.length = None
        self.n_p = None
        self.alg = None

    def set_length(self, value, unit_dimension):
        """
        Set length
        :param value: Length value
        :param unit_dimension: Length unit dimension
        :return: Discretize length in [m]
        """
        if isinstance(value, (list, np.ndarray)):
            if value[0] != 0.:
                length = np.zeros([len(value) + 1], dtype=np.float64)
                length[1:] = value
                self.length = self.uc.convert_to_meter(length, unit_dimension)
            else:
                self.length = self.uc.convert_to_meter(value, unit_dimension)
        else:
            length = self.uc.convert_to_meter(value, unit_dimension)
            self.length = np.linspace(0, length, num=10)

        self.n_p = self.length.size
        self._setup.n_p = self.n_p
        self._setup.length = self.length
        return self.length

    def set_diameter(self, value, unit_dimension):
        """
        Set diameter
        :param value: Diameter value
        :param unit_dimension: Diameter unit dimension
        :return: Diameter in [m]
        """
        diameter = self.uc.convert_to_meter(value, unit_dimension)
        self.area = np.pi * 0.25 * np.square(diameter)
        self._setup.area = self.area
        self._setup.diameter = diameter
        return diameter

    def set_gas_diffusion(self, value):
        """
        Enable/Disable gas diffusion
        :param value: Variable to enable/disable gas diffusion
        :return: Bool for gas diffusion
        """
        self.gas_diffusion = InputParser.true_parser(value)
        self._setup.gas_diffusion = self.gas_diffusion

    def set_mass_flow_rate(self, value, unit_dimension):
        """
        Set mass flow rate
        :param value: Mass flow rate value
        :param unit_dimension: Mass flow rate unit dimension
        :return: Mass flow rate in [kg/s]
        """
        self.inlet_mass_flow_rate = self.uc.convert_to_kg_per_seconds(value, unit_dimension)
        self.inlet_volumetric_flow_rate = 0.
        self.is_mass_flow_rate = True
        self._setup.inlet_mass_flow_rate = self.inlet_mass_flow_rate
        self._setup.inlet_volumetric_flow_rate = self.inlet_volumetric_flow_rate
        self._setup.is_mass_flow_rate = self.is_mass_flow_rate
        return self.inlet_mass_flow_rate

    def set_volumetric_flow_rate(self, value, unit_dimension):
        """
        Set volumetric flow rate
        :param value: Volumetric flow rate value
        :param unit_dimension: Volumetric flow rate unit dimension
        :return: Volumetric flow rate in [m3/s]
        """
        self.inlet_volumetric_flow_rate = self.uc.convert_to_cubic_meter_per_seconds(value, unit_dimension)
        self.inlet_mass_flow_rate = 0.
        self.is_mass_flow_rate = False
        self._setup.inlet_mass_flow_rate = self.inlet_mass_flow_rate
        self._setup.inlet_volumetric_flow_rate = self.inlet_volumetric_flow_rate
        self._setup.is_mass_flow_rate = self.is_mass_flow_rate
        return self.inlet_volumetric_flow_rate

    def set_inlet_mass_fraction(self, value):
        """
        Set inlet mass fraction
        :param value: Mass fraction
        :return: Inlet mass fraction
        """
        self.gas.Y = value
        self.inlet_mass_fraction = self.gas.Y
        self.inlet_mole_fraction = self.gas.X
        self._setup.inlet_mass_fraction = self.inlet_mass_fraction
        self._setup.inlet_mole_fraction = self.inlet_mole_fraction
        return self.inlet_mass_fraction

    def set_inlet_mole_fraction(self, value):
        """
        Set inlet mole fraction
        :param value: Mole fraction
        :return: Inlet mole fraction
        """
        self.gas.X = value
        self.inlet_mass_fraction = self.gas.Y
        self.inlet_mole_fraction = self.gas.X
        self._setup.inlet_mass_fraction = self.inlet_mass_fraction
        self._setup.inlet_mole_fraction = self.inlet_mole_fraction
        return self.inlet_mole_fraction

    def set_inlet_temperature(self, value, unit_dimension):
        """
        Set inlet temperature
        :param value: Temperature value
        :param unit_dimension: Temperature unit dimension
        :return: Temperature in [K]
        """
        self.inlet_temperature = self.uc.convert_to_kelvin(value, unit_dimension)
        self._setup.inlet_temperature = self.inlet_temperature
        return self.inlet_temperature

    def set_inert_specie(self, specie_name):
        """
        Set inert specie
        :param specie_name: Specie name
        :return: Specie index
        """
        self.inert_specie_index = self.gas.species_index(specie_name)
        self._setup.inert_specie_index = self.inert_specie_index
        return self.inert_specie_index

    def set_inert_coverage(self, coverage_name):
        """
        Set inert coverage species
        :param coverage_name: Coverage specie name
        :return: Coverage specie index
        """
        self.inert_coverage_index = self.surf.species_index(coverage_name)
        self._setup.inert_coverage_index = self.inert_coverage_index
        return self.inert_coverage_index
