from asali.reactors.basic import BasicReactor

import numpy as np

from asali.reactors.ph1d_steady_state import SteadyStatePseudoHomogeneous1DReactor
from asali.reactors.ph1d_transient import TransientPseudoHomogeneous1DReactor
from asali.utils.input_parser import InputParser, ReactorType, ResolutionMethod


class PseudoHomogeneous1DReactor(BasicReactor):
    def __init__(self, cantera_input_file, gas_phase_name, surface_phase_name):
        """
        Class representing PseudoHomogeneous 1D reactor model
        :param cantera_input_file: Cantera file path
        :param gas_phase_name: Cantera gas phase name
        :param surface_phase_name: Cantera interface phase name
        """
        super().__init__(cantera_input_file=cantera_input_file, gas_phase_name=gas_phase_name,
                         surface_phase_name=surface_phase_name)

        self.solution_parser.reactor_type = ReactorType.PSEUDOHOMOGENEOUSPFR

        self.is_mass_flow_rate = True
        self.gas_diffusion = False

        self.volume = 0.
        self.area = 0.
        self.inlet_mass_flow_rate = 0.
        self.inlet_volumetric_flow_rate = 0.
        self.inlet_temperature = 0.
        self.inert_specie_index = 0
        self.inert_coverage_index = 0

        self.inlet_mass_fraction = None
        self.inlet_mole_fraction = None
        self.length = None

        self.reactor_equations = None

    def set_resolution_method(self, method):
        """
        Set resolution method
        :param method: Resolution method as string
        :return: ResolutionMethod object
        """
        self.solution_parser.resolution_method = InputParser.resolution_parser(method)
        return self.solution_parser.resolution_method

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
        return diameter

    def set_gas_diffusion(self, value):
        """
        Enable/Disable gas diffusion
        :param value: Variable to enable/disable gas diffusion
        :return: Bool for gas diffusion
        """
        self.gas_diffusion = InputParser.true_parser(value)

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
        return self.inlet_mole_fraction

    def set_inlet_temperature(self, value, unit_dimension):
        """
        Set inlet temperature
        :param value: Temperature value
        :param unit_dimension: Temperature unit dimension
        :return: Temperature in [K]
        """
        self.inlet_temperature = self.uc.convert_to_kelvin(value, unit_dimension)
        return self.inlet_temperature

    def set_inert_specie(self, specie_name):
        """
        Set inert specie
        :param specie_name: Specie name
        :return: Specie index
        """
        self.inert_specie_index = self.gas.species_index(specie_name)
        return self.inert_specie_index

    def set_inert_coverage(self, coverage_name):
        """
        Set inert coverage species
        :param coverage_name: Coverage specie name
        :return: Coverage specie index
        """
        self.inert_coverage_index = self.surf.species_index(coverage_name)
        return self.inert_coverage_index

    def initial_condition_steady_state(self):
        """
        Function creating the initial condition of the Steady State solution
        :return: Matrix representing the initial mass fraction, coverage and temperature
        """
        if not self.is_mass_flow_rate:
            self.gas.TPY = self.inlet_temperature, self.pressure, self.inlet_mass_fraction
            self.inlet_mass_flow_rate = self.inlet_volumetric_flow_rate * self.gas.density
        return np.block([self.inlet_mass_fraction, self.initial_coverage, self.inlet_temperature])

    def initial_condition_transient(self):
        """
        Generate initial conditions for TRANSIENT model
        :return: Vector/Matrix representing the initial conditions
        """
        NP = self.length.size
        NV = self.gas.n_species + self.surf.n_species + 1
        y0_matrix = np.zeros([NP, NV], dtype=np.float64)

        y0_matrix[:, :self.gas.n_species] = self.gas.Y
        y0_matrix[:, self.gas.n_species:self.gas.n_species + self.surf.n_species] = self.surf.coverages
        y0_matrix[:, -1] = self.initial_temperature

        if not self.energy:
            self.inlet_temperature = self.initial_temperature

        y0_matrix[0, :self.gas.n_species] = self.inlet_mass_fraction
        y0_matrix[0, -1] = self.inlet_temperature

        if not self.is_mass_flow_rate:
            self.gas.TPY = self.inlet_temperature, self.pressure, self.inlet_mass_fraction
            self.inlet_mass_flow_rate = self.inlet_volumetric_flow_rate * self.gas.density

        return y0_matrix.flatten()

    def equations(self, t, y):
        pass

    def initial_condition(self):
        """
        Generate initial conditions for the selected model
        :return: Vector/Matrix representing the initial conditions
        """
        if self.solution_parser.resolution_method == ResolutionMethod.STEADYSTATE:
            return self.initial_condition_steady_state()

        return self.initial_condition_transient()

    def solve(self, tspan=None, time_ud=None):
        """
        Solve selected model
        :param tspan: Vector representing the integration time
        :param time_ud: Time unit dimension
        :return: Vector/Matrix representing the results
        """
        if self.solution_parser.resolution_method == ResolutionMethod.STEADYSTATE:
            y0 = self.initial_condition_steady_state()
            reactor_object = SteadyStatePseudoHomogeneous1DReactor(self.gas,
                                                                   self.surf,
                                                                   self.pressure,
                                                                   self.alfa,
                                                                   self.energy,
                                                                   self.inlet_mass_flow_rate,
                                                                   self.area)
            x, y = reactor_object.solve(self.numerical_solver,
                                        self.length,
                                        y0)

            self.solution_parser.x = x
            self.solution_parser.y = y
            self.solution_parser.is_solved = True
            return y

        if self.solution_parser.resolution_method == ResolutionMethod.TRANSIENT:
            y0 = self.initial_condition_transient()
            reactor_object = TransientPseudoHomogeneous1DReactor(self.gas,
                                                                 self.surf,
                                                                 self.pressure,
                                                                 self.alfa,
                                                                 self.energy,
                                                                 self.inlet_mass_flow_rate,
                                                                 self.area,
                                                                 self.gas_diffusion,
                                                                 self.inlet_temperature,
                                                                 self.inlet_mass_fraction,
                                                                 self.length,
                                                                 self.inert_specie_index,
                                                                 self.inert_coverage_index)
            x, y = reactor_object.solve(self.numerical_solver,
                                        self.uc.convert_to_seconds(tspan, time_ud),
                                        y0)

            self.solution_parser.x = x
            self.solution_parser.y = y
            self.solution_parser.length = self.length
            self.solution_parser.is_solved = True
            return y

        return None
