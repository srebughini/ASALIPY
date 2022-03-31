from asali.reactors.basic import BasicReactor

import numpy as np

from asali.reactors.het1d_steady_state import SteadyStateHeterogeneous1DReactor
from asali.reactors.het1d_transient import TransientHeterogeneous1DReactor
from asali.reactors.shapes.honeycomb import HoneyCombReactorShape
from asali.reactors.shapes.packed_bed import PackedBedReactorShape
from asali.utils.input_parser import ReactorModel, ReactorSection, ReactorType, InputParser, ResolutionMethod
from asali.utils.solid_material import SolidMaterial
from asali.reactors.shapes.tubular import TubularReactorShape


class Heterogeneous1DReactor(BasicReactor):
    def __init__(self, cantera_input_file, gas_phase_name, surface_phase_name):
        """
        Class representing Heterogeneous 1D reactor model
        :param cantera_input_file: Cantera file path
        :param gas_phase_name: Cantera gas phase name
        :param surface_phase_name: Cantera interface phase name
        """
        super().__init__(cantera_input_file=cantera_input_file, gas_phase_name=gas_phase_name,
                         surface_phase_name=surface_phase_name)
        self.solution_parser.reactor_type = ReactorType.HETEROGENEOUSPRF
        self.solid = SolidMaterial()
        self.reactor_shape_object = None

        self.is_mass_flow_rate = True
        self.gas_diffusion = False

        self.inlet_mass_flow_rate = 0.
        self.inlet_volumetric_flow_rate = 0.
        self.inlet_temperature = 0.
        self.inert_specie_index = 0
        self.inert_coverage_index = 0
        self.initial_solid_temperature = 0.

        self.alg = None
        self.length = None
        self.reactor_section = None
        self.inlet_mass_fraction = None
        self.inlet_mole_fraction = None

    def set_gas_diffusion(self, value):
        """
        Enable/Disable gas diffusion
        :param value: Variable to enable/disable gas diffusion
        :return: Bool for gas diffusion
        """
        self.gas_diffusion = InputParser.true_parser(value)

    def set_resolution_method(self, method):
        """
        Set resolution method
        :param method: Resolution method as string
        :return: ResolutionMethod object
        """
        self.solution_parser.resolution_method = InputParser.resolution_parser(method)
        return self.solution_parser.resolution_method

    def set_reactor_section(self, method):
        """
        Set reactor section
        :param method: Reactor section as string
        :return: ReactorSection
        """
        self.reactor_section = InputParser.section_parser(method)
        return self.reactor_section

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

    def set_solid_density(self, value, unit_dimension):
        """
        Set solid density
        :param value: Density value
        :param unit_dimension: Density unit dimension
        :return: Solid density in [kg/m3]
        """
        self.solid.density = (value, unit_dimension)
        return self.solid.density

    def set_solid_specific_heat(self, value, unit_dimension):
        """
        Set solid specific heat
        :param value: Specific heat value
        :param unit_dimension: Specific heat unit dimension
        :return: Solid specific heat in [J/kg/K]
        """
        self.solid.specific_heat = (value, unit_dimension)
        return self.solid.specific_heat

    def set_solid_thermal_conductivity(self, value, unit_dimension):
        """
        Set solid thermal conductivity
        :param value: Thermal conductivity value
        :param unit_dimension: Thermal conductivity unit dimension
        :return: Solid thermal conductivity in [W/m/K]
        """
        self.solid.thermal_conductivity = (value, unit_dimension)
        return self.solid.thermal_conductivity

    def set_initial_solid_temperature(self, value, unit_dimension):
        """
        Set initial solid temperature
        :param value: Temperature value
        :param unit_dimension: Temperature unit dimension
        :return: Initial solid temperature in [K]
        """
        self.initial_solid_temperature = self.uc.convert_to_kelvin(value, unit_dimension)
        return self.initial_solid_temperature

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

    def set_tubular_reactor(self, tube_diameter, tube_diameter_ud, wall_thickness, wall_thickness_ud):
        """
        Set Tubular reactor properties
        :param tube_diameter: Tube diameter
        :param tube_diameter_ud: Tube diameter unit dimensions
        :param wall_thickness: Wall thickness
        :param wall_thickness_ud: Wall thickness unit dimensions
        :return:
        """
        self.reactor_shape_object = TubularReactorShape()
        self.reactor_shape_object.set_properties(self.reactor_section,
                                                 tube_diameter,
                                                 tube_diameter_ud,
                                                 wall_thickness,
                                                 wall_thickness_ud)

        return [self.reactor_shape_object.tube_diameter,
                self.reactor_shape_object.void_fraction,
                self.reactor_shape_object.specific_area,
                self.reactor_shape_object.section_area]

    def set_honeycomb_reactor(self, cpsi, wall_thickness, wall_thickness_ud):
        """
        Set HONEYCOMB reactor properties
        :param cpsi: CPSI
        :param wall_thickness: Wall thickness
        :param wall_thickness_ud: Wall thickness unit dimensions
        :return: List of estimated reactor properties
        """
        self.reactor_shape_object = HoneyCombReactorShape()
        self.reactor_shape_object.set_properties(self.reactor_section,
                                                 cpsi,
                                                 wall_thickness,
                                                 wall_thickness_ud)

        return [self.reactor_shape_object.tube_diameter,
                self.reactor_shape_object.void_fraction,
                self.reactor_shape_object.specific_area,
                self.reactor_shape_object.section_area]

    def set_packed_bed_reactor(self, particle_diameter, particle_diameter_ud, tube_diameter, tube_diameter_ud,
                               void_fraction):
        """
        Set PACKEDBED reactor shape properties
        :param particle_diameter: Particle diameter
        :param particle_diameter_ud: Particle diameter unit dimensions
        :param tube_diameter: Tube diameter
        :param tube_diameter_ud: Tube diameter unit dimensions
        :param void_fraction: Void fraction
        :return: List of estimated reactor properties
        """

        self.reactor_shape_object = PackedBedReactorShape()
        self.reactor_shape_object.set_properties(particle_diameter, particle_diameter_ud, tube_diameter,
                                                 tube_diameter_ud,
                                                 void_fraction)

        return [self.reactor_shape_object.tube_diameter,
                self.reactor_shape_object.void_fraction,
                self.reactor_shape_object.particle_diameter,
                self.reactor_shape_object.specific_area,
                self.reactor_shape_object.section_area]

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

    def initial_condition_steady_state(self):
        """
        Function creating the initial condition of the Steady State solution
        :return: Matrix representing the initial mass fraction, coverage and temperature
        """
        n_p = self.length.size
        n_s = self.gas.n_species
        n_surf = self.surf.n_species
        n_v = n_s + n_s + n_surf + 1 + 1

        self.initial_mass_fraction = self.inlet_mass_fraction
        self.initial_temperature = self.inlet_temperature

        y0_matrix = np.zeros([n_p, n_v], dtype=np.float64)

        y0_matrix[:, :n_s] = self.initial_mass_fraction
        y0_matrix[:, n_s:n_s + n_s] = self.initial_mass_fraction
        y0_matrix[:, n_s + n_s:n_s + n_s + n_surf] = self.initial_coverage
        y0_matrix[:, -2] = self.initial_temperature
        y0_matrix[:, -1] = self.initial_solid_temperature

        if not self.energy:
            self.inlet_temperature = self.initial_temperature
            y0_matrix[:, -1] = self.initial_temperature

        y0_matrix[0, :n_s] = self.inlet_mass_fraction
        y0_matrix[0, -1] = self.inlet_temperature

        if not self.is_mass_flow_rate:
            self.gas.TPY = self.inlet_temperature, self.pressure, self.inlet_mass_fraction
            self.inlet_mass_flow_rate = self.inlet_volumetric_flow_rate * self.gas.density

        return y0_matrix.flatten()

    def initial_condition_transient(self):
        """
        Function creating the initial condition of the Transient solution
        :return: Vector/Matrix representing the initial conditions
        """
        n_p = self.length.size
        n_s = self.gas.n_species
        n_surf = self.surf.n_species
        n_v = n_s + n_s + n_surf + 1 + 1

        y0_matrix = np.zeros([n_p, n_v], dtype=np.float64)

        y0_matrix[:, :n_s] = self.initial_mass_fraction
        y0_matrix[:, n_s:n_s + n_s] = self.initial_mass_fraction
        y0_matrix[:, n_s + n_s:n_s + n_s + n_surf] = self.initial_coverage
        y0_matrix[:, -2] = self.initial_temperature
        y0_matrix[:, -1] = self.initial_solid_temperature

        if not self.energy:
            self.inlet_temperature = self.initial_temperature
            y0_matrix[:, -1] = self.initial_temperature

        y0_matrix[0, :n_s] = self.inlet_mass_fraction
        y0_matrix[0, -1] = self.inlet_temperature

        if not self.is_mass_flow_rate:
            self.gas.TPY = self.inlet_temperature, self.pressure, self.inlet_mass_fraction
            self.inlet_mass_flow_rate = self.inlet_volumetric_flow_rate * self.gas.density

        return y0_matrix.flatten()

    def solve(self, tspan=None, time_ud=None):
        """
        Solve selected model
        :param tspan: Vector representing the integration time
        :param time_ud: Time unit dimension
        :return: Vector/Matrix representing the results
        """
        if self.solution_parser.resolution_method == ResolutionMethod.STEADYSTATE:
            y0 = self.initial_condition_steady_state()
            reactor_object = SteadyStateHeterogeneous1DReactor(self.gas,
                                                               self.surf,
                                                               self.pressure,
                                                               self.alfa,
                                                               self.energy,
                                                               self.inlet_mass_flow_rate,
                                                               self.gas_diffusion,
                                                               self.inlet_temperature,
                                                               self.inlet_mass_fraction,
                                                               self.length,
                                                               self.inert_specie_index,
                                                               self.inert_coverage_index,
                                                               self.reactor_shape_object,
                                                               self.solid)
            x, y = reactor_object.solve(self.numerical_solver,
                                        y0)

            self.solution_parser.x = x
            self.solution_parser.y = y
            self.solution_parser.length = self.length
            self.solution_parser.is_solved = True
            return y

        if self.solution_parser.resolution_method == ResolutionMethod.TRANSIENT:
            y0 = self.initial_condition_transient()
            reactor_object = TransientHeterogeneous1DReactor(self.gas,
                                                             self.surf,
                                                             self.pressure,
                                                             self.alfa,
                                                             self.energy,
                                                             self.inlet_mass_flow_rate,
                                                             self.gas_diffusion,
                                                             self.inlet_temperature,
                                                             self.inlet_mass_fraction,
                                                             self.length,
                                                             self.inert_specie_index,
                                                             self.inert_coverage_index,
                                                             self.reactor_shape_object,
                                                             self.solid)
            x, y = reactor_object.solve(self.numerical_solver,
                                        self.uc.convert_to_seconds(tspan, time_ud),
                                        y0)

            self.solution_parser.x = x
            self.solution_parser.y = y
            self.solution_parser.length = self.length
            self.solution_parser.is_solved = True
            return y

        return None

    def get_solid_mass_fraction(self, index=None):
        """
        Get mass fraction of solid phase
        :param index: Index of the species to be extracted
        :return: Vector/Matrix representing the resulting mass fraction
        """
        mass_fraction = self.solution_parser.get_mass_fraction_wall()
        if index is None:
            return mass_fraction

        return mass_fraction[index]

    def get_solid_mole_fraction(self, index=None):
        """
        Get mole fraction of solid phase
        :param index: Index of the species to be extracted
        :return: Vector/Matrix representing the resulting mole fraction
        """
        mole_fraction = self.solution_parser.get_mole_fraction_wall()
        if index is None:
            return mole_fraction

        return mole_fraction[index]

    def get_solid_temperature(self, index=None):
        """
        Get temperature of solid phase
        :param index: Index of the axial/time position
        :return: Vector/Matrix representing the resulting temperature at a fixed axial/time position
        """
        temperature = self.solution_parser.get_temperature_wall()
        if index is None:
            return temperature

        return temperature[index]
