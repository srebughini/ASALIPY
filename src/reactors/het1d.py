from asali.reactors.basic import BasicReactor

import numpy as np

from asali.reactors.shapes.honeycomb import HoneyCombReactorShape
from asali.reactors.shapes.packed_bed import PackedBedReactorShape
from asali.utils.input_parser import ReactorModel, ReactorSection, ReactorType, InputParser
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

    def estimate_mass_transfer_coefficient(self, viscosity, density, diffusivity):
        """
        Estimate mass transfer coefficient based on ReactorModel
        :param viscosity: Gas viscosity in [Pas]
        :param density: Gas density in [kg/m3]
        :param diffusivity: Gas mixture diffusivity in [m2/s]
        :return: Mass transfer coefficient
        """

        self.reactor_shape_object.length = self.length
        return self.reactor_shape_object.estimate_mass_transfer_coefficient(self.inlet_mass_flow_rate,
                                                                            viscosity,
                                                                            density,
                                                                            diffusivity)

    def estimate_heat_transfer_coefficient(self, viscosity, conductivity, specific_heat):
        """
        Estimate heat transfer coefficient based on ReactorModel
        :param viscosity: Gas viscosity in [Pas]
        :param conductivity: Gas thermal conductivity in [W/m/K]
        :param specific_heat: Gas specific heat in [J/kg/K]
        :return: Heat transfer coefficient
        """
        self.reactor_shape_object.length = self.length
        return self.reactor_shape_object.estimate_heat_transfer_coefficient(self.inlet_mass_flow_rate,
                                                                            viscosity,
                                                                            conductivity,
                                                                            specific_heat)

    def ode_equations(self, t, y):
        """
        Function representing the ODE system to estimate the DAE initial conditions
        :param t: Independent variable - Time
        :param y: Dependent variable - Species composition, coverage and temperature as function of reactor length
        :return: Dependent variable variations based on independent variable
        """
        return self.equations(t, y) * np.fabs(np.round(self.alg - 1))

    def algebraic_equations(self):
        """
        Generate the vector describing algebraic (0) and differential (1) equations
        :return: Vector on 0/1 describing algebraic/differential equations
        """
        NP = self.length.size
        NV = self.gas.n_species + self.gas.n_species + self.surf.n_species + 1 + 1

        alg_matrix = np.ones([NP, NV], dtype=int)

        alg_matrix[0, :self.gas.n_species] = 0

        alg_matrix[:, self.gas.n_species:self.gas.n_species + self.gas.n_species] = 0

        alg_matrix[:, self.inert_specie_index] = 0
        alg_matrix[:, self.gas.n_species + self.inert_specie_index] = 0
        alg_matrix[:, self.gas.n_species + self.gas.n_species + self.inert_coverage_index] = 0

        alg_matrix[0, -2] = 0
        alg_matrix[0, -1] = 0
        alg_matrix[-1, -2] = 0
        alg_matrix[-1, -1] = 0

        return alg_matrix.flatten()

    def residuals(self, t, y, dy):
        """
        Residuals required by the DAE solver
        :param t: Independent variable - Time
        :param y: Dependent variable - Species composition, coverage and temperature as function of reactor length
        :param dy: Dependent variable variations based on independent variable
        :return: Residuals - y - dy
        """
        res = self.equations(t, y)
        diff_mask = self.alg == 1
        res[diff_mask] = res[diff_mask] - dy[diff_mask]
        return res

    def equations(self, t, y):
        """
        Function representing the Reactor model equations
        :param t: Independent variable - Time
        :param y: Dependent variable - Species composition, coverage and temperature as function of reactor length
        :return: Dependent variable variations based on independent variable
        """
        NP = self.length.size
        NV = self.gas.n_species + self.gas.n_species + self.surf.n_species + 1 + 1

        # Extraction of geometrical properties
        void_fraction = self.reactor_shape_object.void_fraction
        area = self.reactor_shape_object.section_area
        specific_area = self.reactor_shape_object.specific_area

        # Extraction of solid properties
        solid_k = self.solid.thermal_conductivity
        solid_cp = self.solid.specific_heat
        solid_rho = self.solid.density

        y_matrix = y.reshape(NP, NV)

        omega_bulk = y_matrix[:, :self.gas.n_species]
        omega_wall = y_matrix[:, self.gas.n_species:self.gas.n_species + self.gas.n_species]
        z = y_matrix[:,
            self.gas.n_species + self.gas.n_species:self.gas.n_species + self.gas.n_species + self.surf.n_species]
        T_bulk = y_matrix[:, -2]
        T_wall = y_matrix[:, -1]

        gas_reaction_rates = np.zeros([NP, self.gas.n_species], dtype=np.float64)
        gas_reaction_rates_from_surface = np.zeros([NP, self.gas.n_species], dtype=np.float64)
        coverage_reaction_rates = np.zeros([NP, self.surf.n_species], dtype=np.float64)
        mix_diff_coeffs_mass = np.zeros([NP, self.gas.n_species], dtype=np.float64)
        heat_of_reaction_from_gas = np.zeros([NP], dtype=np.float64)
        heat_from_reaction_from_surface = np.zeros([NP], dtype=np.float64)
        density = np.zeros([NP], dtype=np.float64)
        cp_mass = np.zeros([NP], dtype=np.float64)
        thermal_conductivity = np.zeros([NP], dtype=np.float64)
        viscosity = np.zeros([NP], dtype=np.float64)

        for i in range(0, NP):
            self.gas.TPY = T_bulk[i], self.pressure, omega_bulk[i, :]

            gas_reaction_rates[i, :] = self.gas.net_production_rates * self.gas.molecular_weights

            density[i] = self.gas.density
            viscosity[i] = self.gas.viscosity
            cp_mass[i] = self.gas.cp_mass
            thermal_conductivity[i] = self.gas.thermal_conductivity

            diff_mix = self.gas.mix_diff_coeffs_mass
            diff_mix_zero = diff_mix == 0
            diff_mix[diff_mix_zero] = self.gas.binary_diff_coeffs[diff_mix_zero, diff_mix_zero]
            mix_diff_coeffs_mass[i, :] = diff_mix

            if self.energy:
                if self.gas.n_reactions > 0:
                    heat_of_reaction_from_gas[i] = -np.dot(self.gas.net_rates_of_progress, self.gas.delta_enthalpy)

            self.gas.TPY = T_wall[i], self.pressure, omega_wall[i, :]
            self.surf.TP = T_wall[i], self.pressure
            self.surf.coverages = z[i, :]
            gas_reaction_rates_from_surface[i, :] = self.surf.net_production_rates[
                                                    :self.gas.n_species] * self.gas.molecular_weights
            coverage_reaction_rates[i, :] = self.surf.net_production_rates[-self.surf.n_species:]

            if self.energy:
                if self.surf.n_reactions > 0:
                    heat_from_reaction_from_surface[i] = -np.dot(self.surf.net_rates_of_progress,
                                                                 self.surf.delta_enthalpy)

        k_mat = self.estimate_mass_transfer_coefficient(viscosity, density, mix_diff_coeffs_mass)
        k_heat = self.estimate_heat_transfer_coefficient(viscosity, thermal_conductivity, cp_mass)

        domega_bulk = np.zeros_like(omega_bulk)
        dT_bulk = np.zeros_like(T_bulk)
        dT_wall = np.zeros_like(T_wall)

        domega_bulk[0, :] = self.inlet_mass_fraction - omega_bulk[0, :]  # Inlet conditions
        dT_bulk[0] = self.inlet_temperature - T_bulk[0]  # Inlet conditions
        dT_bulk[-1] = T_bulk[-1] - T_bulk[-2]  # Outlet conditions

        dT_wall[0] = T_wall[1] - T_wall[0]  # Inlet conditions
        dT_wall[-1] = T_wall[-1] - T_wall[-2]  # Outlet conditions

        delta_omega = specific_area * (density * (k_mat * (omega_bulk - omega_wall)).T).T

        derivative_1st_omega_backward = (
                (omega_bulk[1:, :] - omega_bulk[:-1, :]).T / (self.length[1:] - self.length[:-1])).T

        domega_bulk[1:, :] = - (
                derivative_1st_omega_backward.T * (self.inlet_mass_flow_rate / (area * density[1:]))).T
        domega_bulk[1:, :] = domega_bulk[1:, :] + (gas_reaction_rates[1:, :].T / density[1:]).T
        domega_bulk[1:, :] = domega_bulk[1:, :] - delta_omega[1:, :] / void_fraction

        # Equations of wall interface
        domega_wall = delta_omega * void_fraction + self.alfa * void_fraction * gas_reaction_rates_from_surface

        # Inert specie
        domega_bulk[:, self.inert_specie_index] = 1. - np.sum(omega_bulk, axis=1)
        domega_wall[:, self.inert_specie_index] = 1. - np.sum(omega_wall, axis=1)

        # Equations for site fraction
        dz = coverage_reaction_rates / self.surf.site_density

        # Inert specie
        dz[:, self.inert_coverage_index] = 1. - np.sum(z, axis=1)

        if self.energy:
            derivative_1st_temperature_bulk_forward = (T_bulk[2:] - T_bulk[1:-1]) / (
                    self.length[2:] - self.length[1:-1])
            derivative_1st_temperature_bulk_backward = (T_bulk[1:-1] - T_bulk[:-2]) / (
                    self.length[1:-1] - self.length[:-2])

            thermal_coefficient_forward = 0.5 * (thermal_conductivity[2:] + thermal_conductivity[1:-1])
            thermal_coefficient_backward = 0.5 * (thermal_conductivity[1:-1] + thermal_conductivity[:-2])

            delta_T = k_heat * specific_area * (T_bulk - T_wall)

            dT_bulk[1:-1] = -(self.inlet_mass_flow_rate / (
                    area * density[1:-1])) * derivative_1st_temperature_bulk_backward
            dT_bulk[1:-1] = dT_bulk[1:-1] + (
                    thermal_coefficient_forward * derivative_1st_temperature_bulk_forward - thermal_coefficient_backward * derivative_1st_temperature_bulk_backward) / (
                                    0.5 * (self.length[2:] - self.length[:-2]) * density[1:-1] * cp_mass[1:-1])
            dT_bulk[1:-1] = dT_bulk[1:-1] + heat_of_reaction_from_gas[1:-1] / (density[1:-1] * cp_mass[1:-1])
            dT_bulk[1:-1] = dT_bulk[1:-1] - delta_T[1:-1] / (void_fraction * density[1:-1] * cp_mass[1:-1])

            derivative_1st_temperature_wall_forward = (T_wall[2:] - T_wall[1:-1]) / (
                    self.length[2:] - self.length[1:-1])
            derivative_1st_temperature_wall_backward = (T_wall[1:-1] - T_wall[:-2]) / (
                    self.length[1:-1] - self.length[:-2])

            dT_wall[1:-1] = (solid_k / (solid_cp * solid_rho)) * (
                    derivative_1st_temperature_wall_forward - derivative_1st_temperature_wall_backward) / (
                                    0.5 * (self.length[2:] - self.length[:-2]))
            dT_wall[1:-1] = dT_wall[1:-1] + self.alfa * heat_from_reaction_from_surface[1:-1] / (
                    solid_cp * solid_rho * (1 - void_fraction))
            dT_wall[1:-1] = dT_wall[1:-1] + delta_T[1:-1] / (solid_cp * solid_rho * (1 - void_fraction))

        dy_matrix = np.zeros_like(y_matrix)

        dy_matrix[:, :self.gas.n_species] = domega_bulk

        dy_matrix[:, self.gas.n_species:self.gas.n_species + self.gas.n_species] = domega_wall
        dy_matrix[:,
        self.gas.n_species + self.gas.n_species:self.gas.n_species + self.gas.n_species + self.surf.n_species] = dz

        if self.energy:
            dy_matrix[:, -2] = dT_bulk
            dy_matrix[:, -1] = dT_wall

        return dy_matrix.flatten()

    def initial_condition(self):
        """
        Generate initial conditions for Reactor model
        :return: Vector/Matrix representing the initial conditions
        """
        NP = self.length.size
        NV = self.gas.n_species + self.gas.n_species + self.surf.n_species + 1 + 1

        y0_matrix = np.zeros([NP, NV], dtype=np.float64)

        y0_matrix[:, :self.gas.n_species] = self.initial_mass_fraction
        y0_matrix[:, self.gas.n_species:self.gas.n_species + self.gas.n_species] = self.initial_mass_fraction
        y0_matrix[:,
        self.gas.n_species + self.gas.n_species:self.gas.n_species + self.gas.n_species + self.surf.n_species] = self.initial_coverage
        y0_matrix[:, -2] = self.initial_temperature
        y0_matrix[:, -1] = self.initial_solid_temperature

        if not self.energy:
            self.inlet_temperature = self.initial_temperature
            y0_matrix[:, -1] = self.initial_temperature

        y0_matrix[0, :self.gas.n_species] = self.inlet_mass_fraction
        y0_matrix[0, -1] = self.inlet_temperature

        if not self.is_mass_flow_rate:
            self.gas.TPY = self.inlet_temperature, self.pressure, self.inlet_mass_fraction
            self.inlet_mass_flow_rate = self.inlet_volumetric_flow_rate * self.gas.density

        return y0_matrix.flatten()

    def solve(self, tspan, time_ud):
        self.alg = self.algebraic_equations()
        x, y = self.numerical_solver.solve_dae(self.ode_equations,
                                               self.equations,
                                               self.residuals,
                                               self.initial_condition(),
                                               self.uc.convert_to_seconds(tspan, time_ud),
                                               self.alg)

        self.solution_parser.x = x
        self.solution_parser.y = y
        self.solution_parser.is_solved = True
        self.solution_parser.length = self.length
        return y

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
