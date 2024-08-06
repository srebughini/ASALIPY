from abc import ABC

from asali.reactors.basic1d import Basic1DReactor
from asali.reactors.shapes.honeycomb import HoneyCombReactorShape
from asali.reactors.shapes.packed_bed import PackedBedReactorShape
from asali.utils.input_parser import ReactorType, InputParser
from asali.utils.solid_material import SolidMaterial
from asali.reactors.shapes.tubular import TubularReactorShape

import numpy as np


class Heterogeneous1DReactor(Basic1DReactor, ABC):
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
        self.initial_solid_temperature = 0.
        self.reactor_section = None

    def set_reactor_section(self, method):
        """
        Set reactor section
        :param method: Reactor section as string
        :return: ReactorSection
        """
        self.reactor_section = InputParser.section_parser(method)
        self._setup.reactor_section = self.reactor_section
        return self.reactor_section

    def set_solid_density(self, value, unit_dimension):
        """
        Set solid density
        :param value: Density value
        :param unit_dimension: Density unit dimension
        :return: Solid density in [kg/m3]
        """
        self.solid.density = (value, unit_dimension)
        self._setup.solid_density = self.solid.density
        return self.solid.density

    def set_solid_specific_heat(self, value, unit_dimension):
        """
        Set solid specific heat
        :param value: Specific heat value
        :param unit_dimension: Specific heat unit dimension
        :return: Solid specific heat in [J/kg/K]
        """
        self.solid.specific_heat = (value, unit_dimension)
        self._setup.solid_specific_heat = self.solid.specific_heat
        return self.solid.specific_heat

    def set_solid_thermal_conductivity(self, value, unit_dimension):
        """
        Set solid thermal conductivity
        :param value: Thermal conductivity value
        :param unit_dimension: Thermal conductivity unit dimension
        :return: Solid thermal conductivity in [W/m/K]
        """
        self.solid.thermal_conductivity = (value, unit_dimension)
        self._setup.solid_thermal_conductivity = self.solid.thermal_conductivity
        return self.solid.thermal_conductivity

    def set_initial_solid_temperature(self, value, unit_dimension):
        """
        Set initial solid temperature
        :param value: Temperature value
        :param unit_dimension: Temperature unit dimension
        :return: Initial solid temperature in [K]
        """
        self.initial_solid_temperature = self.uc.convert_to_kelvin(value, unit_dimension)
        self._setup.initial_solid_temperature = self.initial_solid_temperature
        return self.initial_solid_temperature

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

        self._setup.tube_diameter = self.reactor_shape_object.tube_diameter
        self._setup.void_fraction = self.reactor_shape_object.void_fraction
        self._setup.specific_area = self.reactor_shape_object.specific_area
        self._setup.section_area = self.reactor_shape_object.section_area

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

        self._setup.tube_diameter = self.reactor_shape_object.tube_diameter
        self._setup.void_fraction = self.reactor_shape_object.void_fraction
        self._setup.specific_area = self.reactor_shape_object.specific_area
        self._setup.section_area = self.reactor_shape_object.section_area

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

        self._setup.tube_diameter = self.reactor_shape_object.tube_diameter
        self._setup.void_fraction = self.reactor_shape_object.void_fraction
        self._setup.specific_area = self.reactor_shape_object.specific_area
        self._setup.section_area = self.reactor_shape_object.section_area
        self._setup.particle_diameter = self.reactor_shape_object.particle_diameter

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

    def equations(self, t, y):
        """
        Function representing the Reactor model equations
        :param t: Independent variable - Time
        :param y: Dependent variable - Species composition, coverage and temperature as function of reactor length
        :return: Dependent variable variations based on independent variable
        """
        # Extraction of geometrical properties
        void_fraction = self.reactor_shape_object.void_fraction
        area = self.reactor_shape_object.section_area
        specific_area = self.reactor_shape_object.specific_area

        # Extraction of solid properties
        solid_k = self.solid.thermal_conductivity
        solid_cp = self.solid.specific_heat
        solid_rho = self.solid.density

        y_matrix = y.reshape(self.n_p, self.n_v)

        omegab = y_matrix[:, :self.n_s]
        omegaw = y_matrix[:, self.n_s:self.n_s + self.n_s]
        z = y_matrix[:, self.n_s + self.n_s:self.n_s + self.n_s + self.n_surf]
        Tb = y_matrix[:, -2]
        Tw = y_matrix[:, -1]

        r_gas = np.zeros([self.n_p, self.n_s], dtype=np.float64)
        r_from_surface = np.zeros([self.n_p, self.n_s], dtype=np.float64)
        r_surface = np.zeros([self.n_p, self.n_surf], dtype=np.float64)
        gas_mix_diff = np.zeros([self.n_p, self.n_s], dtype=np.float64)
        q_from_gas = np.zeros([self.n_p], dtype=np.float64)
        q_from_surface = np.zeros([self.n_p], dtype=np.float64)
        gas_rho = np.zeros([self.n_p], dtype=np.float64)
        gas_cp = np.zeros([self.n_p], dtype=np.float64)
        gas_k = np.zeros([self.n_p], dtype=np.float64)
        gas_mu = np.zeros([self.n_p], dtype=np.float64)

        for i in range(0, self.n_p):
            self.gas.TPY = Tb[i], self.pressure, omegab[i, :]

            r_gas[i, :] = self.get_homogeneous_gas_species_reaction_rates() * self.gas.molecular_weights

            gas_rho[i] = self.gas.density
            gas_mu[i] = self.gas.viscosity
            gas_cp[i] = self.gas.cp_mass
            gas_k[i] = self.gas.thermal_conductivity

            diff_mix = self.gas.mix_diff_coeffs_mass
            diff_mix_zero = diff_mix == 0
            diff_mix[diff_mix_zero] = self.gas.binary_diff_coeffs[diff_mix_zero, diff_mix_zero]
            gas_mix_diff[i, :] = diff_mix

            if self.energy:
                q_from_gas[i] = self.get_homogeneous_heat_of_reaction()

            self.gas.TPY = Tw[i], self.pressure, omegaw[i, :]
            self.surf.TP = Tw[i], self.pressure
            self.surf.coverages = z[i, :]
            r_from_surface[i, :] = self.get_heterogeneous_gas_species_reaction_rates() * self.gas.molecular_weights
            r_surface[i, :] = self.get_surface_species_reaction_rates()

            if self.energy:
                q_from_surface[i] = self.get_heterogeneous_heat_of_reaction()

        k_mat = self.estimate_mass_transfer_coefficient(gas_mu, gas_rho, gas_mix_diff)
        k_heat = self.estimate_heat_transfer_coefficient(gas_mu, gas_k, gas_cp)

        domegab = np.zeros_like(omegab)
        dTb = np.zeros_like(Tb)
        dTw = np.zeros_like(Tw)

        delta_omega = specific_area * (gas_rho * (k_mat * (omegab - omegaw)).T).T
        delta_T = k_heat * specific_area * (Tb - Tw)
        d1st_length_backward = self.length[1:-1] - self.length[:-2]
        d1st_length_forward = self.length[2:] - self.length[1:-1]
        d2nd_length = 0.5 * (self.length[2:] - self.length[:-2])

        # Inlet conditions
        domegab[0, :] = self.inlet_mass_fraction - omegab[0, :]
        if self.energy:
            dTb[0] = self.inlet_temperature - Tb[0]
            dTw[0] = Tw[1] - Tw[0]

        # Outlet conditions
        if self.gas_diffusion:
            domegab[-1, :] = omegab[-1, :] - omegab[-2, :]
        else:
            d1st_omegab_outlet = (omegab[-1, :] - omegab[-2, :]) / (self.length[-1] - self.length[-2])
            domegab[-1, :] = - self.inlet_mass_flow_rate * d1st_omegab_outlet / (area * gas_rho[-1])
            domegab[-1, :] = domegab[-1, :] + r_gas[-1, :] / gas_rho[-1]
            domegab[-1, :] = domegab[-1, :] - delta_omega[-1, :] / (void_fraction * gas_rho[-1])

        if self.energy:
            if self.gas_diffusion:
                dTb[-1] = Tb[-1] - Tb[-2]
            else:
                d1st_Tb_outlet = (Tb[-1] - Tb[-2]) / (self.length[-1] - self.length[-2])
                dTb[-1] = -(self.inlet_mass_flow_rate / (area * gas_rho[-1])) * d1st_Tb_outlet
                dTb[-1] = dTb[-1] + q_from_gas[-1] / (gas_rho[-1] * gas_cp[-1])
                dTb[-1] = dTb[-1] - delta_T[-1] / (void_fraction * gas_rho[-1] * gas_cp[-1])

            dTw[-1] = Tw[-1] - Tw[-2]

        # Equations for BULK mass
        d1st_omegab_backward = (omegab[1:-1, :] - omegab[:-2, :]) / d1st_length_backward[:, np.newaxis]
        domegab[1:-1, :] = - (self.inlet_mass_flow_rate / area) * d1st_omegab_backward
        domegab[1:-1, :] = domegab[1:-1, :] + r_gas[1:-1, :]
        domegab[1:-1, :] = domegab[1:-1, :] - delta_omega[1:-1, :] / void_fraction
        domegab[1:-1, :] = domegab[1:-1, :] / gas_rho[1:-1, np.newaxis]

        if self.gas_diffusion:
            d1st_omegab_forward = (omegab[2:, :] - omegab[1:-1, :]) / d1st_length_forward[:, np.newaxis]
            gas_diff_forward = 0.5 * (gas_mix_diff[2:, :] + gas_mix_diff[1:-1, :])
            gas_diff_backward = 0.5 * (gas_mix_diff[1:-1, :] + gas_mix_diff[:-2, :])
            domegab[1:-1, :] = domegab[1:-1, :] + gas_diff_forward * d1st_omegab_forward / d2nd_length[:,
                                                                                           np.newaxis]
            domegab[1:-1, :] = domegab[1:-1, :] - gas_diff_backward * d1st_omegab_backward / d2nd_length[:,
                                                                                             np.newaxis]

        # Equations of WALL mass
        domega_wall = delta_omega * void_fraction + self.alfa * void_fraction * r_from_surface

        # Inert specie
        domegab[:, self.inert_specie_index] = 1. - np.sum(omegab, axis=1)
        domega_wall[:, self.inert_specie_index] = 1. - np.sum(omegaw, axis=1)

        # Equations for site fraction
        dz = r_surface / self.surf.site_density

        # Inert specie
        dz[:, self.inert_coverage_index] = 1. - np.sum(z, axis=1)

        if self.energy:
            # Equations of BULK energy
            d1st_Tb_backward = (Tb[1:-1] - Tb[:-2]) / d1st_length_backward
            dTb[1:-1] = -(self.inlet_mass_flow_rate / (area * gas_rho[1:-1])) * d1st_Tb_backward
            dTb[1:-1] = dTb[1:-1] + q_from_gas[1:-1] / (gas_rho[1:-1] * gas_cp[1:-1])
            dTb[1:-1] = dTb[1:-1] - delta_T[1:-1] / (void_fraction * gas_rho[1:-1] * gas_cp[1:-1])

            if self.gas_diffusion:
                d1st_Tb_forward = (Tb[2:] - Tb[1:-1]) / d1st_length_forward
                gas_k_forward = 0.5 * (gas_k[2:] + gas_k[1:-1]) / (gas_rho[1:-1] * gas_cp[1:-1])
                gas_k_backward = 0.5 * (gas_k[1:-1] + gas_k[:-2]) / (gas_rho[1:-1] * gas_cp[1:-1])
                dTb[1:-1] = dTb[1:-1] + gas_k_forward * d1st_Tb_forward / d2nd_length
                dTb[1:-1] = dTb[1:-1] - gas_k_backward * d1st_Tb_backward / d2nd_length

            # Equations of WALL energy
            d1st_Tw_forward = (Tw[2:] - Tw[1:-1]) / d1st_length_forward
            d1st_Tw_backward = (Tw[1:-1] - Tw[:-2]) / d1st_length_backward

            dTw[1:-1] = (solid_k / (solid_cp * solid_rho)) * (d1st_Tw_forward - d1st_Tw_backward) / d2nd_length
            dTw[1:-1] = dTw[1:-1] + self.alfa * q_from_surface[1:-1] / (solid_cp * solid_rho * (1 - void_fraction))
            dTw[1:-1] = dTw[1:-1] + delta_T[1:-1] / (solid_cp * solid_rho * (1 - void_fraction))

        dy_matrix = np.zeros(shape=y_matrix.shape, dtype=np.float64)
        dy_matrix[:, :self.n_s] = domegab
        dy_matrix[:, self.n_s:self.n_s + self.n_s] = domega_wall
        dy_matrix[:, self.n_s + self.n_s:self.n_s + self.n_s + self.n_surf] = dz

        if self.energy:
            dy_matrix[:, -2] = dTb
            dy_matrix[:, -1] = dTw

        return dy_matrix.flatten()

    def algebraic_equations(self):
        """
        Generate the vector describing algebraic (0) and differential (1) equations
        :return: Vector on 0/1 describing algebraic/differential equations
        """
        alg_matrix = np.ones([self.n_p, self.n_v], dtype=int)

        # Inlet conditions
        alg_matrix[0, :self.gas.n_species] = 0
        if self.energy:
            alg_matrix[0, -1] = 0

        # Outlet conditions
        if self.gas_diffusion:
            alg_matrix[-1, :self.gas.n_species] = 0
            if self.energy:
                alg_matrix[-1, -1] = 0

        # Inlet conditions
        alg_matrix[0, :self.n_s] = 0
        if self.energy:
            alg_matrix[0, -2] = 0
            alg_matrix[0, -1] = 0

        # Equations of WALL mass
        alg_matrix[:, self.n_s:self.n_s + self.n_s] = 0

        # Inert species
        alg_matrix[:, self.inert_specie_index] = 0
        alg_matrix[:, self.n_s + self.inert_specie_index] = 0
        alg_matrix[:, self.n_s + self.n_s + self.inert_coverage_index] = 0

        # Outlet conditions
        if self.gas_diffusion:
            alg_matrix[-1, :self.n_s] = 0

        if self.energy:
            if self.gas_diffusion:
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

    def ode_equations(self, t, y):
        """
        Function representing the ODE system to estimate the DAE initial conditions
        :param t: Independent variable - Time
        :param y: Dependent variable - Species composition, coverage and temperature as function of reactor length
        :return: Dependent variable variations based on independent variable
        """
        return self.equations(t, y) * np.fabs(np.round(self.alg - 1))

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
