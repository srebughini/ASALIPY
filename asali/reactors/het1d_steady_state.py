import numpy as np

from asali.reactors.het1d import Heterogeneous1DReactor
from asali.utils.input_parser import ResolutionMethod


class SteadyStateHeterogeneous1DReactor(Heterogeneous1DReactor):
    def __init__(self, cantera_input_file, gas_phase_name, surface_phase_name):
        """
        Class representing Steady State pseudoHomogeneous 1D reactor model
        :param cantera_input_file: Cantera file path
        :param gas_phase_name: Cantera gas phase name
        :param surface_phase_name: Cantera interface phase name
        """
        super().__init__(cantera_input_file=cantera_input_file, gas_phase_name=gas_phase_name,
                         surface_phase_name=surface_phase_name)

        self.solution_parser.resolution_method = ResolutionMethod.STEADYSTATE

        self.n_s = self.gas.n_species
        self.n_surf = self.surf.n_species
        self.n_v = self.n_s + self.n_s + self.n_surf + 1 + 1
        self.n_loop = 5

    def estimate_integration_time(self):
        """
        Estimate integration time
        :return: Integration time
        """
        self.gas.TPY = self.inlet_temperature, self.pressure, self.inlet_mass_fraction
        if self.is_mass_flow_rate:
            return self.inlet_mass_flow_rate / (
                    self.gas.density * self.reactor_shape_object.section_area * self.reactor_shape_object.void_fraction)

        self.inlet_mass_flow_rate = self.inlet_volumetric_flow_rate * self.gas.density
        return self.inlet_mass_flow_rate / (
                self.gas.density * self.reactor_shape_object.section_area * self.reactor_shape_object.void_fraction)


    def interpolate_variables_vector(self, xp, fp, x):
        """
        Interpolated vector variables on a new length
        :param xp: Original length vector
        :param fp: Original variables vector
        :param x: New length vector
        :return: Interpolated variables vector
        """
        fp_matrix = fp.reshape(-1, self.n_v)
        f_matrix = np.zeros([len(x), self.n_v], dtype=np.float64)

        for i in range(0, self.n_v):
            f_matrix[:, i] = np.interp(x, xp, fp_matrix[:, i])

        return f_matrix.flatten()

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
            domegab[-1, :] = domegab[-1, :] - delta_omega[-1, :] / void_fraction

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
        domegab[1:-1, :] = d1st_omegab_backward / gas_rho[1:-1, np.newaxis]
        domegab[1:-1, :] = - (self.inlet_mass_flow_rate / area) * domegab[1:-1, :]
        domegab[1:-1, :] = domegab[1:-1, :] + (r_gas[1:-1, :] / gas_rho[1:-1, np.newaxis])
        domegab[1:-1, :] = domegab[1:-1, :] - delta_omega[1:-1, :] / void_fraction

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

    def initial_condition(self):
        """
        Function creating the initial condition of the Steady State solution
        :return: Matrix representing the initial mass fraction, coverage and temperature
        """
        self.initial_mass_fraction = self.inlet_mass_fraction
        self.initial_temperature = self.inlet_temperature

        y0_matrix = np.zeros([self.n_p, self.n_v], dtype=np.float64)

        y0_matrix[:, :self.n_s] = self.initial_mass_fraction
        y0_matrix[:, self.n_s:self.n_s + self.n_s] = self.initial_mass_fraction
        y0_matrix[:, self.n_s + self.n_s:self.n_s + self.n_s + self.n_surf] = self.initial_coverage
        y0_matrix[:, -2] = self.initial_temperature
        y0_matrix[:, -1] = self.initial_solid_temperature

        if not self.energy:
            self.inlet_temperature = self.initial_temperature
            y0_matrix[:, -1] = self.initial_temperature

        y0_matrix[0, :self.n_s] = self.inlet_mass_fraction
        y0_matrix[0, -1] = self.inlet_temperature

        if not self.is_mass_flow_rate:
            self.gas.TPY = self.inlet_temperature, self.pressure, self.inlet_mass_fraction
            self.inlet_mass_flow_rate = self.inlet_volumetric_flow_rate * self.gas.density

        return y0_matrix.flatten()

    def solve(self, tspan=None, time_ud=None):
        """
        Solve selected model
        :param tspan: Not required
        :param time_ud: Not required
        :return: Vector/Matrix representing the results
        """
        tspan = [0, 100 * self.estimate_integration_time()]

        n_p_vector = np.linspace(3, len(self.length), num=self.n_loop, endpoint=True, dtype=int)

        length_matrix = [self.length[np.linspace(0, len(self.length) - 1, num=n_p, dtype=int)] for n_p in n_p_vector]

        y0 = self.interpolate_variables_vector(self.length, self.initial_condition(), length_matrix[0])

        for i, n_p in enumerate(n_p_vector[:-1]):
            self.n_p = n_p
            self.length = length_matrix[i]
            self.alg = self.algebraic_equations()
            _, y = self.numerical_solver.solve_dae(self.ode_equations,
                                                   self.equations,
                                                   self.residuals,
                                                   y0,
                                                   tspan,
                                                   self.alg)

            y0 = self.interpolate_variables_vector(self.length, y[-1, :], length_matrix[i + 1])

        self.length = length_matrix[-1]
        self.n_p = self.length.size
        self.alg = self.algebraic_equations()

        _, y = self.numerical_solver.solve_dae(self.ode_equations,
                                               self.equations,
                                               self.residuals,
                                               y0,
                                               tspan,
                                               self.alg)

        self.solution_parser.x = self.length
        self.solution_parser.y = y[-1, :].reshape(self.n_p, self.n_v)
        self.solution_parser.length = self.length
        self.solution_parser.is_solved = True
        return self.solution_parser.y
