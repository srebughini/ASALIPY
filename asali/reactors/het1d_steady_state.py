import copy

import numpy as np


class SteadyStateHeterogeneous1DReactor:
    def __init__(self,
                 gas,
                 surf,
                 pressure,
                 alfa,
                 energy,
                 inlet_mass_flow_rate,
                 gas_diffusion,
                 inlet_temperature,
                 inlet_mass_fraction,
                 length,
                 inert_specie_index,
                 inert_coverage_index,
                 reactor_shape_object,
                 solid):
        """
        Class representing Heterogeneous 1D reactor model
        :param gas: Cantera gas phase object
        :param surf: Cantera surface phase object
        :param pressure: Reactor pressure
        :param alfa: Reactor catalytic load
        :param energy: Enable/Disable energy balance
        :param inlet_mass_flow_rate: Reactor inlet mass flow rate
        :param gas_diffusion: Enable/Disable gas diffusion
        :param inlet_temperature: Inlet temperature
        :param inlet_mass_fraction: Inlet mass fraction
        :param length: Reactor length
        :param inert_specie_index: Inert specie index
        :param inert_coverage_index: Inert coverage index
        :param reactor_shape_object: Reactor shape class object
        :param solid: Solid class object
        """
        self.gas = gas
        self.surf = surf
        self.pressure = pressure
        self.alfa = alfa
        self.energy = energy
        self.inlet_mass_flow_rate = inlet_mass_flow_rate
        self.gas_diffusion = gas_diffusion

        self.inlet_temperature = inlet_temperature
        self.inlet_mass_fraction = inlet_mass_fraction
        self.length = length

        self.inert_specie_index = inert_specie_index
        self.inert_coverage_index = inert_coverage_index

        self.reactor_shape_object = reactor_shape_object
        self.solid = solid

        self.n_p = self.length.size
        self.n_s = self.gas.n_species
        self.n_surf = self.surf.n_species
        self.n_v = self.n_s + self.n_s + self.n_surf + 1 + 1
        self.alg = None

        self.n_loop = 5

    def estimate_integration_time(self):
        """
        Estimate integration time
        :return: Integration time
        """
        self.gas.TPY = self.inlet_temperature, self.pressure, self.inlet_mass_fraction
        return self.inlet_mass_flow_rate / (
                self.gas.density * self.reactor_shape_object.section_area * self.reactor_shape_object.void_fraction)

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
        q_gas = np.zeros([self.n_p], dtype=np.float64)
        q_surface = np.zeros([self.n_p], dtype=np.float64)
        gas_rho = np.zeros([self.n_p], dtype=np.float64)
        gas_cp = np.zeros([self.n_p], dtype=np.float64)
        gas_k = np.zeros([self.n_p], dtype=np.float64)
        gas_mu = np.zeros([self.n_p], dtype=np.float64)

        for i in range(0, self.n_p):
            self.gas.TPY = Tb[i], self.pressure, omegab[i, :]

            r_gas[i, :] = self.gas.net_production_rates * self.gas.molecular_weights

            gas_rho[i] = self.gas.density
            gas_mu[i] = self.gas.viscosity
            gas_cp[i] = self.gas.cp_mass
            gas_k[i] = self.gas.thermal_conductivity

            diff_mix = self.gas.mix_diff_coeffs_mass
            diff_mix_zero = diff_mix == 0
            diff_mix[diff_mix_zero] = self.gas.binary_diff_coeffs[diff_mix_zero, diff_mix_zero]
            gas_mix_diff[i, :] = diff_mix

            if self.energy:
                if self.gas.n_reactions > 0:
                    q_gas[i] = -np.dot(self.gas.net_rates_of_progress, self.gas.delta_enthalpy)

            self.gas.TPY = Tw[i], self.pressure, omegaw[i, :]
            self.surf.TP = Tw[i], self.pressure
            self.surf.coverages = z[i, :]
            r_from_surface[i, :] = self.surf.net_production_rates[:self.n_s] * self.gas.molecular_weights
            r_surface[i, :] = self.surf.net_production_rates[-self.n_surf:]

            if self.energy:
                if self.surf.n_reactions > 0:
                    q_surface[i] = -np.dot(self.surf.net_rates_of_progress, self.surf.delta_enthalpy)

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
                dTb[-1] = dTb[-1] + q_gas[-1] / (gas_rho[-1] * gas_cp[-1])
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
            dTb[1:-1] = dTb[1:-1] + q_gas[1:-1] / (gas_rho[1:-1] * gas_cp[1:-1])
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
            dTw[1:-1] = dTw[1:-1] + self.alfa * q_surface[1:-1] / (solid_cp * solid_rho * (1 - void_fraction))
            dTw[1:-1] = dTw[1:-1] + delta_T[1:-1] / (solid_cp * solid_rho * (1 - void_fraction))

        dy_matrix = np.zeros_like(y_matrix)

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

    def solve(self, numerical_solver, initial_conditions):
        """
        Solve Transient equations
        :param numerical_solver: Numerical solver object
        :param initial_conditions: Reactor initial conditions
        :return: Time vector,
                 Matrix representing the solution in terms of composition, coverage and temperature as function of time and reactor length
        """
        tspan = [0, 100 * self.estimate_integration_time()]

        n_p_vector = np.linspace(3, len(self.length), num=self.n_loop, endpoint=True, dtype=int)

        length_matrix = [self.length[np.linspace(0, len(self.length) - 1, num=n_p, dtype=int)] for n_p in n_p_vector]

        y0 = self.interpolate_variables_vector(self.length, initial_conditions, length_matrix[0])

        for i, n_p in enumerate(n_p_vector[:-1]):
            self.n_p = n_p
            self.length = length_matrix[i]
            self.alg = self.algebraic_equations()
            _, y = numerical_solver.solve_dae(self.ode_equations,
                                              self.equations,
                                              self.residuals,
                                              y0,
                                              tspan,
                                              self.alg)

            y0 = self.interpolate_variables_vector(self.length, y[-1, :], length_matrix[i + 1])

        self.length = length_matrix[-1]
        self.n_p = self.length.size
        self.alg = self.algebraic_equations()

        _, y = numerical_solver.solve_dae(self.ode_equations,
                                          self.equations,
                                          self.residuals,
                                          y0,
                                          tspan,
                                          self.alg)

        return self.length, y[-1, :].reshape(self.n_p, self.n_v)
