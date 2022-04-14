import numpy as np


class TransientPseudoHomogeneous1DReactor:
    def __init__(self,
                 gas,
                 surf,
                 pressure,
                 alfa,
                 energy,
                 inlet_mass_flow_rate,
                 area,
                 gas_diffusion,
                 inlet_temperature,
                 inlet_mass_fraction,
                 length,
                 inert_specie_index,
                 inert_coverage_index):
        """
        Class representing PseudoHomogeneous 1D reactor model
        :param gas: Cantera gas phase object
        :param surf: Cantera surface phase object
        :param pressure: Reactor pressure
        :param alfa: Reactor catalytic load
        :param energy: Enable/Disable energy balance
        :param inlet_mass_flow_rate: Reactor inlet mass flow rate
        :param area: Reactor section area
        :param gas_diffusion: Enable/Disable gas diffusion
        :param inlet_temperature: Inlet temperature
        :param inlet_mass_fraction: Inlet mass fraction
        :param length: Reactor length
        :param inert_specie_index: Inert specie index
        :param inert_coverage_index: Inert coverage index
        """
        self.gas = gas
        self.surf = surf
        self.pressure = pressure
        self.alfa = alfa
        self.energy = energy
        self.inlet_mass_flow_rate = inlet_mass_flow_rate
        self.area = area
        self.gas_diffusion = gas_diffusion

        self.inlet_temperature = inlet_temperature
        self.inlet_mass_fraction = inlet_mass_fraction
        self.length = length

        self.inert_specie_index = inert_specie_index
        self.inert_coverage_index = inert_coverage_index

        self.n_p = self.length.size
        self.n_v = self.gas.n_species + self.surf.n_species + 1
        self.alg = None

    def equations(self, t, y):
        """
        Function representing the Transient equations
        :param t: Independent variable - Time
        :param y: Dependent variable - Species composition, coverage and temperature as function of reactor length
        :return: Dependent variable variations based on independent variable
        """
        y_matrix = y.reshape(self.n_p, self.n_v)

        omega = y_matrix[:, :self.gas.n_species]
        z = y_matrix[:, self.gas.n_species:self.gas.n_species + self.surf.n_species]
        T = y_matrix[:, -1]

        r_gas = np.zeros([self.n_p, self.gas.n_species], dtype=np.float64)
        r_from_surface = np.zeros([self.n_p, self.gas.n_species], dtype=np.float64)
        r_surface = np.zeros([self.n_p, self.surf.n_species], dtype=np.float64)
        gas_mix_diff = np.zeros([self.n_p, self.gas.n_species], dtype=np.float64)
        q_gas = np.zeros([self.n_p], dtype=np.float64)
        q_from_surface = np.zeros([self.n_p], dtype=np.float64)
        gas_rho = np.zeros([self.n_p], dtype=np.float64)
        gas_cp = np.zeros([self.n_p], dtype=np.float64)
        gas_k = np.zeros([self.n_p], dtype=np.float64)

        for i in range(0, self.n_p):
            self.gas.TPY = T[i], self.pressure, np.maximum(0.0, omega[i, :])
            self.surf.TP = T[i], self.pressure
            self.surf.coverages = np.maximum(0.0, z[i, :])

            if self.gas.n_reactions > 0:
                r_gas[i, :] = self.gas.net_production_rates * self.gas.molecular_weights

            if self.surf.n_reactions > 0:
                r_from_surface[i, :] = self.surf.net_production_rates[:self.gas.n_species] * self.gas.molecular_weights
                r_surface[i, :] = self.surf.net_production_rates[-self.surf.n_species:]

            gas_rho[i] = self.gas.density
            gas_cp[i] = self.gas.cp_mass
            gas_k[i] = self.gas.thermal_conductivity / (self.gas.density * self.gas.cp_mass)

            diff_mix = self.gas.mix_diff_coeffs_mass
            diff_mix_zero = diff_mix == 0
            diff_mix[diff_mix_zero] = self.gas.binary_diff_coeffs[diff_mix_zero, diff_mix_zero]
            gas_mix_diff[i, :] = diff_mix

            if self.energy:
                if self.gas.n_reactions > 0:
                    q_gas[i] = -np.dot(self.gas.net_rates_of_progress, self.gas.delta_enthalpy)

                if self.surf.n_reactions > 0:
                    q_from_surface[i] = -np.dot(self.surf.net_rates_of_progress, self.surf.delta_enthalpy)

        domega = np.zeros_like(omega)
        dT = np.zeros_like(T)
        d1st_length_backward = self.length[1:-1] - self.length[:-2]
        d1st_length_forward = self.length[2:] - self.length[1:-1]
        d2nd_length = 0.5 * (self.length[2:] - self.length[:-2])

        # Inlet conditions
        domega[0, :] = self.inlet_mass_fraction - omega[0, :]
        if self.energy:
            dT[0] = self.inlet_temperature - T[0]

        # Outlet conditions
        if self.gas_diffusion:
            domega[-1, :] = omega[-1, :] - omega[-2, :]
        else:
            domega_outlet_derivative = (omega[-1, :] - omega[-2, :]) / (self.length[-1] - self.length[-2])
            domega[-1, :] = - self.inlet_mass_flow_rate * domega_outlet_derivative / (self.area * gas_rho[-1])
            domega[-1, :] = domega[-1, :] + r_gas[-1, :] / gas_rho[-1]
            domega[-1, :] = domega[-1, :] + self.alfa * r_from_surface[-1, :] / gas_rho[-1]

        if self.energy:
            if self.gas_diffusion:
                dT[-1] = T[-1] - T[-2]
            else:
                dT_outlet_derivative = (T[-1] - T[-2]) / (self.length[-1] - self.length[-2])
                dT[-1] = -(self.inlet_mass_flow_rate / (self.area * gas_rho[-1])) * dT_outlet_derivative
                dT[-1] = dT[-1] + q_gas[-1] / (gas_rho[-1] * gas_cp[-1])
                dT[-1] = dT[-1] + self.alfa * q_from_surface[-1] / (gas_rho[-1] * gas_cp[-1])

        # Equations for mass
        d1st_omega_backward = (omega[1:-1, :] - omega[:-2, :]) / d1st_length_backward[:, np.newaxis]
        domega[1:-1, :] = d1st_omega_backward / gas_rho[1:-1, np.newaxis]
        domega[1:-1, :] = - (self.inlet_mass_flow_rate / self.area) * domega[1:-1, :]
        domega[1:-1, :] = domega[1:-1, :] + (r_gas[1:-1, :] / gas_rho[1:-1, np.newaxis])
        domega[1:-1, :] = domega[1:-1, :] + self.alfa * (r_from_surface[1:-1, :] / gas_rho[1:-1, np.newaxis])

        if self.gas_diffusion:
            d1st_omega_forward = (omega[2:, :] - omega[1:-1, :]) / d1st_length_forward[:, np.newaxis]
            gas_diff_forward = 0.5 * (gas_mix_diff[2:, :] + gas_mix_diff[1:-1, :])
            gas_diff_backward = 0.5 * (gas_mix_diff[1:-1, :] + gas_mix_diff[:-2, :])
            domega[1:-1, :] = domega[1:-1, :] + gas_diff_forward * d1st_omega_forward / d2nd_length[:, np.newaxis]
            domega[1:-1, :] = domega[1:-1, :] - gas_diff_backward * d1st_omega_backward / d2nd_length[:, np.newaxis]

        # Inert specie
        domega[:, self.inert_specie_index] = 1. - np.sum(omega, axis=1)

        # Equations for site fraction
        dz = r_surface / self.surf.site_density

        # Inert specie
        dz[:, self.inert_coverage_index] = 1. - np.sum(z, axis=1)

        # Equations for temperature
        if self.energy:
            d1st_T_backward = (T[1:-1] - T[:-2]) / d1st_length_backward
            dT[1:-1] = -(self.inlet_mass_flow_rate / (self.area * gas_rho[1:-1])) * d1st_T_backward
            dT[1:-1] = dT[1:-1] + q_gas[1:-1] / (gas_rho[1:-1] * gas_cp[1:-1])
            dT[1:-1] = dT[1:-1] + self.alfa * q_from_surface[1:-1] / (gas_rho[1:-1] * gas_cp[1:-1])

            if self.gas_diffusion:
                d1st_T_forward = (T[2:] - T[1:-1]) / d1st_length_forward
                gas_k_forward = 0.5 * (gas_k[2:] + gas_k[1:-1]) / (gas_rho[1:-1] * gas_cp[1:-1])
                gas_k_backward = 0.5 * (gas_k[1:-1] + gas_k[:-2]) / (gas_rho[1:-1] * gas_cp[1:-1])
                dT[1:-1] = dT[1:-1] + gas_k_forward * d1st_T_forward / d2nd_length
                dT[1:-1] = dT[1:-1] - gas_k_backward * d1st_T_backward / d2nd_length

        dy_matrix = np.zeros(shape=y_matrix.shape, dtype=np.float64)

        dy_matrix[:, :self.gas.n_species] = domega
        dy_matrix[:, self.gas.n_species:self.gas.n_species + self.surf.n_species] = dz
        dy_matrix[:, -1] = dT

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

        # Inert species
        alg_matrix[:, self.inert_specie_index] = 0
        alg_matrix[:, self.gas.n_species + self.inert_coverage_index] = 0

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

    def solve(self, numerical_solver, tspan, initial_conditions):
        """
        Solve Transient equations
        :param numerical_solver: Numerical solver object
        :param tspan: Integration time
        :param initial_conditions: Reactor initial conditions
        :return: Time vector,
                 Matrix representing the solution in terms of composition, coverage and temperature as function of time and reactor length
        """
        self.alg = self.algebraic_equations()
        return numerical_solver.solve_dae(self.ode_equations,
                                          self.equations,
                                          self.residuals,
                                          initial_conditions,
                                          tspan,
                                          self.alg)
