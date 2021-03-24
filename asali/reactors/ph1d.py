from asali.reactors.basic import ReactorType, BasicReactor, ResolutionMethod

import numpy as np


class PseudoHomogeneous1DReactor(BasicReactor):
    def __init__(self, cantera_input_file, gas_phase_name, surface_phase_name):
        super().__init__(cantera_input_file=cantera_input_file, gas_phase_name=gas_phase_name, surface_phase_name=surface_phase_name)
        self.reactor_type = ReactorType.PSEUDOHOMOGENEOUSPFR
        self.alg = None

    def _equations_steady_state(self, t, y):
        dy = np.zeros_like(y)

        omega = y[:self.gas.n_species]
        z = y[self.gas.n_species:self.gas.n_species + self.surf.n_species]
        T = y[-1]

        self.gas.TPY = T, self.pressure, omega
        self.surf.TP = T, self.pressure
        self.surf.coverages = z

        gas_reaction_rates = self.gas.net_production_rates
        gas_reaction_rates_from_surface = self.surf.net_production_rates[:self.gas.n_species]
        coverage_reaction_rates = self.surf.net_production_rates[-self.surf.n_species:]

        dy[:self.gas.n_species] = (self.gas.molecular_weights * gas_reaction_rates + self.alfa * gas_reaction_rates_from_surface * self.gas.molecular_weights) * self.area / self.m_in

        dy[self.gas.n_species:self.gas.n_species + self.surf.n_species] = 1e03 * (coverage_reaction_rates / self.surf.site_density)

        if self.energy:
            if self.gas.n_reactions > 0:
                heat_of_reaction_from_gas = -np.dot(self.gas.net_rates_of_progress, self.gas.delta_enthalpy)
            else:
                heat_of_reaction_from_gas = 0.

            if self.surf.n_reactions > 0:
                heat_from_reaction_from_surface = -np.dot(self.surf.net_rates_of_progress, self.surf.delta_enthalpy)
            else:
                heat_from_reaction_from_surface = 0.

            dy[-1] = (heat_of_reaction_from_gas + self.alfa * heat_from_reaction_from_surface) * self.area / (self.m_in * self.gas.cp_mass)

        return dy

    def _solve_steady_state(self):
        self.length, self.sol = self._solve_ode(self._equations_steady_state, self._initial_conditions_steady_state(), self.length, atol=self.atol, rtol=self.rtol)

        self.y_sol = self.sol[:, :self.gas.n_species]
        self.x_sol = np.zeros_like(self.y_sol)

        for i in range(0, len(self.length)):
            self.x_sol[i, :self.gas.n_species] = self._convert_mass_fraction_to_mole_fraction(self.y_sol[i, :self.gas.n_species])

        self.coverage_sol = self.sol[:, self.gas.n_species:self.gas.n_species + self.surf.n_species]
        self.temperature_sol = self.sol[:, -1]

        self.is_solved = True
        return self.sol

    def _initial_conditions_steady_state(self):
        if not self.is_mass_flow_rate:
            self.gas.TPY = self.T_in, self.pressure, self.inlet_mass_fraction
            self.m_in = self.inlet_volumetric_flow_rate * self.gas.density
        return np.block([self.inlet_mass_fraction, self.surf.coverages, self.T_in])

    def _equations_transient(self, t, y):
        NP = self.length.size
        NV = self.gas.n_species + self.surf.n_species + 1

        y_matrix = y.reshape(NP, NV)

        omega = y_matrix[:, :self.gas.n_species]
        z = y_matrix[:, self.gas.n_species:self.gas.n_species + self.surf.n_species]
        T = y_matrix[:, -1]

        gas_reaction_rates = np.zeros([NP, self.gas.n_species], dtype=np.float64)
        gas_reaction_rates_from_surface = np.zeros([NP, self.gas.n_species], dtype=np.float64)
        coverage_reaction_rates = np.zeros([NP, self.surf.n_species], dtype=np.float64)
        mix_diff_coeffs_mass = np.zeros([NP, self.surf.n_species], dtype=np.float64)
        heat_of_reaction_from_gas = np.zeros([NP], dtype=np.float64)
        heat_from_reaction_from_surface = np.zeros([NP], dtype=np.float64)
        density = np.zeros([NP], dtype=np.float64)
        cp_mass = np.zeros([NP], dtype=np.float64)
        thermal_conductivity = np.zeros([NP], dtype=np.float64)

        for i in range(0, NP):
            self.gas.TPY = T[i], self.pressure, omega[i, :]
            self.surf.TP = T[i], self.pressure
            self.surf.coverages = z[i, :]

            gas_reaction_rates[i, :] = self.gas.net_production_rates * self.gas.molecular_weights
            gas_reaction_rates_from_surface[i, :] = self.surf.net_production_rates[:self.gas.n_species] * self.gas.molecular_weights
            coverage_reaction_rates[i, :] = self.surf.net_production_rates[-self.surf.n_species:]

            density[i] = self.gas.density
            cp_mass[i] = self.gas.cp_mass
            thermal_conductivity[i] = self.gas.thermal_conductivity / (self.gas.density * self.gas.cp_mass)

            diff_mix = self.gas.mix_diff_coeffs_mass
            diff_mix_zero = diff_mix == 0
            diff_mix[diff_mix_zero] = self.gas.binary_diff_coeffs[diff_mix_zero, diff_mix_zero]
            mix_diff_coeffs_mass[i, :] = diff_mix

            if self.energy:
                if self.gas.n_reactions > 0:
                    heat_of_reaction_from_gas[i] = -np.dot(self.gas.net_rates_of_progress, self.gas.delta_enthalpy)

                if self.surf.n_reactions > 0:
                    heat_from_reaction_from_surface[i] = -np.dot(self.surf.net_rates_of_progress, self.surf.delta_enthalpy)

        domega = np.zeros_like(omega)
        dT = np.zeros_like(T)

        domega[0, :] = self.inlet_mass_fraction - omega[0, :]  # Outlet conditions
        domega[-1, :] = omega[-1, :] - omega[-2, :]  # Outlet conditions

        dT[0] = self.T_in - T[0]  # Outlet conditions
        dT[-1] = T[-1] - T[-2]  # Outlet conditions

        derivative_1st_omega_forward = ((omega[2:, :] - omega[1:-1, :]).T / (self.length[2:] - self.length[1:-1])).T
        derivative_1st_omega_backward = ((omega[1:-1, :] - omega[:-2, :]).T / (self.length[1:-1] - self.length[:-2])).T

        diffusion_coefficient_forward = 0.5 * (mix_diff_coeffs_mass[2:, :] + mix_diff_coeffs_mass[1:-1, :])
        diffusion_coefficient_backward = 0.5 * (mix_diff_coeffs_mass[1:-1, :] + mix_diff_coeffs_mass[:-2, :])

        domega[1:-1, :] = - (derivative_1st_omega_backward.T * (self.m_in / (self.area * density[1:-1]))).T
        domega[1:-1, :] = domega[1:-1, :] + ((diffusion_coefficient_forward * derivative_1st_omega_forward - diffusion_coefficient_backward * derivative_1st_omega_backward).T / (0.5 * (self.length[2:] - self.length[:-2]))).T
        domega[1:-1, :] = domega[1:-1, :] + (gas_reaction_rates[1:-1, :].T / density[1:-1]).T
        domega[1:-1, :] = domega[1:-1, :] + self.alfa * (gas_reaction_rates_from_surface[1:-1, :].T / density[1:-1]).T

        dz = coverage_reaction_rates / self.surf.density

        if self.energy:
            derivative_1st_temperature_forward = (T[2:] - T[1:-1]) / (self.length[2:] - self.length[1:-1])
            derivative_1st_temperature_backward = (T[1:-1] - T[:-2]) / (self.length[1:-1] - self.length[:-2])

            thermal_coefficient_forward = 0.5 * (thermal_conductivity[2:] + thermal_conductivity[1:-1])
            thermal_coefficient_backward = 0.5 * (thermal_conductivity[1:-1] + thermal_conductivity[:-2])

            dT[1:-1] = -(self.m_in / (self.area * density[1:-1])) * derivative_1st_temperature_backward
            dT[1:-1] = dT[1:-1] + (thermal_coefficient_forward * derivative_1st_temperature_forward - thermal_coefficient_backward * derivative_1st_temperature_backward) / (0.5 * (self.length[2:] - self.length[:-2]))
            dT[1:-1] = dT[1:-1] + heat_of_reaction_from_gas[1:-1] / (density[1:-1] * cp_mass[1:-1])
            dT[1:-1] = dT[1:-1] + self.alfa * heat_from_reaction_from_surface[1:-1] / (density[1:-1] * cp_mass[1:-1])

        dy_matrix = np.zeros_like(y_matrix)

        dy_matrix[:, :self.gas.n_species] = domega
        dy_matrix[:, self.gas.n_species:self.gas.n_species + self.surf.n_species] = dz
        dy_matrix[:, -1] = dT

        return dy_matrix.flatten()

    def _ode_equations_for_transient(self, t, y):
        return self._equations_transient(t, y) * np.fabs(np.round(self.alg - 1))

    def _algebraic_equations(self):
        NP = self.length.size
        NV = self.gas.n_species + self.surf.n_species + 1

        alg_matrix = np.ones([NP, NV], dtype=int)

        alg_matrix[0, :self.gas.n_species] = 0
        alg_matrix[0, -1] = 0
        alg_matrix[-1, :self.gas.n_species] = 0
        alg_matrix[-1, -1] = 0

        return alg_matrix.flatten()

    def _residuals_transient(self, t, y, dy):
        res = self._equations_transient(t, y)
        diff_mask = self.alg == 1
        res[diff_mask] = res[diff_mask] - dy[diff_mask]
        return res

    def _solve_transient(self, tspan, time_ud):
        self.alg = self._algebraic_equations()
        self.tspan, self.sol = self._solve_dae(self._ode_equations_for_transient,
                                               self._equations_transient,
                                               self._residuals_transient,
                                               self._initial_conditions_transient(),
                                               self.uc.convert_to_seconds(tspan, time_ud),
                                               self.alg)

        self.y_sol = np.zeros([self.tspan.size], dtype=np.ndarray)
        self.x_sol = np.zeros([self.tspan.size], dtype=np.ndarray)
        self.coverage_sol = np.zeros([self.tspan.size], dtype=np.ndarray)
        self.temperature_sol = np.zeros([self.tspan.size], dtype=np.ndarray)

        NP = self.length.size
        NV = self.gas.n_species + self.surf.n_species + 1

        for i in range(0, self.sol.shape[0]):
            sol_for_time = self.sol[i, :].reshape(NP, NV)
            self.y_sol[i] = sol_for_time[:, :self.gas.n_species]
            self.coverage_sol[i] = sol_for_time[:, self.gas.n_species:self.gas.n_species + self.surf.n_species]
            self.temperature_sol[i] = sol_for_time[:, -1]
            self.x_sol[i] = np.zeros_like(self.y_sol[i])
            for j in range(0, NP):
                self.x_sol[i][j, :] = self._convert_mass_fraction_to_mole_fraction(self.y_sol[i][j, :])

        self.is_solved = True
        return self.sol

    def _initial_conditions_transient(self):
        NP = self.length.size
        NV = self.gas.n_species + self.surf.n_species + 1

        y0_matrix = np.zeros([NP, NV], dtype=np.float64)

        y0_matrix[:, :self.gas.n_species] = self.gas.Y
        y0_matrix[:, self.gas.n_species:self.gas.n_species + self.surf.n_species] = self.surf.coverages
        y0_matrix[:, -1] = self.temperature

        if not self.energy:
            self.T_in = self.temperature

        y0_matrix[0, :self.gas.n_species] = self.inlet_mass_fraction
        y0_matrix[0, -1] = self.T_in

        if not self.is_mass_flow_rate:
            self.gas.TPY = self.T_in, self.pressure, self.inlet_mass_fraction
            self.m_in = self.inlet_volumetric_flow_rate * self.gas.density

        return y0_matrix.flatten()

    def equations(self, t, y):
        if self.resolution_method == ResolutionMethod.STEADYSTATE:
            return self._equations_steady_state(t, y)

        return self._equations_transient(t, y)

    def initial_condition(self):
        if self.resolution_method == ResolutionMethod.STEADYSTATE:
            return self._initial_conditions_steady_state()

        return self._initial_conditions_transient()

    def solve(self, tspan=None, time_ud=None):
        if self.resolution_method == ResolutionMethod.STEADYSTATE:
            return self._solve_steady_state()

        return self._solve_transient(tspan, time_ud)
