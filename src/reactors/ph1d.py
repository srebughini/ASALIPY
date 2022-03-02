from asali.reactors.basic import ReactorType, BasicReactor, ResolutionMethod

import numpy as np


class PseudoHomogeneous1DReactor(BasicReactor):
    def __init__(self, cantera_input_file, gas_phase_name, surface_phase_name):
        super().__init__(cantera_input_file=cantera_input_file, gas_phase_name=gas_phase_name,
                         surface_phase_name=surface_phase_name)
        self.reactor_type = ReactorType.PSEUDOHOMOGENEOUSPFR
        self.alg = None
        self.gas_diffusion = False

    def _equations_steady_state(self, t, y):
        """
        Function representing the Steady State equations
        :param t: Independent variable - Reactor length
        :param y: Dependent variable - Species composition, coverage and temperature
        :return:
        """
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

        dy[:self.gas.n_species] = (
                                          self.gas.molecular_weights * gas_reaction_rates + self.alfa * gas_reaction_rates_from_surface * self.gas.molecular_weights) * self.area / self.inlet_mass_flow_rate

        dy[self.gas.n_species:self.gas.n_species + self.surf.n_species] = 1e03 * (
                coverage_reaction_rates / self.surf.site_density)

        if self.energy:
            if self.gas.n_reactions > 0:
                heat_of_reaction_from_gas = -np.dot(self.gas.net_rates_of_progress, self.gas.delta_enthalpy)
            else:
                heat_of_reaction_from_gas = 0.

            if self.surf.n_reactions > 0:
                heat_from_reaction_from_surface = -np.dot(self.surf.net_rates_of_progress, self.surf.delta_enthalpy)
            else:
                heat_from_reaction_from_surface = 0.

            dy[-1] = (heat_of_reaction_from_gas + self.alfa * heat_from_reaction_from_surface) * self.area / (
                    self.inlet_mass_flow_rate * self.gas.cp_mass)

        return dy

    def _solve_steady_state(self):
        """
        Solve Steady State equations
        :return: Matrix representing the solution in terms of composition, coverage and temperature as function of reactor length
        """
        self.length, self.sol = self._solve_ode(self._equations_steady_state,
                                                self._initial_conditions_steady_state(),
                                                self.length,
                                                atol=self.atol,
                                                rtol=self.rtol,
                                                verbosity=self.verbosity)

        self.y_sol = self.sol[:, :self.gas.n_species]
        self.x_sol = np.zeros_like(self.y_sol)

        for i in range(0, len(self.length)):
            self.x_sol[i, :self.gas.n_species] = self._convert_mass_fraction_to_mole_fraction(
                self.y_sol[i, :self.gas.n_species])

        self.coverage_sol = self.sol[:, self.gas.n_species:self.gas.n_species + self.surf.n_species]
        self.temperature_sol = self.sol[:, -1]

        self.is_solved = True
        return self.sol

    def _initial_conditions_steady_state(self):
        """
        Function creating the initial condition of the Steady State solution
        :return: Matrix representing the initial mass fraction, coverage and temperature
        """
        if not self.is_mass_flow_rate:
            self.gas.TPY = self.inlet_temperature, self.pressure, self.inlet_mass_fraction
            self.inlet_mass_flow_rate = self.inlet_volumetric_flow_rate * self.gas.density
        return np.block([self.inlet_mass_fraction, self.initial_coverage, self.inlet_temperature])

    def _inlet_conditions(self, omega, T):
        """
        Function estimating the inlet condition for the Transient solution
        :param omega: Matrix representing the species mass fraction
        :param T: Vector representing the gas temperature
        :return: Inlet conditions for species composition and temperature
        """
        return self.inlet_mass_fraction - omega[0, :], self.inlet_temperature - T[0]

    def _outlet_conditions(self, omega, T, density):
        """
        Function estimating the outlet condition for the Transient solution
        :param omega: Matrix representing the species mass fraction
        :param T: Vector representing the gas temperature
        :param density: Vector representing the gas density
        :return: Outlet conditions for species composition and temperature
        """
        if self.gas_diffusion:
            if self.energy:
                return omega[-1, :] - omega[-2, :], T[-1] - T[-2]

            return omega[-1, :] - omega[-2, :], 0.

        derivative_1st_omega = (omega[-1, :] - omega[-2, :]) / (self.length[-1] - self.length[-2])
        omega_outlet = - self.inlet_mass_flow_rate * derivative_1st_omega / (self.area * density[-1])

        if self.energy:
            derivative_1st_temperature = (T[-1] - T[-2]) / (self.length[-1] - self.length[-2])
            temperature_outlet = - self.inlet_mass_flow_rate * derivative_1st_temperature / (self.area * density[-1])
            return omega_outlet, temperature_outlet

        return omega_outlet, 0.

    def _equations_transient(self, t, y):
        """
        Function representing the Transient equations
        :param t: Independent variable - Time
        :param y: Dependent variable - Species composition, coverage and temperature as function of reactor length
        :return:
        """
        NP = self.length.size
        NV = self.gas.n_species + self.surf.n_species + 1

        y_matrix = y.reshape(NP, NV)

        omega = y_matrix[:, :self.gas.n_species]
        z = y_matrix[:, self.gas.n_species:self.gas.n_species + self.surf.n_species]
        T = y_matrix[:, -1]

        gas_reaction_rates = np.zeros([NP, self.gas.n_species], dtype=np.float64)
        gas_reaction_rates_from_surface = np.zeros([NP, self.gas.n_species], dtype=np.float64)
        coverage_reaction_rates = np.zeros([NP, self.surf.n_species], dtype=np.float64)
        mix_diff_coeffs_mass = np.zeros([NP, self.gas.n_species], dtype=np.float64)
        heat_of_reaction_from_gas = np.zeros([NP], dtype=np.float64)
        heat_from_reaction_from_surface = np.zeros([NP], dtype=np.float64)
        density = np.zeros([NP], dtype=np.float64)
        cp_mass = np.zeros([NP], dtype=np.float64)
        thermal_conductivity = np.zeros([NP], dtype=np.float64)

        for i in range(0, NP):
            self.gas.TPY = T[i], self.pressure, np.maximum(0.0, omega[i, :])
            self.surf.TP = T[i], self.pressure
            self.surf.coverages = np.maximum(0.0, z[i, :])

            if self.gas.n_reactions > 0:
                gas_reaction_rates[i, :] = self.gas.net_production_rates * self.gas.molecular_weights

            if self.surf.n_reactions > 0:
                gas_reaction_rates_from_surface[i, :] = self.surf.net_production_rates[
                                                        :self.gas.n_species] * self.gas.molecular_weights
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
                    heat_from_reaction_from_surface[i] = -np.dot(self.surf.net_rates_of_progress,
                                                                 self.surf.delta_enthalpy)

        domega = np.zeros_like(omega)
        dT = np.zeros_like(T)

        # Inlet and outlet conditions
        domega[0, :], dT[0] = self._inlet_conditions(omega, T)
        domega[-1, :], dT[-1] = self._outlet_conditions(omega, T, density)

        # Equations for mass
        d1st_omega_backward = (omega[1:-1, :] - omega[:-2, :]) / (
                self.length[1:-1, np.newaxis] - self.length[:-2, np.newaxis])
        domega[1:-1, :] = d1st_omega_backward / density[1:-1, np.newaxis]
        domega[1:-1, :] = - (self.inlet_mass_flow_rate / self.area) * domega[1:-1, :]
        domega[1:-1, :] = domega[1:-1, :] + (gas_reaction_rates[1:-1, :] / density[1:-1, np.newaxis])
        domega[1:-1, :] = domega[1:-1, :] + self.alfa * (
                gas_reaction_rates_from_surface[1:-1, :] / density[1:-1, np.newaxis])

        if self.gas_diffusion:
            d1st_omega_forward = (omega[2:, :] - omega[1:-1, :]) / (
                    self.length[2:, np.newaxis] - self.length[1:-1, np.newaxis])
            diff_coeff_forward = 0.5 * (mix_diff_coeffs_mass[2:, :] + mix_diff_coeffs_mass[1:-1, :])
            diff_coeff_backward = 0.5 * (mix_diff_coeffs_mass[1:-1, :] + mix_diff_coeffs_mass[:-2, :])
            domega[1:-1, :] = domega[1:-1, :] + (
                    (diff_coeff_forward * d1st_omega_forward - diff_coeff_backward * d1st_omega_backward) / (
                    0.5 * (self.length[2:, np.newaxis] - self.length[:-2, np.newaxis])))

        # Inert specie
        domega[:, self.inert_specie_index] = 1. - np.sum(omega, axis=1)

        # Equations for site fraction
        dz = coverage_reaction_rates / self.surf.site_density

        # Equations for temperature
        if self.energy:
            d1st_temperature_backward = (T[1:-1] - T[:-2]) / (self.length[1:-1] - self.length[:-2])
            dT[1:-1] = -(self.inlet_mass_flow_rate / (self.area * density[1:-1])) * d1st_temperature_backward
            dT[1:-1] = dT[1:-1] + heat_of_reaction_from_gas[1:-1] / (density[1:-1] * cp_mass[1:-1])
            dT[1:-1] = dT[1:-1] + self.alfa * heat_from_reaction_from_surface[1:-1] / (density[1:-1] * cp_mass[1:-1])

            if self.gas_diffusion:
                d1st_temperature_forward = (T[2:] - T[1:-1]) / (self.length[2:] - self.length[1:-1])
                thermal_coeff_forward = 0.5 * (thermal_conductivity[2:] + thermal_conductivity[1:-1])
                thermal_coeff_backward = 0.5 * (thermal_conductivity[1:-1] + thermal_conductivity[:-2])
                dT[1:-1] = dT[1:-1] + (
                        thermal_coeff_forward * d1st_temperature_forward - thermal_coeff_backward * d1st_temperature_backward) / (
                                   0.5 * (self.length[2:] - self.length[:-2]))

        dy_matrix = np.zeros_like(y_matrix)

        dy_matrix[:, :self.gas.n_species] = domega
        dy_matrix[:, self.gas.n_species:self.gas.n_species + self.surf.n_species] = dz

        if self.energy:
            dy_matrix[:, -1] = dT

        return dy_matrix.flatten()

    def _ode_equations_for_transient(self, t, y):
        return self._equations_transient(t, y) * np.fabs(np.round(self.alg - 1))

    def _algebraic_equations(self):
        """
        Generate the vector describing algebraic (0) and differential (1) equations
        :return: Vector on 0/1 describing algebraic/differential equations
        """
        NP = self.length.size
        NV = self.gas.n_species + self.surf.n_species + 1

        alg_matrix = np.ones([NP, NV], dtype=int)

        alg_matrix[0, :self.gas.n_species] = 0
        if self.energy:
            alg_matrix[0, -1] = 0

        if self.gas_diffusion:
            alg_matrix[-1, :self.gas.n_species] = 0
            if self.energy:
                alg_matrix[-1, -1] = 0

        alg_matrix[:, self.inert_specie_index] = 0

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
                                               self.alg,
                                               atol=self.atol,
                                               rtol=self.rtol,
                                               verbosity=self.verbosity)

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
            self.inlet_temperature = self.temperature

        y0_matrix[0, :self.gas.n_species] = self.inlet_mass_fraction
        y0_matrix[0, -1] = self.inlet_temperature

        if not self.is_mass_flow_rate:
            self.gas.TPY = self.inlet_temperature, self.pressure, self.inlet_mass_fraction
            self.inlet_mass_flow_rate = self.inlet_volumetric_flow_rate * self.gas.density

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

    def set_gas_diffusion(self, gas_diffusion):
        self.gas_diffusion = self._true_parser(gas_diffusion)
