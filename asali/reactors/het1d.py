################################################################################################
#                                                                                              #
#     #############       #############       #############       ####                ####     #
#    #             #     #             #     #             #     #    #              #    #    #
#    #    #####    #     #    #########      #    #####    #     #    #              #    #    #
#    #    #   #    #     #    #              #    #   #    #     #    #              #    #    #
#    #    #####    #     #    #              #    #####    #     #    #              #    #    #
#    #             #     #    #########      #             #     #    #              #    #    #
#    #             #     #             #     #             #     #    #              #    #    #
#    #    #####    #      #########    #     #    #####    #     #    #              #    #    #
#    #    #   #    #              #    #     #    #   #    #     #    #              #    #    #
#    #    #   #    #      #########    #     #    #   #    #     #    #########      #    #    #
#    #    #   #    #     #             #     #    #   #    #     #             #     #    #    #
#     ####     ####       #############       ####     ####       #############       ####     #
#                                                                                              #
#   Author: Stefano Rebughini <ste.rebu@outlook.it>                                            #
#                                                                                              #
################################################################################################
#                                                                                              #
#   License                                                                                    #
#                                                                                              #
#   This file is part of ASALI.                                                                #
#                                                                                              #
#   ASALI is free software: you can redistribute it and/or modify                              #
#   it under the terms of the GNU General Public License as published by                       #
#   the Free Software Foundation, either version 3 of the License, or                          #
#   (at your option) any later version.                                                        #
#                                                                                              #
#   ASALI is distributed in the hope that it will be useful,                                   #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of                             #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                              #
#   GNU General Public License for more details.                                               #
#                                                                                              #
#   You should have received a copy of the GNU General Public License                          #
#   along with ASALI. If not, see <http://www.gnu.org/licenses/>.                              #
#                                                                                              #
################################################################################################
from asali.reactors.basic import ReactorType, BasicReactor, ReactorSection
from enum import IntEnum

import numpy as np


class ReactorModel(IntEnum):
    TUBULAR = 0
    PACKEDBED = 1
    HONEYCOMB = 2


class Heterogeneous1DReactor(BasicReactor):
    def __init__(self, cantera_input_file, gas_phase_name, surface_phase_name):
        super().__init__(cantera_input_file=cantera_input_file, gas_phase_name=gas_phase_name, surface_phase_name=surface_phase_name)
        self.reactor_type = ReactorType.HETEROGENEOUSPRF
        self.reactor_model = None
        self.alg = None

        self.Nusselt = 0.
        self.Sherwood = 0.
        self.void_fraction = 0.
        self.tube_diameter = 0.
        self.specific_area = 0.
        self.particle_diameter = 0.

        self.y_sol_wall = None
        self.x_sol_wall = None
        self.temperature_sol_wall = None

    def _estimate_mass_transfer_coefficient(self, viscosity, density, diffusivity):
        if self.reactor_model == ReactorModel.PACKEDBED:
            reynolds_real = self.m_in * self.particle_diameter / (viscosity * 0.25 * np.square(self.tube_diameter))
            reynolds = reynolds_real / ((1. - self.void_fraction) * 6)

            laminar_flow = np.logical_and(reynolds < 50, reynolds > 0)
            turbulent_flow = reynolds >= 50

            j_mass = np.zeros_like(reynolds_real)
            j_mass[laminar_flow] = 0.91 / (np.power(reynolds[laminar_flow], 0.51))
            j_mass[turbulent_flow] = 0.61 / (np.power(reynolds[turbulent_flow], 0.41))

            Sc = ((viscosity / density) / diffusivity.T).T

            return ((j_mass * reynolds_real) * (np.power(Sc, 1. / 3.) * diffusivity).T).T / self.particle_diameter

        reynolds_real = self.m_in * self.tube_diameter / (viscosity * self.void_fraction * 0.25 * np.square(self.tube_diameter))
        Sc = ((viscosity / density) / diffusivity.T).T

        z_star = np.fabs(np.maximum(1e-06, self.length)) / (self.tube_diameter * reynolds_real)

        z_star = (z_star / Sc.T).T

        return ((self.Sherwood + 6.874 * np.power(1000 * z_star, -0.488) * np.exp(-57.2 * z_star)) * diffusivity) / self.tube_diameter

    def _estimate_heat_transfer_coefficient(self, viscosity, conductivity, specific_heat):
        if self.reactor_model == ReactorModel.PACKEDBED:
            reynolds_real = self.m_in * self.particle_diameter / (viscosity * 0.25 * np.square(self.tube_diameter))
            reynolds = reynolds_real / ((1. - self.void_fraction) * 6)

            laminar_flow = np.logical_and(reynolds < 50, reynolds > 0)
            turbulent_flow = reynolds >= 50

            j_heat = np.zeros_like(reynolds_real)
            j_heat[laminar_flow] = 0.91 / (np.power(reynolds[laminar_flow], 0.51))
            j_heat[turbulent_flow] = 0.61 / (np.power(reynolds[turbulent_flow], 0.41))

            Pr = (specific_heat * viscosity) / conductivity

            return j_heat * reynolds_real * np.power(Pr, 1. / 3.) * conductivity / self.particle_diameter

        reynolds_real = self.m_in * self.tube_diameter / (viscosity * self.void_fraction * 0.25 * np.square(self.tube_diameter))
        Pr = specific_heat * viscosity / conductivity
        z_star = np.fabs(np.maximum(1e-06, self.length)) / (self.tube_diameter * reynolds_real * Pr)
        return ((self.Nusselt + 8.827 * np.power(1000 * z_star, -0.545) * np.exp(-48.2 * z_star)) * conductivity) / self.tube_diameter

    def _ode_equations(self, t, y):
        return self.equations(t, y) * np.fabs(np.round(self.alg - 1))

    def _algebraic_equations(self):
        NP = self.length.size
        NV = self.gas.n_species + self.gas.n_species + self.surf.n_species + 1 + 1

        alg_matrix = np.ones([NP, NV], dtype=int)

        alg_matrix[0, :self.gas.n_species] = 0

        alg_matrix[:, self.gas.n_species:self.gas.n_species + self.gas.n_species] = 0

        alg_matrix[0, -2] = 0
        alg_matrix[0, -1] = 0
        alg_matrix[-1, -2] = 0
        alg_matrix[-1, -1] = 0

        return alg_matrix.flatten()

    def _residuals(self, t, y, dy):
        res = self.equations(t, y)
        diff_mask = self.alg == 1
        res[diff_mask] = res[diff_mask] - dy[diff_mask]
        return res

    def equations(self, t, y):
        NP = self.length.size
        NV = self.gas.n_species + self.gas.n_species + self.surf.n_species + 1 + 1

        y_matrix = y.reshape(NP, NV)

        omega_bulk = y_matrix[:, :self.gas.n_species]
        omega_wall = y_matrix[:, self.gas.n_species:self.gas.n_species + self.gas.n_species]
        z = y_matrix[:, self.gas.n_species + self.gas.n_species:self.gas.n_species + self.gas.n_species + self.surf.n_species]
        T_bulk = y_matrix[:, -2]
        T_wall = y_matrix[:, -1]

        gas_reaction_rates = np.zeros([NP, self.gas.n_species], dtype=np.float64)
        gas_reaction_rates_from_surface = np.zeros([NP, self.gas.n_species], dtype=np.float64)
        coverage_reaction_rates = np.zeros([NP, self.surf.n_species], dtype=np.float64)
        mix_diff_coeffs_mass = np.zeros([NP, self.surf.n_species], dtype=np.float64)
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
            gas_reaction_rates_from_surface[i, :] = self.surf.net_production_rates[:self.gas.n_species] * self.gas.molecular_weights
            coverage_reaction_rates[i, :] = self.surf.net_production_rates[-self.surf.n_species:]

            if self.energy:
                if self.surf.n_reactions > 0:
                    heat_from_reaction_from_surface[i] = -np.dot(self.surf.net_rates_of_progress, self.surf.delta_enthalpy)

        k_mat = self._estimate_mass_transfer_coefficient(viscosity, density, mix_diff_coeffs_mass)
        k_heat = self._estimate_heat_transfer_coefficient(viscosity, thermal_conductivity, cp_mass)

        domega_bulk = np.zeros_like(omega_bulk)
        dT_bulk = np.zeros_like(T_bulk)
        dT_wall = np.zeros_like(T_wall)

        domega_bulk[0, :] = self.y_in - omega_bulk[0, :]  # Inlet conditions
        dT_bulk[0] = self.T_in - T_bulk[0]  # Inlet conditions
        dT_bulk[-1] = T_bulk[-1] - T_bulk[-2]  # Outlet conditions

        dT_wall[0] = T_wall[1] - T_wall[0]  # Inlet conditions
        dT_wall[-1] = T_wall[-1] - T_wall[-2]  # Outlet conditions

        delta_omega = self.specific_area * (density * (k_mat * (omega_bulk - omega_wall)).T).T

        derivative_1st_omega_backward = ((omega_bulk[1:, :] - omega_bulk[:-1, :]).T / (self.length[1:] - self.length[:-1])).T

        domega_bulk[1:, :] = - (derivative_1st_omega_backward.T * (self.m_in / (self.area * density[1:]))).T
        domega_bulk[1:, :] = domega_bulk[1:, :] + (gas_reaction_rates[1:, :].T / density[1:]).T
        domega_bulk[1:, :] = domega_bulk[1:, :] - delta_omega[1:, :] / self.void_fraction

        domega_wall = delta_omega * self.void_fraction + self.alfa * self.void_fraction * gas_reaction_rates_from_surface

        dz = coverage_reaction_rates / self.surf.density

        if self.energy:
            derivative_1st_temperature_bulk_forward = (T_bulk[2:] - T_bulk[1:-1]) / (self.length[2:] - self.length[1:-1])
            derivative_1st_temperature_bulk_backward = (T_bulk[1:-1] - T_bulk[:-2]) / (self.length[1:-1] - self.length[:-2])

            thermal_coefficient_forward = 0.5 * (thermal_conductivity[2:] + thermal_conductivity[1:-1])
            thermal_coefficient_backward = 0.5 * (thermal_conductivity[1:-1] + thermal_conductivity[:-2])

            delta_T = k_heat * self.specific_area * (T_bulk - T_wall)

            dT_bulk[1:-1] = -(self.m_in / (self.area * density[1:-1])) * derivative_1st_temperature_bulk_backward
            dT_bulk[1:-1] = dT_bulk[1:-1] + (thermal_coefficient_forward * derivative_1st_temperature_bulk_forward - thermal_coefficient_backward * derivative_1st_temperature_bulk_backward) / (
                    0.5 * (self.length[2:] - self.length[:-2]) * density[1:-1] * cp_mass[1:-1])
            dT_bulk[1:-1] = dT_bulk[1:-1] + heat_of_reaction_from_gas[1:-1] / (density[1:-1] * cp_mass[1:-1])
            dT_bulk[1:-1] = dT_bulk[1:-1] - delta_T[1:-1] / (self.void_fraction * density[1:-1] * cp_mass[1:-1])

            derivative_1st_temperature_wall_forward = (T_wall[2:] - T_wall[1:-1]) / (self.length[2:] - self.length[1:-1])
            derivative_1st_temperature_wall_backward = (T_wall[1:-1] - T_wall[:-2]) / (self.length[1:-1] - self.length[:-2])

            dT_wall[1:-1] = (self.solid_k / (self.solid_cp * self.solid_rho)) * (derivative_1st_temperature_wall_forward - derivative_1st_temperature_wall_backward) / (0.5 * (self.length[2:] - self.length[:-2]))
            dT_wall[1:-1] = dT_wall[1:-1] + self.alfa * heat_from_reaction_from_surface[1:-1] / (self.solid_cp * self.solid_rho * (1 - self.void_fraction))
            dT_wall[1:-1] = dT_wall[1:-1] + delta_T[1:-1] / (self.solid_cp * self.solid_rho * (1 - self.void_fraction))

        dy_matrix = np.zeros_like(y_matrix)

        dy_matrix[:, :self.gas.n_species] = domega_bulk

        dy_matrix[:, self.gas.n_species:self.gas.n_species + self.surf.n_species] = domega_wall
        dy_matrix[:, self.gas.n_species + self.surf.n_species:self.gas.n_species + self.surf.n_species + self.surf.n_species] = dz
        dy_matrix[:, -2] = dT_bulk
        dy_matrix[:, -1] = dT_wall

        return dy_matrix.flatten()

    def initial_condition(self):
        NP = self.length.size
        NV = self.gas.n_species + self.gas.n_species + self.surf.n_species + 1 + 1

        y0_matrix = np.zeros([NP, NV], dtype=np.float64)

        y0_matrix[:, :self.gas.n_species] = self.gas.Y
        y0_matrix[:, self.gas.n_species:self.gas.n_species + self.gas.n_species] = self.gas.Y
        y0_matrix[:, self.gas.n_species + self.gas.n_species:self.gas.n_species + self.gas.n_species + self.surf.n_species] = self.surf.coverages
        y0_matrix[:, -2] = self.temperature
        y0_matrix[:, -1] = self.solid_temperature

        if not self.energy:
            self.T_in = self.temperature
            y0_matrix[:, -1] = self.temperature

        y0_matrix[0, :self.gas.n_species] = self.y_in
        y0_matrix[0, -1] = self.T_in

        if not self.is_mass_flow_rate:
            self.gas.TPY = self.T_in, self.pressure, self.y_in
            self.m_in = self.Q_in * self.gas.density

        return y0_matrix.flatten()

    def solve(self, tspan, time_ud):
        self.alg = self._algebraic_equations()
        self.tspan, self.sol = self._solve_dae(self._ode_equations,
                                               self.equations,
                                               self._residuals,
                                               self.initial_condition(),
                                               self.uc.convert_to_seconds(tspan, time_ud),
                                               self.alg)

        sol = np.copy(self.sol)
        self.y_sol = np.zeros([self.tspan.size], dtype=np.ndarray)
        self.x_sol = np.zeros([self.tspan.size], dtype=np.ndarray)
        self.y_sol_wall = np.zeros([self.tspan.size], dtype=np.ndarray)
        self.x_sol_wall = np.zeros([self.tspan.size], dtype=np.ndarray)
        self.coverage_sol = np.zeros([self.tspan.size], dtype=np.ndarray)
        self.temperature_sol = np.zeros([self.tspan.size], dtype=np.ndarray)
        self.temperature_sol_wall = np.zeros([self.tspan.size], dtype=np.ndarray)

        NP = self.length.size
        NV = self.gas.n_species + self.gas.n_species + self.surf.n_species + 1 + 1

        for i in range(0, self.sol.shape[0]):
            sol_for_time = sol[i, :].reshape(NP, NV)
            self.y_sol[i] = sol_for_time[:, :self.gas.n_species]
            self.y_sol_wall[i] = sol_for_time[:, self.gas.n_species:self.gas.n_species + self.gas.n_species]
            self.coverage_sol[i] = sol_for_time[:, self.gas.n_species + self.gas.n_species:self.gas.n_species + self.gas.n_species + self.surf.n_species]
            self.temperature_sol_wall[i] = sol_for_time[:, -2]
            self.temperature_sol[i] = sol_for_time[:, -1]
            self.x_sol[i] = np.zeros_like(self.y_sol[i])
            self.x_sol_wall[i] = np.zeros_like(self.y_sol_wall[i])
            for j in range(0, NP):
                self.x_sol[i][j, :] = self._convert_mass_fraction_to_mole_fraction(self.y_sol[i][j, :])
                self.x_sol_wall[i][j, :] = self._convert_mass_fraction_to_mole_fraction(self.y_sol_wall[i][j, :])

        self.is_solved = True
        return self.sol

    def set_tubular_reactor(self, tube_diameter, tube_diameter_ud, wall_thickness, wall_thickness_ud):
        Dt = self.uc.convert_to_meter(tube_diameter, tube_diameter_ud)
        tw = self.uc.convert_to_meter(wall_thickness, wall_thickness_ud)
        self.tube_diameter = Dt - 2. * tw
        self.specific_area = 4. / self.tube_diameter
        if self.reactor_section == ReactorSection.CIRCLE:
            self.Nusselt = 3.659
            self.Sherwood = 3.659
        elif self.reactor_section == ReactorSection.SQUARE:
            self.Nusselt = 2.977
            self.Sherwood = 2.977
        elif self.reactor_section == ReactorSection.TRIANGLE:
            self.Nusselt = 2.494
            self.Sherwood = 2.494

        self.void_fraction = np.square(self.tube_diameter / Dt)
        self.area = np.pi * 0.25 * np.square(self.tube_diameter)
        self.reactor_model = ReactorModel.TUBULAR

    def set_honeycomb_reactor(self, cpsi, wall_thickness, wall_thickness_ud):
        tw = self.uc.convert_to_meter(wall_thickness, wall_thickness_ud)

        if self.reactor_section == ReactorSection.CIRCLE:
            self.Nusselt = 3.659
            self.Sherwood = 3.659
        elif self.reactor_section == ReactorSection.SQUARE:
            self.Nusselt = 2.977
            self.Sherwood = 2.977
        elif self.reactor_section == ReactorSection.TRIANGLE:
            self.Nusselt = 2.494
            self.Sherwood = 2.494

        self.tube_diameter = np.sqrt(1. / cpsi) * 2.54 * 1e-02 - tw
        self.void_fraction = np.square(self.tube_diameter) / np.square(self.tube_diameter + tw)
        self.specific_area = 4. * self.void_fraction / self.tube_diameter
        self.area = np.pi * 0.25 * np.square(self.tube_diameter)
        self.reactor_model = ReactorModel.HONEYCOMB
        return [self.tube_diameter, self.void_fraction, self.particle_diameter, self.specific_area, self.area]

    def set_packed_bed_reactor(self, particle_diameter, particle_diameter_ud, tube_diameter, tube_diameter_ud, void_fraction):
        self.tube_diameter = self.uc.convert_to_meter(tube_diameter, tube_diameter_ud)
        self.void_fraction = void_fraction
        self.particle_diameter = self.uc.convert_to_meter(particle_diameter, particle_diameter_ud)
        self.specific_area = 6. * (1. - self.void_fraction) / self.particle_diameter
        self.area = np.pi * 0.25 * np.square(self.tube_diameter)
        self.reactor_model = ReactorModel.PACKEDBED
        return [self.tube_diameter, self.void_fraction, self.particle_diameter, self.specific_area, self.area]

    def get_solid_mass_fraction(self, index=None):
        if index is None:
            return self.y_sol_wall

        return self.y_sol_wall[index]

    def get_solid_mole_fraction(self, index=None):
        if index is None:
            return self.x_sol_wall

        return self.x_sol_wall[index]

    def get_solid_temperature(self, index=None):
        if index is None:
            return self.temperature_sol_wall

        return self.temperature_sol_wall[index]
