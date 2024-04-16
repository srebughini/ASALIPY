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

    def estimate_integration_time(self):
        """
        Estimate integration time
        :return: Integration time
        """
        self.gas.TPY = self.inlet_temperature, self.pressure, self.inlet_mass_fraction
        if self.is_mass_flow_rate:
            return 100000.0 * self.inlet_mass_flow_rate / (
                    self.gas.density * self.reactor_shape_object.section_area * self.reactor_shape_object.void_fraction)

        self.inlet_mass_flow_rate = self.inlet_volumetric_flow_rate * self.gas.density
        return 100000.0 * self.inlet_mass_flow_rate / (
                self.gas.density * self.reactor_shape_object.section_area * self.reactor_shape_object.void_fraction)

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
        tspan = [0, self.estimate_integration_time()]

        self.alg = self.algebraic_equations()

        _, y = self.numerical_solver.solve_dae(self.ode_equations,
                                               self.equations,
                                               self.residuals,
                                               self.initial_condition(),
                                               tspan,
                                               self.alg)

        self.solution_parser.x = self.length
        self.solution_parser.y = y[-1, :].reshape(self.n_p, self.n_v)
        self.solution_parser.length = self.length
        self.solution_parser.is_solved = True
        return self.solution_parser.y
