import numpy as np

from asali.reactors.het1d import Heterogeneous1DReactor
from asali.utils.input_parser import ResolutionMethod


class TransientHeterogeneous1DReactor(Heterogeneous1DReactor):
    def __init__(self, cantera_input_file, gas_phase_name, surface_phase_name):
        """
        Class representing Steady State pseudoHomogeneous 1D reactor model
        :param cantera_input_file: Cantera file path
        :param gas_phase_name: Cantera gas phase name
        :param surface_phase_name: Cantera interface phase name
        """
        super().__init__(cantera_input_file=cantera_input_file, gas_phase_name=gas_phase_name,
                         surface_phase_name=surface_phase_name)

        self.solution_parser.resolution_method = ResolutionMethod.TRANSIENT

        self.n_s = self.gas.n_species
        self.n_surf = self.surf.n_species
        self.n_v = self.n_s + self.n_s + self.n_surf + 1 + 1

    def initial_condition(self):
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

    def solve(self, tspan, time_ud):
        """
        Solve model
        :param tspan: Vector representing the integration time
        :param time_ud: Time unit dimension
        :return: Vector/Matrix representing the results
        :return: Time vector,
                 Matrix representing the solution in terms of composition, coverage and temperature as function of time and reactor length
        """
        self.alg = self.algebraic_equations()
        x, y = self.numerical_solver.solve_dae(self.ode_equations,
                                               self.equations,
                                               self.residuals,
                                               self.initial_condition(),
                                               self.uc.convert_to_seconds(tspan, time_ud),
                                               self.alg)

        self.solution_parser.x = x
        self.solution_parser.y = y
        self.solution_parser.length = self.length
        self.solution_parser.is_solved = True
        return y

