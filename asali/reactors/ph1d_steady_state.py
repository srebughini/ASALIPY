import numpy as np

from asali.reactors.ph1d import PseudoHomogeneous1DReactor
from asali.utils.input_parser import ResolutionMethod


class SteadyStatePseudoHomogeneous1DReactor(PseudoHomogeneous1DReactor):
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

    def equations(self, t, y):
        """
        Function representing the Steady State equations
        :param t: Independent variable - Reactor length
        :param y: Dependent variable - Species composition, coverage and temperature
        :return:
        """
        dy = np.zeros(shape=y.shape, dtype=np.float64)

        omega = y[:self.gas.n_species]
        z = y[self.gas.n_species:self.gas.n_species + self.surf.n_species]
        T = y[-1]

        self.gas.TPY = T, self.pressure, omega
        self.surf.TP = T, self.pressure
        self.surf.coverages = z

        r_gas = self.get_homogeneous_gas_species_reaction_rates()
        r_from_surface = self.get_heterogeneous_gas_species_reaction_rates()
        r_surface = self.get_surface_species_reaction_rates()

        # Equation of mass
        domega = self.gas.molecular_weights * r_gas
        domega = domega + self.alfa * r_from_surface * self.gas.molecular_weights
        domega = domega * self.area / self.inlet_mass_flow_rate

        # Equation of coverage
        dz = 1e03 * (r_surface / self.surf.site_density)

        dT = 0.0
        if self.energy:
            q_from_gas = self.get_homogeneous_heat_of_reaction()
            q_from_surface = self.get_heterogeneous_heat_of_reaction()

            # Equation of energy
            dT = (q_from_gas + self.alfa * q_from_surface) * self.area / (self.inlet_mass_flow_rate * self.gas.cp_mass)

        dy[:self.gas.n_species] = domega
        dy[self.gas.n_species:self.gas.n_species + self.surf.n_species] = dz
        dy[-1] = dT

        return dy

    def initial_condition(self):
        """
        Generate initial conditions
        :return: Vector/Matrix representing the initial conditions
        """
        if not self.is_mass_flow_rate:
            self.gas.TPY = self.inlet_temperature, self.pressure, self.inlet_mass_fraction
            self.inlet_mass_flow_rate = self.inlet_volumetric_flow_rate * self.gas.density
        return np.block([self.inlet_mass_fraction, self.initial_coverage, self.inlet_temperature])

    def solve(self, tspan=None, time_ud=None):
        """
        Solve selected model
        :param tspan: Not required
        :param time_ud: Not required
        :return: Vector/Matrix representing the results
        """
        x, y = self.numerical_solver.solve_ode(self.equations,
                                               self.initial_condition(),
                                               self.length)

        self.solution_parser.x = x
        self.solution_parser.y = y
        self.solution_parser.is_solved = True
        return y
