from asali.reactors.basic import BasicReactor
import numpy as np

from asali.utils.input_parser import ReactorType


class BatchReactor(BasicReactor):
    def __init__(self, cantera_input_file, gas_phase_name, surface_phase_name):
        """
        Class representing BATCH reactor model
        :param cantera_input_file: Cantera file path
        :param gas_phase_name: Cantera gas phase name
        :param surface_phase_name: Cantera interface phase name
        """
        super().__init__(cantera_input_file=cantera_input_file, gas_phase_name=gas_phase_name,
                         surface_phase_name=surface_phase_name)
        self.solution_parser.reactor_type = ReactorType.BATCH
        self.volume = 0.

    def set_volume(self, value, unit_dimension):
        """
        Set volume
        :param value: Volume value
        :param unit_dimension: Volume unit dimension
        :return: Volume in [m3]
        """
        self.volume = self.uc.convert_to_cubic_meter(value, unit_dimension)
        return self.volume

    def equations(self, t, y):
        """
        Function representing the model
        :param t: Independent variable - Time
        :param y: Dependent variable - Species composition, coverage and temperature
        :return:
        """
        dy = np.zeros(shape=y.shape, dtype=np.float64)

        omega = y[:self.gas.n_species]
        mass = y[self.gas.n_species]
        z = y[self.gas.n_species + 1:self.gas.n_species + 1 + self.surf.n_species]
        T = y[-1]

        self.gas.TPY = T, self.pressure, omega
        self.surf.TP = T, self.pressure
        self.surf.coverages = z

        gas_reaction_rates = self.gas.net_production_rates
        gas_reaction_rates_from_surface = self.surf.net_production_rates[:self.gas.n_species]
        coverage_reaction_rates = self.surf.net_production_rates[-self.surf.n_species:]

        dmass = self.volume * self.alfa * np.dot(gas_reaction_rates_from_surface, self.gas.molecular_weights)

        domega = self.gas.molecular_weights * gas_reaction_rates / self.gas.density #+ (
                #-omega * dy[self.gas.n_species] + (
                #mass / self.gas.density) * self.alfa * gas_reaction_rates_from_surface * self.gas.molecular_weights) / mass

        domega = domega - omega * dy[self.gas.n_species]/mass
        domega = domega + self.alfa * gas_reaction_rates_from_surface * self.gas.molecular_weights/ self.gas.density

        dz = coverage_reaction_rates / self.surf.site_density

        dT = 0
        if self.energy:
            if self.gas.n_reactions > 0:
                heat_of_reaction_from_gas = -np.dot(self.gas.net_rates_of_progress, self.gas.delta_enthalpy)
            else:
                heat_of_reaction_from_gas = 0.

            if self.surf.n_reactions > 0:
                heat_from_reaction_from_surface = -np.dot(self.surf.net_rates_of_progress, self.surf.delta_enthalpy)
            else:
                heat_from_reaction_from_surface = 0.

            dT = (heat_of_reaction_from_gas + self.alfa * heat_from_reaction_from_surface) / (
                    self.gas.density * self.gas.cp_mass)

        dy[:self.gas.n_species] = domega
        dy[self.gas.n_species] = dmass
        dy[self.gas.n_species + 1:self.gas.n_species + 1 + self.surf.n_species] = dz
        dy[-1] = dT

        return dy

    def initial_condition(self):
        """
        Generate initial conditions
        :return: Vector/Matrix representing the initial conditions
        """
        return np.block([self.initial_mass_fraction,
                         self.gas.density * self.volume,
                         self.initial_coverage,
                         self.initial_temperature])

    def solve(self, tspan, time_ud):
        """
        Solve model
        :param tspan: Vector representing the integration time
        :param time_ud: Time unit dimension
        :return: Vector/Matrix representing the results
        """
        x, y = self.numerical_solver.solve_ode(self.equations,
                                               self.initial_condition(),
                                               self.uc.convert_to_seconds(tspan, time_ud))

        self.solution_parser.x = x
        self.solution_parser.y = y
        self.solution_parser.is_solved = True
        return y
