from asali.reactors.basic import ReactorType, BasicReactor
import numpy as np


class BatchReactor(BasicReactor):
    def __init__(self, cantera_input_file, gas_phase_name, surface_phase_name):
        super().__init__(cantera_input_file=cantera_input_file, gas_phase_name=gas_phase_name, surface_phase_name=surface_phase_name)
        self.mass_sol = None
        self.reactor_type = ReactorType.BATCH
        self.initial_mass = None

    def equations(self, t, y):
        dy = np.zeros_like(y)

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

        dy[self.gas.n_species] = self.volume * self.alfa * np.dot(gas_reaction_rates_from_surface, self.gas.molecular_weights)

        dy[:self.gas.n_species] = self.gas.molecular_weights * gas_reaction_rates / self.gas.density + (
                -omega * dy[self.gas.n_species] + (mass / self.gas.density) * self.alfa * gas_reaction_rates_from_surface * self.gas.molecular_weights) / mass

        dy[self.gas.n_species + 1:self.gas.n_species + 1 + self.surf.n_species] = coverage_reaction_rates / self.surf.site_density

        if self.energy:
            if self.gas.n_reactions > 0:
                heat_of_reaction_from_gas = -np.dot(self.gas.net_rates_of_progress, self.gas.delta_enthalpy)
            else:
                heat_of_reaction_from_gas = 0.

            if self.surf.n_reactions > 0:
                heat_from_reaction_from_surface = -np.dot(self.surf.net_rates_of_progress, self.surf.delta_enthalpy)
            else:
                heat_from_reaction_from_surface = 0.

            dy[-1] = (heat_of_reaction_from_gas + self.alfa * heat_from_reaction_from_surface) / (self.gas.density * self.gas.cp_mass)

        return dy

    def initial_condition(self):
        if self.initial_mass is None:
            self.initial_mass = self.gas.density * self.volume

        return np.block([self.initial_mass_fraction, self.initial_mass, self.initial_coverage, self.initial_temperature])

    def solve(self, tspan, time_ud):
        self.tspan, self.sol = self._solve_ode(self.equations,
                                               self.initial_condition(),
                                               self.uc.convert_to_seconds(tspan, time_ud),
                                               atol=self.atol,
                                               rtol=self.rtol,
                                               verbosity=self.verbosity)

        self.y_sol = self.sol[:, :self.gas.n_species]
        self.x_sol = np.zeros_like(self.y_sol)

        for i in range(0, len(self.tspan)):
            self.x_sol[i, :self.gas.n_species] = self._convert_mass_fraction_to_mole_fraction(self.y_sol[i, :self.gas.n_species])

        self.mass_sol = self.sol[:, self.gas.n_species]
        self.coverage_sol = self.sol[:, self.gas.n_species + 1:self.gas.n_species + 1 + self.surf.n_species]
        self.temperature_sol = self.sol[:, -1]

        self.is_solved = True
        return self.sol
