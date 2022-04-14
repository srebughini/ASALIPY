import numpy as np


class SteadyStatePseudoHomogeneous1DReactor:
    def __init__(self, gas, surf, pressure, alfa, energy, inlet_mass_flow_rate, area):
        """
        Class representing PseudoHomogeneous 1D reactor model
        :param gas: Cantera gas phase object
        :param surf: Cantera surface phase object
        :param pressure: Reactor pressure
        :param alfa: Reactor catalytic load
        :param energy: Enable/Disable energy balance
        :param inlet_mass_flow_rate: Reactor inlet mass flow rate
        :param area: Reactor section area
        """
        self.gas = gas
        self.surf = surf
        self.pressure = pressure
        self.alfa = alfa
        self.energy = energy
        self.inlet_mass_flow_rate = inlet_mass_flow_rate
        self.area = area

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

        r_gas = self.gas.net_production_rates
        r_from_surface = self.surf.net_production_rates[:self.gas.n_species]
        r_surface = self.surf.net_production_rates[-self.surf.n_species:]

        # Equation of mass
        dy[:self.gas.n_species] = self.gas.molecular_weights * r_gas
        dy[:self.gas.n_species] = dy[:self.gas.n_species] + self.alfa * r_from_surface * self.gas.molecular_weights
        dy[:self.gas.n_species] = dy[:self.gas.n_species] * self.area / self.inlet_mass_flow_rate

        # Equation of coverage
        dy[self.gas.n_species:self.gas.n_species + self.surf.n_species] = 1e03 * (r_surface / self.surf.site_density)

        if self.energy:
            if self.gas.n_reactions > 0:
                q_gas = -np.dot(self.gas.net_rates_of_progress, self.gas.delta_enthalpy)
            else:
                q_gas = 0.

            if self.surf.n_reactions > 0:
                q_from_surface = -np.dot(self.surf.net_rates_of_progress, self.surf.delta_enthalpy)
            else:
                q_from_surface = 0.

            # Equation of energy
            dy[-1] = (q_gas + self.alfa * q_from_surface) * self.area / (self.inlet_mass_flow_rate * self.gas.cp_mass)

        return dy

    def solve(self, numerical_solver, length, initial_conditions):
        """
        Solve Steady State equations
        :param numerical_solver: Numerical solver object
        :param length: Reactor length
        :param initial_conditions: Reactor initial conditions
        :return: Reactor length vector,
                 Matrix representing the solution in terms of composition, coverage and temperature as function of reactor length
        """
        return numerical_solver.solve_ode(self.equations,
                                          initial_conditions,
                                          length)
