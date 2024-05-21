from asali.reactors.basic import BasicReactor

import numpy as np

from asali.utils.input_parser import ReactorType


class CstrReactor(BasicReactor):
    def __init__(self, cantera_input_file, gas_phase_name, surface_phase_name):
        """
        Class representing CSTR reactor model
        :param cantera_input_file: Cantera file path
        :param gas_phase_name: Cantera gas phase name
        :param surface_phase_name: Cantera interface phase name
        """
        super().__init__(cantera_input_file=cantera_input_file, gas_phase_name=gas_phase_name,
                         surface_phase_name=surface_phase_name)
        self.solution_parser.reactor_type = ReactorType.CSTR

        self.is_mass_flow_rate = True

        self.volume = 0.
        self.inlet_mass_flow_rate = 0.
        self.inlet_volumetric_flow_rate = 0.
        self.inlet_temperature = 0.

        self.inlet_mass_fraction = None
        self.inlet_mole_fraction = None

    def set_volume(self, value, unit_dimension):
        """
        Set volume
        :param value: Volume value
        :param unit_dimension: Volume unit dimension
        :return: Volume in [m3]
        """
        self.volume = self.uc.convert_to_cubic_meter(value, unit_dimension)
        return self.volume

    def set_mass_flow_rate(self, value, unit_dimension):
        """
        Set mass flow rate
        :param value: Mass flow rate value
        :param unit_dimension: Mass flow rate unit dimension
        :return: Mass flow rate in [kg/s]
        """
        self.inlet_mass_flow_rate = self.uc.convert_to_kg_per_seconds(value, unit_dimension)
        self.inlet_volumetric_flow_rate = 0.
        self.is_mass_flow_rate = True
        return self.inlet_mass_flow_rate

    def set_volumetric_flow_rate(self, value, unit_dimension):
        """
        Set volumetric flow rate
        :param value: Volumetric flow rate value
        :param unit_dimension: Volumetric flow rate unit dimension
        :return: Volumetric flow rate in [m3/s]
        """
        self.inlet_volumetric_flow_rate = self.uc.convert_to_cubic_meter_per_seconds(value, unit_dimension)
        self.inlet_mass_flow_rate = 0.
        self.is_mass_flow_rate = False
        return self.inlet_volumetric_flow_rate

    def set_inlet_mass_fraction(self, value):
        """
        Set inlet mass fraction
        :param value: Mass fraction
        :return: Inlet mass fraction
        """
        self.gas.Y = value
        self.inlet_mass_fraction = self.gas.Y
        self.inlet_mole_fraction = self.gas.X
        return self.inlet_mass_fraction

    def set_inlet_mole_fraction(self, value):
        """
        Set inlet mole fraction
        :param value: Mole fraction
        :return: Inlet mole fraction
        """
        self.gas.X = value
        self.inlet_mass_fraction = self.gas.Y
        self.inlet_mole_fraction = self.gas.X
        return self.inlet_mole_fraction

    def set_inlet_temperature(self, value, unit_dimension):
        """
        Set inlet temperature
        :param value: Temperature value
        :param unit_dimension: Temperature unit dimension
        :return: Temperature in [K]
        """
        self.inlet_temperature = self.uc.convert_to_kelvin(value, unit_dimension)
        return self.inlet_temperature

    def equations(self, t, y):
        """
        Function representing the model
        :param t: Independent variable - Time
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

        domega = (self.inlet_mass_flow_rate / self.volume) * (self.inlet_mass_fraction - omega)
        domega = domega + self.gas.molecular_weights * r_gas
        domega = domega + self.alfa * r_from_surface * self.gas.molecular_weights
        domega = domega / self.gas.density

        dz = r_surface / self.surf.site_density

        dT = 0.0
        if self.energy:
            q_from_gas = self.get_homogeneous_heat_of_reaction()
            q_from_surface = self.get_heterogeneous_heat_of_reaction()

            dT = (self.inlet_mass_flow_rate / self.volume) * (self.inlet_temperature - T) / self.gas.density
            dT = dT + (q_from_gas + self.alfa * q_from_surface) / (self.gas.density * self.gas.cp_mass)

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

        return np.block([self.initial_mass_fraction,
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
