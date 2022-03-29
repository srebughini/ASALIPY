from abc import abstractmethod, ABC

from asali.utils.input_parser import InputParser
from asali.utils.numerical_solvers import NumericalSolvers
from asali.utils.solution_parser import SolutionParser
from asali.utils.unit_converter import UnitConverter

import cantera as ct


class BasicReactor(ABC):
    def __init__(self, cantera_input_file, gas_phase_name, surface_phase_name):
        """
        Abstract class that represents reactor models
        :param cantera_input_file: Cantera file path
        :param gas_phase_name: Cantera gas phase name
        :param surface_phase_name: Cantera interface phase name
        """
        self.uc = UnitConverter()
        self.cantera_input_file = cantera_input_file
        self.gas_phase_name = gas_phase_name
        self.surface_phase_name = surface_phase_name
        self.gas = ct.Solution(cantera_input_file, gas_phase_name)
        self.surf = ct.Interface(cantera_input_file, surface_phase_name, [self.gas])
        self.numerical_solver = NumericalSolvers()
        self.solution_parser = SolutionParser()
        self.solution_parser.gas = self.gas
        self.solution_parser.surf = self.surf

        self.pressure = 0.
        self.alfa = 0.

        self.initial_mass_fraction = None
        self.initial_mole_fraction = None
        self.initial_temperature = None
        self.initial_coverage = None

        self.energy = False

    @abstractmethod
    def equations(self, t, y):
        pass

    @abstractmethod
    def initial_condition(self):
        pass

    @abstractmethod
    def solve(self, tspan, time_ud):
        pass

    def set_energy(self, value):
        """
        Enable/Disable energy balance
        :param value: Variable to enable/disable energy balance
        :return: Bool for energy balance
        """
        self.energy = InputParser.true_parser(value)
        return self.energy

    def set_catalytic_load(self, value, unit_dimension):
        """
        Set catalytic load
        :param value: Catalytic load value
        :param unit_dimension: Catalytic load unit dimension
        :return: Catalytic load in [1/m]
        """
        self.alfa = self.uc.convert_to_one_over_meter(value, unit_dimension)
        return self.alfa

    def set_pressure(self, value, unit_dimension):
        """
        Set pressure
        :param value: Pressure value
        :param unit_dimension: Pressure unit dimension
        :return: Pressure in [Pa]
        """
        self.pressure = self.uc.convert_to_pascal(value, unit_dimension)
        return self.pressure

    def set_initial_mass_fraction(self, value):
        """
        Set initial mass fraction
        :param value: Mass fraction
        :return: Initial mole fraction
        """
        self.gas.Y = value
        self.initial_mass_fraction = self.gas.Y
        self.initial_mole_fraction = self.gas.X
        return self.gas.Y

    def set_initial_mole_fraction(self, value):
        """
        Set initial mole fraction
        :param value: mole fraction
        :return: Initial mass fraction
        """
        self.gas.X = value
        self.initial_mass_fraction = self.gas.Y
        self.initial_mole_fraction = self.gas.X
        return self.gas.X

    def set_initial_coverage(self, value):
        """
        Set initial coverage
        :param value: coverage
        :return: Initial coverage
        """
        self.surf.coverages = value
        self.initial_coverage = self.surf.coverages
        return self.surf.coverages

    def set_initial_temperature(self, value, unit_dimension):
        """
        Set initial temperature
        :param value: Temperature value
        :param unit_dimension: Temperature unit dimension
        :return: Initial temperature in [K]
        """
        self.initial_temperature = self.uc.convert_to_kelvin(value, unit_dimension)
        return self.initial_temperature

    def set_integration_parameters(self, atol, rtol, verbosity):
        """
        Set resolution parameters
        :param atol: Absolute tolerance
        :param rtol: Relative tolerance
        :param verbosity: Solver verbosity
        :return: Set values
        """
        return self.set_absolute_tolerance(atol), self.set_relative_tolerance(rtol), self.set_verbosity(verbosity)

    def set_absolute_tolerance(self, atol):
        """
        Set absolute tolerance
        :param atol: Absolute tolerance
        :return: Set Absolute tolerance
        """
        self.numerical_solver.atol = atol
        return self.numerical_solver.atol

    def set_relative_tolerance(self, rtol):
        """
        Set relative tolerance
        :param rtol: Relative tolerance
        :return: Set Relative tolerance
        """
        self.numerical_solver.rtol = rtol
        return self.numerical_solver.rtol

    def set_verbosity(self, verbosity):
        """
        Set solver verbosity
        :param verbosity: Solver verbosity
        :return:
        """
        if InputParser.true_parser(verbosity):
            self.numerical_solver.verbosity = 20
            return self.numerical_solver.verbosity

        return self.numerical_solver.verbosity

    def get_time(self, ud):
        """
        Get time vector
        :param ud: Unit dimension
        :return: Time vector in the request unit dimension
        """
        return self.uc.convert_from_seconds(self.solution_parser.x, ud)

    def get_results(self):
        """
        Get results
        :return: Vector/Matrix representing the results
        """
        return self.solution_parser.y

    def get_mass_fraction(self, index=None):
        """
        Get mass fraction
        :param index: Index of the species to be extracted
        :return: Vector/Matrix representing the resulting mass fraction
        """
        mass_fraction = self.solution_parser.get_mass_fraction()
        if index is None:
            return mass_fraction

        return mass_fraction[index]

    def get_mole_fraction(self, index=None):
        """
        Get mole fraction
        :param index: Index of the species to be extracted
        :return: Vector/Matrix representing the resulting mole fraction
        """
        mole_fraction = self.solution_parser.get_mole_fraction()
        if index is None:
            return mole_fraction

        return mole_fraction[index]

    def get_temperature(self, index=None):
        """
        Get temperature
        :param index: Index of the axial/time position
        :return: Vector/Matrix representing the resulting temperature at a fixed axial/time position
        """
        temperature = self.solution_parser.get_temperature()
        if index is None:
            return temperature

        return temperature[index]

    def get_coverage(self, index=None):
        """
        Get coverage
        :param index: Index of the coverage to be extracted
        :return: Vector/Matrix representing the resulting coverage
        """
        coverage = self.solution_parser.get_coverage()
        if index is None:
            return coverage

        return coverage[index]
