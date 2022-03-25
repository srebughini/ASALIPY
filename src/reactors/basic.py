from abc import abstractmethod, ABC
from enum import IntEnum
from asali.utils.unit_converter import UnitConverter
from assimulo.problem import Implicit_Problem
from assimulo.problem import Explicit_Problem
from assimulo.solvers import IDA
from assimulo.solvers import CVode

import cantera as ct
import numpy as np
import warnings


class ReactorType(IntEnum):
    """
    Reactor type
    """
    BATCH = 0
    CSTR = 1
    PSEUDOHOMOGENEOUSPFR = 2
    HETEROGENEOUSPRF = 3


class ResolutionMethod(IntEnum):
    """
    Resolution method
    """
    STEADYSTATE = 0
    TRANSIENT = 1


class ReactorSection(IntEnum):
    """
    Reactor section
    """
    CIRCLE = 0
    SQUARE = 1
    TRIANGLE = 2


class BasicReactor(ABC):
    def __init__(self, cantera_input_file, gas_phase_name, surface_phase_name):
        """
        Abstract class that represents reactor models
        :param cantera_input_file: Cantera file path
        :param gas_phase_name: Cantera gas phase name
        :param surface_phase_name: Cantera interface phase name
        """
        self.uc = UnitConverter()
        self.gas = ct.Solution(cantera_input_file, gas_phase_name)
        self.surf = ct.Interface(cantera_input_file, surface_phase_name, [self.gas])

        self.energy = False
        self.is_mass_flow_rate = True
        self.is_solved = False
        self.alfa = 0.
        self.volume = 0.
        self.pressure = 0.
        self.temperature = 0.
        self.area = 0.
        self.inlet_mass_flow_rate = 0.
        self.inlet_volumetric_flow_rate = 0.
        self.inlet_temperature = 0.
        self.solid_rho = 0.
        self.solid_k = 0.
        self.solid_cp = 0.
        self.solid_temperature = 0.

        self.inert_specie_index = None
        self.inlet_mass_fraction = None
        self.inlet_mole_fraction = None
        self.sol = None
        self.tspan = None
        self.x_sol = None
        self.y_sol = None
        self.temperature_sol = None
        self.coverage_sol = None
        self.length = None
        self.initial_mass_fraction = None
        self.initial_mole_fraction = None
        self.initial_temperature = None
        self.initial_coverage = None

        self.atol = 1e-06
        self.rtol = 1e-06
        self.verbosity = 50

        self.reactor_type = None
        self.resolution_method = None
        self.reactor_section = None

    @abstractmethod
    def equations(self, t, y):
        pass

    @abstractmethod
    def initial_condition(self):
        pass

    @abstractmethod
    def solve(self, tspan, time_ud):
        pass

    @staticmethod
    def _print_warning(message, verbosity):
        """
        Print warning of screen
        ::param message: Message to be printed
        :return:
        """
        warning_format = 'ASALI::WARNING::{}'
        if verbosity > 20:
            warnings.warn(warning_format.format(message))

    @staticmethod
    def _raise_error(message):
        """
        Raise error
        :param message: Message to be printed
        :return:
        """
        error_format = 'ASALI::ERROR::{}'
        raise Exception(error_format.format(message))

    @staticmethod
    def _true_parser(v):
        """
        Parser of True
        :param v: Variable to be parsed
        :return: Bool rapresenting the variable
        """
        return v in ['yes', 'y', '1', 1, True, 'true', 'True']

    @staticmethod
    def _resolution_parser(method):
        """
        Parser of the resolution method
        :param method: Resolution method as string
        :return: ResolutionMethod object
        """
        if method in ["steadystate", "ss", "SteadyState", "STEADYSTATE", "steady state", "Steady State"]:
            return ResolutionMethod.STEADYSTATE

        if method in ["transient", "tt", "Transient", "TRANSIENT"]:
            return ResolutionMethod.TRANSIENT

        BasicReactor._raise_error("Unknown resolution method")

    @staticmethod
    def _section_parser(method):
        """
        Parser of the reactor section
        :param method: Reactor section as string
        :return: ReactorSection
        """
        if method in ["circle", "Circle", "CIRCLE"]:
            return ReactorSection.CIRCLE

        if method in ["square", "Square", "SQUARE"]:
            return ReactorSection.SQUARE

        if method in ["triangle", "Triangle", "TRIANGLE"]:
            return ReactorSection.TRIANGLE

        BasicReactor._raise_error("Unknown section shape")

    @staticmethod
    def _solve_ode(f, y0, tspan, atol, rtol, verbosity):
        """
        Solve ODE system
        :param f: Function describing the ODE system
        :param y0: Initial condition vector
        :param tspan: Integration vector
        :param atol: Absolute tolerance
        :param rtol: Relative tolerance
        :param verbosity: Verbosity for ASSIMULO
        :return: tspan: Independent variable vector
                 sol: Solved dependent variable matrix
        """
        ode = Explicit_Problem(f, y0)

        sim_ode = CVode(ode)
        sim_ode.atol = atol
        sim_ode.rtol = rtol
        sim_ode.verbosity = verbosity

        _, sol = sim_ode.simulate(tspan[-1], 0, tspan)
        return tspan, sol

    @staticmethod
    def _solve_dae(ode_equations, dae_equations, residuals, y0, tspan, alg, atol, rtol, verbosity):
        """
        Solve DAE system
        :param ode_equations: Function describing the ODE system for initial conditions estimation
        :param dae_equations: Function describing the DAE system
        :param residuals: Function describing the DAE system residuals
        :param y0: Initial condition vector
        :param tspan: Integration vector
        :param alg: Algebraic/Differential equation vector
        :param atol: Absolute tolerance
        :param rtol: Relative tolerance
        :param verbosity: Verbosity for ASSIMULO
        :return: tspan: Independent variable vector
                 sol: Solved dependent variable matrix
        """
        ode = Explicit_Problem(ode_equations, y0)

        sim_ode = CVode(ode)
        sim_ode.atol = atol
        sim_ode.rtol = rtol
        sim_ode.verbosity = verbosity
        sim_ode.linear_solver = 'SPGMR'

        t, y_ode = sim_ode.simulate(1e06, 2)

        dae = Implicit_Problem(residuals, y_ode[-1, :], dae_equations(0., y_ode[-1, :]), 0.0)

        dae.algvar = alg

        sim_dae = IDA(dae)
        sim_dae.atol = atol
        sim_dae.rtol = rtol
        sim_dae.verbosity = verbosity
        sim_dae.suppress_alg = True
        sim_dae.dqtype = 'FORWARD'
        sim_dae.make_consistent('IDA_YA_YDP_INIT')

        _, y, _ = sim_dae.simulate(tspan[-1], 0, tspan)

        sol = np.asarray(y)
        sol[0, :] = y0

        return tspan, sol

    def _convert_mass_fraction_to_mole_fraction(self, y):
        """
        Convert mass fraction to mole fraction
        :param y: Mass fraction
        :return: Vector representing the mole fraction
        """
        self.gas.Y = y
        return self.gas.X

    def _convert_mole_fraction_to_mass_fraction(self, x):
        """
        Convert mole fraction to mass fraction
        :param x: Mole fraction
        :return: Vector representing the mass fraction
        """
        self.gas.X = x
        return self.gas.Y

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

    def set_catalytic_load(self, value, unit_dimension):
        """
        Set catalytic load
        :param value: Catalytic load value
        :param unit_dimension: Catalytic load unit dimension
        :return: Catalytic load in [1/m]
        """
        self.alfa = self.uc.convert_to_one_over_meter(value, unit_dimension)
        return self.alfa

    def set_volume(self, value, unit_dimension):
        """
        Set volume
        :param value: Volume value
        :param unit_dimension: Volume unit dimension
        :return: Volume in [m3]
        """
        self.volume = self.uc.convert_to_cubic_meter(value, unit_dimension)
        return self.volume

    def set_pressure(self, value, unit_dimension):
        """
        Set pressure
        :param value: Pressure value
        :param unit_dimension: Pressure unit dimension
        :return: Pressure in [Pa]
        """
        self.pressure = self.uc.convert_to_pascal(value, unit_dimension)
        return self.pressure

    def set_diameter(self, value, unit_dimension):
        """
        Set diameter
        :param value: Diameter value
        :param unit_dimension: Diameter unit dimension
        :return: Diameter in [m]
        """
        diameter = self.uc.convert_to_meter(value, unit_dimension)
        self.area = np.pi * 0.25 * np.square(diameter)
        return diameter

    def set_length(self, value, unit_dimension):
        """
        Set length
        :param value: Length value
        :param unit_dimension: Length unit dimension
        :return: Discretize length in [m]
        """
        if isinstance(value, (list, np.ndarray)):
            if value[0] != 0.:
                length = np.zeros([len(value) + 1], dtype=np.float64)
                length[1:] = value
                self.length = self.uc.convert_to_meter(length, unit_dimension)
            else:
                self.length = self.uc.convert_to_meter(value, unit_dimension)
        else:
            length = self.uc.convert_to_meter(value, unit_dimension)
            self.length = np.linspace(0, length, num=10)

        return self.length

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
        self.temperature = self.uc.convert_to_kelvin(value, unit_dimension)
        self.initial_temperature = self.uc.convert_to_kelvin(value, unit_dimension)
        return self.temperature

    def set_initial_solid_temperature(self, value, unit_dimension):
        """
        Set initial solid temperature
        :param value: Temperature value
        :param unit_dimension: Temperature unit dimension
        :return: Initial solid temperature in [K]
        """
        self.solid_temperature = self.uc.convert_to_kelvin(value, unit_dimension)
        return self.solid_temperature

    def set_energy(self, value):
        """
        Enable/Disable energy balance
        :param value: Variable to enable/disable energy balance
        :return: Bool for energy balance
        """
        self.energy = self._true_parser(value)
        return self.energy

    def set_resolution_method(self, method):
        """
        Set resolution method
        :param method: Resolution method as string
        :return: ResolutionMethod object
        """
        self.resolution_method = self._resolution_parser(method)
        return self.resolution_method

    def set_reactor_section(self, method):
        """
        Set reactor section
        :param method: Reactor section as string
        :return: ReactorSection
        """
        self.reactor_section = self._section_parser(method)
        return self.reactor_section

    def set_solid_density(self, value, unit_dimension):
        """
        Set solid density
        :param value: Density value
        :param unit_dimension: Density unit dimension
        :return: Solid density in [kg/m3]
        """
        self.solid_rho = self.uc.convert_to_kg_per_cubic_meter(value, unit_dimension)
        return self.solid_rho

    def set_solid_specific_heat(self, value, unit_dimension):
        """
        Set solid specific heat
        :param value: Specific heat value
        :param unit_dimension: Specific heat unit dimension
        :return: Solid specific heat in [J/kg/K]
        """
        self.solid_cp = self.uc.convert_to_joule_per_kg_per_kelvin(value, unit_dimension)
        return self.solid_cp

    def set_solid_thermal_conductivity(self, value, unit_dimension):
        """
        Set solid thermal conductivity
        :param value: Thermal conductivity value
        :param unit_dimension: Thermal conductivity unit dimension
        :return: Solid thermal conductivity in [W/m/K]
        """
        self.solid_k = self.uc.convert_to_watt_per_meter_per_kelvin(value, unit_dimension)
        return self.solid_k

    def set_inert_specie(self, specie_name):
        """
        Set inert specie
        :param specie_name: Specie name
        :return: Specie index
        """
        self.inert_specie_index = self.gas.species_index(specie_name)
        return self.inert_specie_index

    def set_integration_parameters(self, atol, rtol, verbosity):
        """
        Set resolution parameters
        :param atol: Absolute tolerance
        :param rtol: Relative tolerance
        :param verbosity: Solver verbosity
        :return:
        """
        self.set_absolute_tolerance(atol)
        self.set_relative_tolerance(rtol)
        self.set_verbosity(verbosity)

    def set_absolute_tolerance(self, atol):
        """
        Set absolute tolerance
        :param atol: Absolute tolerance
        :return:
        """
        self.atol = atol

    def set_relative_tolerance(self, rtol):
        """
        Set relative tolerance
        :param rtol: Relative tolerance
        :return:
        """
        self.rtol = rtol

    def set_verbosity(self, verbosity):
        """
        Set solver verbosity
        :param verbosity: Solver verbosity
        :return:
        """
        if self._true_parser(verbosity):
            self.verbosity = 20
        else:
            self.verbosity = 50

    def get_time(self, ud):
        """
        Get time vector
        :param ud: Unit dimension
        :return: Time vector in the request unit dimension
        """
        return self.uc.convert_from_seconds(self.tspan, ud)

    def get_results(self):
        """
        Get results
        :return: Vector/Matrix representing the results
        """
        return self.sol

    def get_mass_fraction(self, index=None):
        """
        Get mass fraction
        :param index: Index of the species to be extracted
        :return: Vector/Matrix representing the resulting mass fraction
        """
        if index is None:
            return self.y_sol

        return self.y_sol[index]

    def get_mole_fraction(self, index=None):
        """
        Get mole fraction
        :param index: Index of the species to be extracted
        :return: Vector/Matrix representing the resulting mole fraction
        """
        if index is None:
            return self.x_sol

        return self.x_sol[index]

    def get_temperature(self, index=None):
        """
        Get temperature
        :param index: Index of the axial/time position
        :return: Vector/Matrix representing the resulting temperature at a fixed axial/time position
        """
        if index is None:
            return self.temperature_sol

        return self.temperature_sol[index]

    def get_coverage(self, index=None):
        """
        Get coverage
        :param index: Index of the coverage to be extracted
        :return: Vector/Matrix representing the resulting coverage
        """
        if index is None:
            return self.coverage_sol

        return self.coverage_sol[index]
