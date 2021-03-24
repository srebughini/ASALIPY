from abc import abstractmethod, ABC
from enum import IntEnum
from asali.utils.unit_converter import UnitConverter
from assimulo.problem import Implicit_Problem
from assimulo.problem import Explicit_Problem
from assimulo.solvers import IDA
from assimulo.solvers import CVode

import cantera as ct
import numpy as np


class ReactorType(IntEnum):
    BATCH = 0
    CSTR = 1
    PSEUDOHOMOGENEOUSPFR = 2
    HETEROGENEOUSPRF = 3


class ResolutionMethod(IntEnum):
    STEADYSTATE = 0
    TRANSIENT = 1


class ReactorSection(IntEnum):
    CIRCLE = 0
    SQUARE = 1
    TRIANGLE = 2


class BasicReactor(ABC):
    def __init__(self, cantera_input_file, gas_phase_name, surface_phase_name):
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

        self.atol = 1.e-12
        self.rtol = 1.e-07

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
    def _true_parser(v):
        return v in ['yes', 'y', '1', 1, True, 'true', 'True']

    @staticmethod
    def _resolution_parser(method):
        if method in ["steadystate", "ss", "SteadyState", "STEADYSTATE", "steady state", "Steady State"]:
            return ResolutionMethod.STEADYSTATE

        if method in ["transient", "tt", "Transient", "TRANSIENT"]:
            return ResolutionMethod.TRANSIENT

        raise Exception("ASALI::ERROR::Unknown resolution method")

    @staticmethod
    def _section_parser(method):
        if method in ["circle", "Circle", "CIRCLE"]:
            return ReactorSection.CIRCLE

        if method in ["square", "Square", "SQUARE"]:
            return ReactorSection.SQUARE

        if method in ["triangle", "Triangle", "TRIANGLE"]:
            return ReactorSection.TRIANGLE

        raise Exception("ASALI::ERROR::Unknown resolution method")

    @staticmethod
    def _solve_ode(f, y0, tspan, atol=1e-12, rtol=1e-7, verbosity=50):
        ode = Explicit_Problem(f, y0)

        sim_ode = CVode(ode)
        sim_ode.atol = atol
        sim_ode.rtol = rtol
        sim_ode.verbosity = verbosity

        sol = np.zeros([len(tspan), y0.size], dtype=np.float64)

        sol[0, :] = y0
        for i, t in enumerate(tspan):
            _, y = sim_ode.simulate(t, 2)
            sol[i, :] = y[-1, :]

        return tspan, sol

    @staticmethod
    def _solve_dae(ode_equations, dae_equations, residuals, y0, tspan, alg, t_ode=1e06, atol=1e-12, rtol=1e-7, verbosity=50):
        ode = Explicit_Problem(ode_equations, y0)

        sim_ode = CVode(ode)
        sim_ode.atol = atol
        sim_ode.rtol = rtol
        sim_ode.verbosity = verbosity

        t, y_ode = sim_ode.simulate(t_ode, 2)

        dae = Implicit_Problem(residuals, y_ode[-1, :], dae_equations(0., y_ode[-1, :]), 0.0)

        dae.algvar = alg

        sim_dae = IDA(dae)
        sim_dae.atol = atol
        sim_dae.rtol = rtol
        sim_dae.verbosity = verbosity

        sol = np.zeros([len(tspan), y_ode[-1, :].size], dtype=np.float64)

        sol[0, :] = y0
        for i, t in enumerate(tspan):
            _, y, _ = sim_dae.simulate(t, 2)
            sol[i, :] = y[-1, :]

        return tspan, sol

    def _convert_mass_fraction_to_mole_fraction(self, y):
        self.gas.Y = y
        return self.gas.X

    def _convert_mole_fraction_to_mass_fraction(self, x):
        self.gas.X = x
        return self.gas.Y

    def set_mass_flow_rate(self, value, unit_dimension):
        if self.reactor_type == ReactorType.BATCH:
            print("ASALI::WARNING::Mass flow rate ignored for BATCH reactor")
        self.inlet_mass_flow_rate = self.uc.convert_to_kg_per_seconds(value, unit_dimension)
        self.inlet_volumetric_flow_rate = 0.
        self.is_mass_flow_rate = True
        return self.inlet_mass_flow_rate

    def set_volumetric_flow_rate(self, value, unit_dimension):
        if self.reactor_type == ReactorType.BATCH:
            print("ASALI::WARNING::Volumetric flow rate ignored for BATCH reactor")
        self.inlet_volumetric_flow_rate = self.uc.convert_to_cubic_meter_per_seconds(value, unit_dimension)
        self.inlet_mass_flow_rate = 0.
        self.is_mass_flow_rate = False
        return self.inlet_volumetric_flow_rate

    def set_inlet_mass_fraction(self, value):
        if self.reactor_type == ReactorType.BATCH:
            print("ASALI::WARNING::Inlet mass ignored for BATCH reactor")
        self.gas.Y = value
        self.inlet_mass_fraction = self.gas.Y
        self.inlet_mole_fraction = self.gas.X
        return self.inlet_mass_fraction

    def set_inlet_mole_fraction(self, value):
        if self.reactor_type == ReactorType.BATCH:
            print("ASALI::WARNING::Inlet mole ignored for BATCH reactor")
        self.gas.X = value
        self.inlet_mass_fraction = self.gas.Y
        self.inlet_mole_fraction = self.gas.X
        return self.inlet_mole_fraction

    def set_inlet_temperature(self, value, unit_dimension):
        if self.reactor_type == ReactorType.BATCH:
            print("ASALI::WARNING::Inlet temperature ignored for BATCH reactor")
        self.inlet_temperature = self.uc.convert_to_kelvin(value, unit_dimension)
        return self.inlet_temperature

    def set_catalytic_load(self, value, unit_dimension):
        self.alfa = self.uc.convert_to_one_over_meter(value, unit_dimension)
        return self.alfa

    def set_energy(self, value):
        self.energy = self._true_parser(value)
        return self.energy

    def set_integration_parameters(self, atol, rtol):
        self.atol = atol
        self.rtol = rtol

    def set_volume(self, value, unit_dimension):
        if self.reactor_type == ReactorType.PSEUDOHOMOGENEOUSPFR:
            print("ASALI::WARNING::Volume ignored for PSEUDOHOMOGENEOUS PLUG FLOW reactor")

        if self.reactor_type == ReactorType.HETEROGENEOUSPRF:
            print("ASALI::WARNING::Volume ignored for HETEROGENEOUS PLUG FLOW reactor")

        self.volume = self.uc.convert_to_cubic_meter(value, unit_dimension)
        return self.volume

    def set_pressure(self, value, unit_dimension):
        self.pressure = self.uc.convert_to_pascal(value, unit_dimension)
        return self.pressure

    def set_diameter(self, value, unit_dimension):
        if self.reactor_type == ReactorType.BATCH:
            print("ASALI::WARNING::Diameter ignored for BATCH reactor")

        if self.reactor_type == ReactorType.CSTR:
            print("ASALI::WARNING::Diameter ignored for CSTR reactor")

        diameter = self.uc.convert_to_meter(value, unit_dimension)
        self.area = np.pi * 0.25 * diameter
        return diameter

    def set_length(self, value, unit_dimension):
        if self.reactor_type == ReactorType.BATCH:
            print("ASALI::WARNING::Length ignored for BATCH reactor")

        if self.reactor_type == ReactorType.CSTR:
            print("ASALI::WARNING::Length ignored for CSTR reactor")

        if not isinstance(value, list):
            length = self.uc.convert_to_meter(value, unit_dimension)
            self.length = np.linspace(0, length, num=10)
        else:
            if value[0] != 0.:
                length = np.zeros([len(value) + 1], dtype=np.float64)
                length[1:] = value
                self.length = self.uc.convert_to_meter(length, unit_dimension)
            else:
                self.length = self.uc.convert_to_meter(value, unit_dimension)

        return self.length

    def set_initial_mass_fraction(self, value):
        self.gas.Y = value
        self.initial_mass_fraction = self.gas.Y
        self.initial_mole_fraction = self.gas.X
        return self.gas.Y

    def set_initial_mole_fraction(self, value):
        self.gas.X = value
        self.initial_mass_fraction = self.gas.Y
        self.initial_mole_fraction = self.gas.X
        return self.gas.X

    def set_initial_coverage(self, value):
        self.surf.coverages = value
        self.initial_coverage = self.surf.coverages
        return self.surf.coverages

    def set_initial_temperature(self, value, unit_dimension):
        self.temperature = self.uc.convert_to_kelvin(value, unit_dimension)
        self.initial_temperature = self.uc.convert_to_kelvin(value, unit_dimension)
        return self.temperature

    def set_initial_solid_temperature(self, value, unit_dimension):
        if self.reactor_type == ReactorType.BATCH:
            print("ASALI::WARNING::Solid temperature ignored for BATCH reactor")

        if self.reactor_type == ReactorType.CSTR:
            print("ASALI::WARNING::Solid temperature ignored for CSTR reactor")

        if self.reactor_type == ReactorType.PSEUDOHOMOGENEOUSPFR:
            print("ASALI::WARNING::Solid temperature ignored for PSEUDOHOMOGENEOUS PLUG FLOW reactor")

        self.solid_temperature = self.uc.convert_to_kelvin(value, unit_dimension)
        return self.solid_temperature

    def set_resolution_method(self, method):
        if self.reactor_type == ReactorType.BATCH:
            print("ASALI::WARNING::Resolution method ignored for BATCH reactor")

        if self.reactor_type == ReactorType.CSTR:
            print("ASALI::WARNING::Resolution method ignored for CSTR reactor")

        if self.reactor_type == ReactorType.HETEROGENEOUSPRF:
            print("ASALI::WARNING::Resolution method ignored for HETEROGENEOUSPRF PLUG FLOW reactor")

        self.resolution_method = self._resolution_parser(method)
        return self.resolution_method

    def set_reactor_section(self, method):
        if self.reactor_type == ReactorType.BATCH:
            print("ASALI::WARNING::Reactor section ignored for BATCH reactor")

        if self.reactor_type == ReactorType.CSTR:
            print("ASALI::WARNING::Reactor section ignored for CSTR reactor")

        if self.reactor_type == ReactorType.PSEUDOHOMOGENEOUSPFR:
            print("ASALI::WARNING::Reactor section ignored for PSEUDOHOMOGENEOUS PLUG FLOW reactor")

        self.reactor_section = self._section_parser(method)
        return self.reactor_section

    def set_solid_density(self, value, unit_dimension):
        if self.reactor_type == ReactorType.BATCH:
            print("ASALI::WARNING::Solid density ignored for BATCH reactor")

        if self.reactor_type == ReactorType.CSTR:
            print("ASALI::WARNING::Solid density ignored for CSTR reactor")

        if self.reactor_type == ReactorType.PSEUDOHOMOGENEOUSPFR:
            print("ASALI::WARNING::Solid density ignored for PSEUDOHOMOGENEOUS PLUG FLOW reactor")

        self.solid_rho = self.uc.convert_to_kg_per_cubic_meter(value, unit_dimension)
        return self.solid_rho

    def set_solid_specific_heat(self, value, unit_dimension):
        if self.reactor_type == ReactorType.BATCH:
            print("ASALI::WARNING::Solid specific heat ignored for BATCH reactor")

        if self.reactor_type == ReactorType.CSTR:
            print("ASALI::WARNING::Solid specific heat ignored for CSTR reactor")

        if self.reactor_type == ReactorType.PSEUDOHOMOGENEOUSPFR:
            print("ASALI::WARNING::Solid specific heat ignored for PSEUDOHOMOGENEOUS PLUG FLOW reactor")

        self.solid_cp = self.uc.convert_to_joule_per_kg_per_kelvin(value, unit_dimension)
        return self.solid_cp

    def set_solid_thermal_conductivity(self, value, unit_dimension):
        if self.reactor_type == ReactorType.BATCH:
            print("ASALI::WARNING::Solid thermal conductivity ignored for BATCH reactor")

        if self.reactor_type == ReactorType.CSTR:
            print("ASALI::WARNING::Solid thermal conductivity ignored for CSTR reactor")

        if self.reactor_type == ReactorType.PSEUDOHOMOGENEOUSPFR:
            print("ASALI::WARNING::Solid thermal conductivity ignored for PSEUDOHOMOGENEOUS PLUG FLOW reactor")

        self.solid_k = self.uc.convert_to_watt_per_meter_per_kelvin(value, unit_dimension)
        return self.solid_k

    def get_time(self, ud):
        return self.uc.convert_from_seconds(self.tspan, ud)

    def get_results(self):
        return self.sol

    def get_mass_fraction(self, index=None):
        if index is None:
            return self.y_sol

        return self.y_sol[index]

    def get_mole_fraction(self, index=None):
        if index is None:
            return self.x_sol

        return self.x_sol[index]

    def get_temperature(self, index=None):
        if index is None:
            return self.temperature_sol

        return self.temperature_sol[index]

    def get_coverage(self, index=None):
        if index is None:
            return self.coverage_sol

        return self.coverage_sol[index]
