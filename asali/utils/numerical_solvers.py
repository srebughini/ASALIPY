import numpy as np
from assimulo.problem import Implicit_Problem
from assimulo.problem import Explicit_Problem
from assimulo.solvers import IDA
from assimulo.solvers import CVode


class NumericalSolvers:
    def __init__(self):
        """
        Class to solve ODE and DAE systems using SUNDIALS/ASSIMULO
        """
        self._atol = 1.e-06
        self._rtol = 1.e-06
        self._verbosity = 50

    def get_atol(self):
        """
        Absolute tolerance getter
        :return: Absolute tolerance
        """
        return self._atol

    def set_atol(self, value):
        """
        Absolute tolerance setter
        :param value: Absolute tolerance
        :return:
        """
        self._atol = value

    # Creating a property object for Absolute tolerance
    atol = property(get_atol, set_atol)

    def get_rtol(self):
        """
        Relative tolerance getter
        :return: Relative tolerance
        """
        return self._rtol

    def set_rtol(self, value):
        """
        Relative tolerance setter
        :param value: Relative tolerance
        :return:
        """
        self._rtol = value

    # Creating a property object for Relative tolerance
    rtol = property(get_rtol, set_rtol)

    def get_verbosity(self):
        """
        Verbosity getter
        :return: Verbosity
        """
        return self._verbosity

    def set_verbosity(self, value):
        """
        Verbosity setter
        :param value: Verbosity
        :return:
        """
        self._verbosity = value

    # Creating a property object for Verbosity
    verbosity = property(get_verbosity, set_verbosity)

    def solve_ode(self, f, y0, tspan):
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
        sim_ode.atol = self._atol
        sim_ode.rtol = self._rtol
        sim_ode.verbosity = self._verbosity

        t, y = sim_ode.simulate(tspan[-1])
        y = np.asarray(y)
        t = np.asarray(t)

        sol = np.zeros([len(tspan), len(y0)], dtype=np.float64)
        for i in range(0, len(y0)):
            sol[:, i] = NumericalSolvers.interpolation(tspan, t, y[:,i])

        return tspan, sol

    def solve_dae(self, ode_equations, dae_equations, residuals, y0, tspan, alg):
        """
        Solve DAE system
        :param ode_equations: Function describing the ODE system for initial conditions estimation
        :param dae_equations: Function describing the DAE system
        :param residuals: Function describing the DAE system residuals
        :param y0: Initial condition vector
        :param tspan: Integration vector
        :param alg: Algebraic/Differential equation vector
        :return: tspan: Independent variable vector
                 sol: Solved dependent variable matrix
        """
        ode = Explicit_Problem(ode_equations, y0)

        sim_ode = CVode(ode)
        sim_ode.atol = self._atol
        sim_ode.rtol = self._rtol
        sim_ode.verbosity = self._verbosity
        sim_ode.linear_solver = 'SPGMR'

        t, y_ode = sim_ode.simulate(1e06, 2)

        dae = Implicit_Problem(residuals, y_ode[-1, :], dae_equations(0., y_ode[-1, :]), 0.0)

        dae.algvar = alg

        sim_dae = IDA(dae)
        sim_dae.atol = self._atol
        sim_dae.rtol = self._rtol
        sim_dae.verbosity = self._verbosity
        sim_dae.suppress_alg = True
        sim_dae.dqtype = 'FORWARD'
        sim_dae.make_consistent('IDA_YA_YDP_INIT')

        t, y, _ = sim_dae.simulate(tspan[-1])
        y = np.asarray(y)
        t = np.asarray(t)

        sol = np.zeros([len(tspan), len(y0)], dtype=np.float64)
        for i in range(0, len(y0)):
            sol[:, i] = NumericalSolvers.interpolation(tspan, t, y[:,i])

        sol[0, :] = y0

        return tspan, sol

    @staticmethod
    def interpolation(x_target, x, y):
        """
        Numerical interpolation
        :param x_target: Independent variable target
        :param x: Independent variable
        :param y: Dependent variable
        :return: Dependent variable interpolated
        """
        return np.interp(x_target, x, y)
