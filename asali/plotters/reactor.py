import matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl

from asali.plotters.batch_and_cstr import BatchAndCstrPlotter
from asali.plotters.het1d import Heterogeneous1DReactorPlotter
from asali.plotters.ph1d import PseudoHomogeneous1DReactorPlotter
from asali.utils.input_parser import ReactorType


class ReactorPlotter:
    def __init__(self, cls, colormap="Blues", style='default'):
        """
        Class to plot reactor results
        :param cls: Reactor class object
        :param colormap: String representing the color map
        :param style: String representing the matplotlib style sheet
        :param toolbar: Enable/Disable the toolbar before plotting
        """

        if cls.solution_parser.reactor_type not in [r for r in ReactorType]:
            raise Exception("ASALI::ERROR::Unknown class type: ", str(type(cls)))

        if not cls.solution_parser.is_solved:
            raise Exception("ASALI::ERROR::Nothing to plot, no solution found")

        self.plotter = ReactorPlotter.get_plotter(cls, colormap)

        plt.style.use(style)

    @staticmethod
    def get_plotter(cls, colormap):
        """
        Function to get the correct plotter based on Reactor Type
        :param cls: Reactor class object
        :param colormap: Selected color map
        :return: Plotter object
        """
        if cls.solution_parser.reactor_type == ReactorType.BATCH:
            return BatchAndCstrPlotter(cls, colormap)

        if cls.solution_parser.reactor_type == ReactorType.CSTR:
            return BatchAndCstrPlotter(cls, colormap)

        if cls.solution_parser.reactor_type == ReactorType.PSEUDOHOMOGENEOUSPFR:
            return PseudoHomogeneous1DReactorPlotter(cls, colormap)

        if cls.solution_parser.reactor_type == ReactorType.HETEROGENEOUSPRF:
            return Heterogeneous1DReactorPlotter(cls, colormap)

    @staticmethod
    def set_rc_params(params):
        """
        Set rc parameters for matplotlib
        :param params: Parameters as dict
        :return:
        """
        plt.rcParams.update(params)

    @staticmethod
    def set_default_rc_params():
        """
        Set default rc parameters for matplotlib
        :return:
        """
        mpl.rcdefaults()

    def plot_species_mass_fraction(self, species_names):
        """
        Plotting mass fraction
        :param species_names: List of species to be plotted
        :return:
        """
        self.plotter.plot_species_mass_fraction(plt, species_names)

    def plot_species_mole_fraction(self, species_names):
        """
        Plotting mole fraction
        :param species_names: List of species to be plotted
        :return:
        """
        self.plotter.plot_species_mole_fraction(plt, species_names)

    def plot_coverage(self, coverage_names):
        """
        Plotting coverage
        :param coverage_names: List of coverage species to be plotted
        :return:
        """
        self.plotter.plot_coverage(plt, coverage_names)

    def plot_temperature(self):
        """
        Plotting temperature
        :return:
        """
        self.plotter.plot_temperature(plt)

    @staticmethod
    def show():
        """
        Show plots
        :return:
        """
        plt.show()
