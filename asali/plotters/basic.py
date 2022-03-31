from abc import ABC

import matplotlib

import numpy as np


class BasicPlotter(ABC):
    def __init__(self, cls, colormap):
        """
        Abstract class that can be used to plot results
        :param cls: Object of class to be plotted
        :param colormap: String representing the color map
        """
        self.cls = cls
        self.nFig = 0

        cmap = matplotlib.cm.get_cmap(colormap)
        self.colors = cmap(np.linspace(0.2, 1., num=cls.solution_parser.x.size))

        self._mass_fraction = None
        self._mole_fraction = None
        self._coverage = None
        self._temperature = None

        self._mass_fraction_wall = None
        self._mole_fraction_wall = None
        self._temperature_wall = None

        self._x = None
        self._length = None

    def get_mass_fraction(self):
        """
        Mass fraction getter
        :return: Mass fraction
        """
        return self._mass_fraction

    def set_mass_fraction(self, value):
        """
        Mass fraction setter
        :param value: Mass fraction
        :return:
        """
        self._mass_fraction = value

    # Creating a property object for Mass fraction
    mass_fraction = property(get_mass_fraction, set_mass_fraction)

    def get_mass_fraction_wall(self):
        """
        Solid mass fraction getter
        :return: Mass fraction
        """
        return self._mass_fraction_wall

    def set_mass_fraction_wall(self, value):
        """
        Solid mass fraction setter
        :param value: Mass fraction
        :return:
        """
        self._mass_fraction_wall = value

    # Creating a property object for Solid mass fraction
    mass_fraction_wall = property(get_mass_fraction_wall, set_mass_fraction_wall)

    def get_mole_fraction(self):
        """
        Mole fraction getter
        :return: Mole fraction
        """
        return self._mole_fraction

    def set_mole_fraction(self, value):
        """
        Mole fraction setter
        :param value: Mole fraction
        :return:
        """
        self._mole_fraction = value

    # Creating a property object for Mole fraction
    mole_fraction = property(get_mole_fraction, set_mole_fraction)

    def get_mole_fraction_wall(self):
        """
        Solid mole fraction getter
        :return: Mole fraction
        """
        return self._mole_fraction_wall

    def set_mole_fraction_wall(self, value):
        """
        Solid mole fraction setter
        :param value: Mole fraction
        :return:
        """
        self._mole_fraction_wall = value

    # Creating a property object for Solid mole fraction
    mole_fraction_wall = property(get_mole_fraction_wall, set_mole_fraction_wall)

    def get_coverage(self):
        """
        Site fraction getter
        :return: Site fraction
        """
        return self._coverage

    def set_coverage(self, value):
        """
        Site fraction setter
        :param value: Site fraction
        :return:
        """
        self._coverage = value

    # Creating a property object for Site fraction
    coverage = property(get_coverage, set_coverage)

    def get_temperature(self):
        """
        Temperature getter
        :return: Temperature
        """
        return self._temperature

    def set_temperature(self, value):
        """
        Temperature setter
        :param value: Temperature
        :return:
        """
        self._temperature = value

    # Creating a property object for Temperature
    temperature = property(get_temperature, set_temperature)

    def get_temperature_wall(self):
        """
        Solid temperature getter
        :return: Temperature
        """
        return self._temperature_wall

    def set_temperature_wall(self, value):
        """
        Solid temperature setter
        :param value: Temperature
        :return:
        """
        self._temperature_wall = value

    # Creating a property object for Solid temperature
    temperature_wall = property(get_temperature_wall, set_temperature_wall)

    def get_x(self):
        """
        Independent variable getter
        :return: Independent variable
        """
        return self._x

    def set_x(self, value):
        """
        Independent variable setter
        :param value: Independent variable
        :return:
        """
        self._x = value

    # Creating a property object for Independent variable
    x = property(get_x, set_x)

    def get_length(self):
        """
        Length getter
        :return: Length
        """
        return self._length

    def set_length(self, value):
        """
        Length setter
        :param value: Length
        :return:
        """
        self._length = value

    # Creating a property object for Length
    length = property(get_length, set_length)

    def plot_x_vector_and_y_vector(self, plt, xv, yv, xlabel, ylabel):
        """
        Plotting x and y
        :param plt: Matplotlib object
        :param xv: Independent variable vector
        :param yv: Dependent variable vector
        :param xlabel: Independent variable label
        :param ylabel: Dependent variable label
        :return:
        """
        self.nFig = self.nFig + 1
        plt.figure(self.nFig)
        plt.plot(xv, yv)

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

    def plot_x_vector_and_y_matrix(self, plt, xv, ym, xlabel, ylabel, legend, title=None, colors=None):
        """
        Plotting x and array of y by adding legend
        :param plt: Matplotlib object
        :param xv: Independent variable vector
        :param ym: Dependent variable matrix
        :param xlabel: Independent variable label
        :param ylabel: Dependent variable label
        :param legend: Legend list
        :param title: Plot title
        :param colors: List of colors
        :return:
        """
        self.nFig = self.nFig + 1
        plt.figure(self.nFig)
        if colors is None:
            for yv in ym:
                plt.plot(xv, yv)
        else:
            for i, yv in enumerate(ym):
                plt.plot(xv, yv, color=colors[i])

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.legend(legend, loc='best').get_frame().set_linewidth(0.0)

        if title is not None:
            plt.title(title, loc='center')
