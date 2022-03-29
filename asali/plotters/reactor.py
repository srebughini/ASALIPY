from asali.plotters.basic import BasicPlotter

import matplotlib.pyplot as plt
import matplotlib
import numpy as np

from asali.utils.input_parser import ReactorType, ResolutionMethod


class ReactorPlotter(BasicPlotter):
    def __init__(self, reactor_cls, colormap="Blues"):
        """
        Class to plot reactor results
        :param reactor_cls: Reactor class object
        :param colormap: String representing the color map
        """
        super().__init__(cls=reactor_cls, colormap=colormap)

        if self.cls.solution_parser.reactor_type not in [r for r in ReactorType]:
            raise Exception("ASALI::ERROR::Unknown class type: ", str(type(self.cls)))

        if not self.cls.solution_parser.is_solved:
            raise Exception("ASALI::ERROR::Nothing to plot, no solution found")

        cmap = matplotlib.cm.get_cmap(self.colormap)

        self.colors = cmap(np.linspace(0.2, 1., num=self.cls.solution_parser.x.size))

        self.mass_fraction = self.cls.solution_parser.get_mass_fraction()
        self.mole_fraction = self.cls.solution_parser.get_mole_fraction()
        self.coverage = self.cls.solution_parser.get_coverage()
        self.temperature = self.cls.solution_parser.get_temperature()

        self.mass_fraction_wall = self.cls.solution_parser.get_mass_fraction_wall()
        self.mole_fraction_wall = self.cls.solution_parser.get_mole_fraction_wall()
        self.temperature_wall = self.cls.solution_parser.get_temperature_wall()

        self.x = self.cls.solution_parser.x
        self.length = self.cls.solution_parser.length

    def _plot_species_mass_fraction_without_discretization(self, species_names, xlabel):
        """
        Plotting mass fraction of BATCH and CSTR reactor
        :param species_names: List of species to be plotted
        :param xlabel: x-axis label
        :return:
        """
        self.nFig = self.nFig + 1
        plt.figure(self.nFig)
        for s in species_names:
            plt.plot(self.x, self.mass_fraction[:, self.cls.gas.species_index(s)])

        plt.ylabel('Mass fraction')
        plt.xlabel(xlabel)
        plt.legend(species_names, loc='best').get_frame().set_linewidth(0.0)

    def _plot_species_mole_fraction_without_discretization(self, species_names, xlabel):
        """
        Plotting mole fraction of BATCH and CSTR reactor
        :param species_names: List of species to be plotted
        :param xlabel: x-axis label
        :return:
        """
        self.nFig = self.nFig + 1
        plt.figure(self.nFig)

        for s in species_names:
            plt.plot(self.x, self.mole_fraction[:, self.cls.gas.species_index(s)])

        plt.ylabel('Mole fraction')
        plt.xlabel(xlabel)
        plt.legend(species_names, loc='best').get_frame().set_linewidth(0.0)

    def _plot_coverage_without_discretization(self, coverage_names, xlabel):
        """
        Plotting coverage of BATCH and CSTR reactor
        :param coverage_names: List of coverage species to be plotted
        :param xlabel: x-axis label
        :return:
        """
        self.nFig = self.nFig + 1
        plt.figure(self.nFig)

        for s in coverage_names:
            plt.plot(self.x, self.coverage[:, self.cls.surf.species_index(s)])

        plt.ylabel('Site fraction')
        plt.xlabel(xlabel)
        plt.legend(coverage_names, loc='best').get_frame().set_linewidth(0.0)

    def _plot_temperature_without_discretization(self, xlabel):
        """
        Plotting temperature of BATCH and CSTR reactor
        :param xlabel: x-axis label
        :return:
        """
        self.nFig = self.nFig + 1
        plt.figure(self.nFig)
        plt.plot(self.x, self.temperature)

        plt.ylabel('Temperature [K]')
        plt.xlabel(xlabel)

    def _plot_species_mass_fraction_pseudo_homogeneous_transient(self, species_names):
        """
        Plotting mass fraction of PSEUDO HOMOGENEOUS reactor TRANSIENT
        :param species_names: List of species to be plotted
        :return:
        """
        for s in species_names:
            self.nFig = self.nFig + 1
            plt.figure(self.nFig)
            legend = list()
            for k, t in enumerate(self.x):
                plt.plot(self.length, self.mass_fraction[k][:, self.cls.gas.species_index(s)], color=self.colors[k])
                legend.append("Time: " + str(round(t, 3)) + " s")

            plt.ylabel('Mass fraction')
            plt.xlabel('Length [m]')
            plt.legend(legend, loc='best').get_frame().set_linewidth(0.0)
            plt.title(s, loc='center')

    def _plot_species_mole_fraction_pseudo_homogeneous_transient(self, species_names):
        """
        Plotting mole fraction of PSEUDO HOMOGENEOUS reactor TRANSIENT
        :param species_names: List of species to be plotted
        :return:
        """
        for s in species_names:
            self.nFig = self.nFig + 1
            plt.figure(self.nFig)
            legend = list()
            for k, t in enumerate(self.x):
                plt.plot(self.length, self.mole_fraction[k][:, self.cls.gas.species_index(s)], color=self.colors[k])
                legend.append("Time: " + str(round(t, 3)) + " s")

            plt.ylabel('Mole fraction')
            plt.xlabel('Length [m]')
            plt.legend(legend, loc='best').get_frame().set_linewidth(0.0)
            plt.title(s, loc='center')

    def _plot_coverage_pseudo_homogeneous_transient_and_heterogeneous(self, coverage_names):
        """
        Plotting coverage of PSEUDO HOMOGENEOUS reactor TRANSIENT and HETEROGENEOUS reactor
        :param coverage_names: List of coverage species to be plotted
        :return:
        """
        for s in coverage_names:
            self.nFig = self.nFig + 1
            plt.figure(self.nFig)
            legend = list()
            for k, t in enumerate(self.x):
                plt.plot(self.length, self.coverage[k][:, self.cls.surf.species_index(s)],
                         color=self.colors[k])
                legend.append("Time: " + str(round(t, 3)) + " s")

            plt.ylabel('Site fraction')
            plt.xlabel('Length [m]')
            plt.legend(legend, loc='best').get_frame().set_linewidth(0.0)
            plt.title(s, loc='center')

    def _plot_temperature_pseudo_homogeneous_transient(self):
        """
        Plotting temperature of PSEUDO HOMOGENEOUS reactor TRANSIENT
        :return:
        """
        self.nFig = self.nFig + 1
        plt.figure(self.nFig)
        legend = list()
        for k, t in enumerate(self.x):
            plt.plot(self.length, self.temperature[k], color=self.colors[k])
            legend.append("Time: " + str(round(t, 3)) + " s")

        plt.ylabel('Temperature [K]')
        plt.xlabel('Length [m]')
        plt.legend(legend, loc='best').get_frame().set_linewidth(0.0)

    def _plot_species_mass_fraction_heterogeneous(self, species_names):
        """
        Plotting mass fraction of HETEROGENEOUS reactor
        :param species_names: List of species to be plotted
        :return:
        """
        for s in species_names:
            self.nFig = self.nFig + 1
            plt.figure(self.nFig)
            legend = list()

            for k, t in enumerate(self.x):
                plt.plot(self.length, self.mass_fraction[k][:, self.cls.gas.species_index(s)], color=self.colors[k])
                legend.append("Time: " + str(round(t, 3)) + " s")

            plt.ylabel('Gas mass fraction')
            plt.xlabel('Length [m]')
            plt.legend(legend, loc='best').get_frame().set_linewidth(0.0)
            plt.title(s, loc='center')

            self.nFig = self.nFig + 1
            plt.figure(self.nFig)
            legend = list()
            for k, t in enumerate(self.x):
                plt.plot(self.length, self.mass_fraction_wall[k][:, self.cls.gas.species_index(s)],
                         color=self.colors[k])
                legend.append("Time: " + str(round(t, 3)) + " s")

            plt.ylabel('Solid mass fraction')
            plt.xlabel('Length [m]')
            plt.legend(legend, loc='best').get_frame().set_linewidth(0.0)
            plt.title(s, loc='center')

    def _plot_species_mole_fraction_heterogeneous(self, species_names):
        """
        Plotting mole fraction of HETEROGENEOUS reactor
        :param species_names: List of species to be plotted
        :return:
        """
        for s in species_names:
            self.nFig = self.nFig + 1
            plt.figure(self.nFig)
            legend = list()
            for k, t in enumerate(self.x):
                plt.plot(self.length, self.mole_fraction[k][:, self.cls.gas.species_index(s)], color=self.colors[k])
                legend.append("Time: " + str(round(t, 3)) + " s")

            plt.ylabel('Gas mole fraction')
            plt.xlabel('Length [m]')
            plt.legend(legend, loc='best').get_frame().set_linewidth(0.0)
            plt.title(s, loc='center')

            self.nFig = self.nFig + 1
            plt.figure(self.nFig)
            legend = list()
            for k, t in enumerate(self.x):
                plt.plot(self.length, self.mole_fraction_wall[k][:, self.cls.gas.species_index(s)],
                         color=self.colors[k])
                legend.append("Time: " + str(round(t, 3)) + " s")

            plt.ylabel('Solid mole fraction')
            plt.xlabel('Length [m]')
            plt.legend(legend, loc='best').get_frame().set_linewidth(0.0)
            plt.title(s, loc='center')

    def _plot_temperature_heterogeneous(self):
        """
        Plotting temperature of HETEROGENEOUS reactor
        :return:
        """
        self.nFig = self.nFig + 1
        plt.figure(self.nFig)
        legend = list()
        for k, t in enumerate(self.x):
            plt.plot(self.length, self.temperature[k], color=self.colors[k])
            legend.append("Time: " + str(round(t, 3)) + " s")

        plt.ylabel('Gas temperature [K]')
        plt.xlabel('Length [m]')
        plt.legend(legend, loc='best').get_frame().set_linewidth(0.0)

        self.nFig = self.nFig + 1
        plt.figure(self.nFig)
        legend = list()
        for k, t in enumerate(self.x):
            plt.plot(self.length, self.temperature_wall[k], color=self.colors[k])
            legend.append("Time: " + str(round(t, 3)) + " s")

        plt.ylabel('Solid temperature [K]')
        plt.xlabel('Length [m]')
        plt.legend(legend, loc='best').get_frame().set_linewidth(0.0)

    def plot_species_mass_fraction(self, species_names):
        """
        Plotting mass fraction
        :param species_names: List of species to be plotted
        :return:
        """
        if self.cls.solution_parser.reactor_type == ReactorType.BATCH or self.cls.solution_parser.reactor_type == ReactorType.CSTR:
            self._plot_species_mass_fraction_without_discretization(species_names, 'Time [s]')
        elif self.cls.solution_parser.reactor_type == ReactorType.PSEUDOHOMOGENEOUSPFR:
            if self.cls.solution_parser.resolution_method == ResolutionMethod.STEADYSTATE:
                self._plot_species_mass_fraction_without_discretization(species_names, 'Length [m]')
            elif self.cls.solution_parser.resolution_method == ResolutionMethod.TRANSIENT:
                self._plot_species_mass_fraction_pseudo_homogeneous_transient(species_names)
        elif self.cls.solution_parser.reactor_type == ReactorType.HETEROGENEOUSPRF:
            self._plot_species_mass_fraction_heterogeneous(species_names)

    def plot_species_mole_fraction(self, species_names):
        """
        Plotting mole fraction
        :param species_names: List of species to be plotted
        :return:
        """
        if self.cls.solution_parser.reactor_type == ReactorType.BATCH or self.cls.solution_parser.reactor_type == ReactorType.CSTR:
            self._plot_species_mole_fraction_without_discretization(species_names, 'Time [s]')
        elif self.cls.solution_parser.reactor_type == ReactorType.PSEUDOHOMOGENEOUSPFR:
            if self.cls.solution_parser.resolution_method == ResolutionMethod.STEADYSTATE:
                self._plot_species_mole_fraction_without_discretization(species_names, 'Length [m]')
            elif self.cls.solution_parser.resolution_method == ResolutionMethod.TRANSIENT:
                self._plot_species_mole_fraction_pseudo_homogeneous_transient(species_names)
        elif self.cls.solution_parser.reactor_type == ReactorType.HETEROGENEOUSPRF:
            self._plot_species_mole_fraction_heterogeneous(species_names)

    def plot_coverage(self, coverage_names):
        """
        Plotting coverage
        :param coverage_names: List of coverage species to be plotted
        :return:
        """
        if self.cls.solution_parser.reactor_type == ReactorType.BATCH or self.cls.solution_parser.reactor_type == ReactorType.CSTR:
            self._plot_coverage_without_discretization(coverage_names, 'Time [s]')
        elif self.cls.solution_parser.reactor_type == ReactorType.PSEUDOHOMOGENEOUSPFR:
            if self.cls.solution_parser.resolution_method == ResolutionMethod.STEADYSTATE:
                self._plot_coverage_without_discretization(coverage_names, 'Length [m]')
            elif self.cls.solution_parser.resolution_method == ResolutionMethod.TRANSIENT:
                self._plot_coverage_pseudo_homogeneous_transient_and_heterogeneous(coverage_names)
        elif self.cls.solution_parser.reactor_type == ReactorType.HETEROGENEOUSPRF:
            self._plot_coverage_pseudo_homogeneous_transient_and_heterogeneous(coverage_names)

    def plot_temperature(self):
        """
        Plotting temperature
        :return:
        """
        if self.cls.solution_parser.reactor_type == ReactorType.BATCH or self.cls.solution_parser.reactor_type == ReactorType.CSTR:
            self._plot_temperature_without_discretization('Time [s]')
        elif self.cls.solution_parser.reactor_type == ReactorType.PSEUDOHOMOGENEOUSPFR:
            if self.cls.solution_parser.resolution_method == ResolutionMethod.STEADYSTATE:
                self._plot_temperature_without_discretization('Length [m]')
            elif self.cls.solution_parser.resolution_method == ResolutionMethod.TRANSIENT:
                self._plot_temperature_pseudo_homogeneous_transient()
        elif self.cls.solution_parser.reactor_type == ReactorType.HETEROGENEOUSPRF:
            self._plot_temperature_heterogeneous()

    @staticmethod
    def show():
        """
        Show plots
        :return:
        """
        plt.show()
