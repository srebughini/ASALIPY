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

    def plot_species_mass_fraction(self, species_names):
        """
        Plotting mass fraction
        :param species_names: List of species to be plotted
        :return:
        """
        if self.cls.solution_parser.reactor_type == ReactorType.BATCH or self.cls.solution_parser.reactor_type == ReactorType.CSTR:
            ym = [self.mass_fraction[:, self.cls.gas.species_index(s)] for s in species_names]
            self.plot_x_vector_and_y_matrix(plt,
                                            self.x,
                                            ym,
                                            'Time [s]',
                                            'Mass fraction',
                                            species_names)

        if self.cls.solution_parser.reactor_type == ReactorType.PSEUDOHOMOGENEOUSPFR:
            if self.cls.solution_parser.resolution_method == ResolutionMethod.STEADYSTATE:
                ym = [self.mass_fraction[:, self.cls.gas.species_index(s)] for s in species_names]
                self.plot_x_vector_and_y_matrix(plt,
                                                self.x,
                                                ym,
                                                'Length [m]',
                                                'Mass fraction',
                                                species_names)

            if self.cls.solution_parser.resolution_method == ResolutionMethod.TRANSIENT:
                legend = ["Time: " + str(round(t, 3)) + " s" for t in self.x]
                colors = [self.colors[k] for k, _ in enumerate(self.x)]
                for s in species_names:
                    ym = [self.mass_fraction[k][:, self.cls.gas.species_index(s)] for k, _ in enumerate(self.x)]
                    self.plot_x_vector_and_y_matrix(plt,
                                                    self.length,
                                                    ym,
                                                    'Length [m]',
                                                    'Mass fraction',
                                                    legend,
                                                    title=s,
                                                    colors=colors)

        if self.cls.solution_parser.reactor_type == ReactorType.HETEROGENEOUSPRF:
            if self.cls.solution_parser.resolution_method == ResolutionMethod.STEADYSTATE:
                ym = [self.mass_fraction[:, self.cls.gas.species_index(s)] for s in species_names]
                self.plot_x_vector_and_y_matrix(plt,
                                                self.x,
                                                ym,
                                                'Length [m]',
                                                'Gas mass fraction',
                                                species_names)

                ym = [self.mass_fraction_wall[:, self.cls.gas.species_index(s)] for s in species_names]
                self.plot_x_vector_and_y_matrix(plt,
                                                self.x,
                                                ym,
                                                'Length [m]',
                                                'Solid mass fraction',
                                                species_names)

            if self.cls.solution_parser.resolution_method == ResolutionMethod.TRANSIENT:
                legend = ["Time: " + str(round(t, 3)) + " s" for t in self.x]
                colors = [self.colors[k] for k, _ in enumerate(self.x)]
                for s in species_names:
                    ym = [self.mass_fraction[k][:, self.cls.gas.species_index(s)] for k, _ in enumerate(self.x)]
                    self.plot_x_vector_and_y_matrix(plt,
                                                    self.length,
                                                    ym,
                                                    'Length [m]',
                                                    'Gas mass fraction',
                                                    legend,
                                                    title=s,
                                                    colors=colors)

                for s in species_names:
                    ym = [self.mass_fraction_wall[k][:, self.cls.gas.species_index(s)] for k, _ in enumerate(self.x)]
                    colors = [self.colors[k] for k, _ in enumerate(self.x)]
                    self.plot_x_vector_and_y_matrix(plt,
                                                    self.length,
                                                    ym,
                                                    'Length [m]',
                                                    'Solid mass fraction',
                                                    legend,
                                                    title=s,
                                                    colors=colors)

    def plot_species_mole_fraction(self, species_names):
        """
        Plotting mole fraction
        :param species_names: List of species to be plotted
        :return:
        """
        if self.cls.solution_parser.reactor_type == ReactorType.BATCH or self.cls.solution_parser.reactor_type == ReactorType.CSTR:
            ym = [self.mole_fraction[:, self.cls.gas.species_index(s)] for s in species_names]
            self.plot_x_vector_and_y_matrix(plt,
                                            self.x,
                                            ym,
                                            'Time [s]',
                                            'Mole fraction',
                                            species_names)

        if self.cls.solution_parser.reactor_type == ReactorType.PSEUDOHOMOGENEOUSPFR:
            if self.cls.solution_parser.resolution_method == ResolutionMethod.STEADYSTATE:
                ym = [self.mole_fraction[:, self.cls.gas.species_index(s)] for s in species_names]
                self.plot_x_vector_and_y_matrix(plt,
                                                self.x,
                                                ym,
                                                'Length [m]',
                                                'Mole fraction',
                                                species_names)

            if self.cls.solution_parser.resolution_method == ResolutionMethod.TRANSIENT:
                legend = ["Time: " + str(round(t, 3)) + " s" for t in self.x]
                colors = [self.colors[k] for k, _ in enumerate(self.x)]
                for s in species_names:
                    ym = [self.mole_fraction[k][:, self.cls.gas.species_index(s)] for k, _ in enumerate(self.x)]
                    self.plot_x_vector_and_y_matrix(plt,
                                                    self.length,
                                                    ym,
                                                    'Length [m]',
                                                    'Mole fraction',
                                                    legend,
                                                    title=s,
                                                    colors=colors)

        if self.cls.solution_parser.reactor_type == ReactorType.HETEROGENEOUSPRF:
            if self.cls.solution_parser.resolution_method == ResolutionMethod.STEADYSTATE:
                ym = [self.mole_fraction[:, self.cls.gas.species_index(s)] for s in species_names]
                self.plot_x_vector_and_y_matrix(plt,
                                                self.x,
                                                ym,
                                                'Length [m]',
                                                'Gas mole fraction',
                                                species_names)

                ym = [self.mole_fraction_wall[:, self.cls.gas.species_index(s)] for s in species_names]
                self.plot_x_vector_and_y_matrix(plt,
                                                self.x,
                                                ym,
                                                'Length [m]',
                                                'Solid mole fraction',
                                                species_names)

            if self.cls.solution_parser.resolution_method == ResolutionMethod.TRANSIENT:
                legend = ["Time: " + str(round(t, 3)) + " s" for t in self.x]
                colors = [self.colors[k] for k, _ in enumerate(self.x)]
                for s in species_names:
                    ym = [self.mole_fraction[k][:, self.cls.gas.species_index(s)] for k, _ in enumerate(self.x)]
                    self.plot_x_vector_and_y_matrix(plt,
                                                    self.length,
                                                    ym,
                                                    'Length [m]',
                                                    'Gas mole fraction',
                                                    legend,
                                                    title=s,
                                                    colors=colors)

                for s in species_names:
                    ym = [self.mole_fraction_wall[k][:, self.cls.gas.species_index(s)] for k, _ in enumerate(self.x)]
                    self.plot_x_vector_and_y_matrix(plt,
                                                    self.length,
                                                    ym,
                                                    'Length [m]',
                                                    'Solid mole fraction',
                                                    legend,
                                                    title=s,
                                                    colors=colors)

    def plot_coverage(self, coverage_names):
        """
        Plotting coverage
        :param coverage_names: List of coverage species to be plotted
        :return:
        """
        if self.cls.solution_parser.reactor_type == ReactorType.BATCH or self.cls.solution_parser.reactor_type == ReactorType.CSTR:
            ym = [self.coverage[:, self.cls.surf.species_index(s)] for s in coverage_names]
            self.plot_x_vector_and_y_matrix(plt,
                                            self.x,
                                            ym,
                                            'Time [s]',
                                            'Site fraction',
                                            coverage_names)

        if self.cls.solution_parser.reactor_type == ReactorType.PSEUDOHOMOGENEOUSPFR or self.cls.solution_parser.reactor_type == ReactorType.HETEROGENEOUSPRF:
            if self.cls.solution_parser.resolution_method == ResolutionMethod.STEADYSTATE:
                ym = [self.coverage[:, self.cls.surf.species_index(s)] for s in coverage_names]
                self.plot_x_vector_and_y_matrix(plt,
                                                self.x,
                                                ym,
                                                'Length [m]',
                                                'Site fraction',
                                                coverage_names)

            if self.cls.solution_parser.resolution_method == ResolutionMethod.TRANSIENT:
                legend = ["Time: " + str(round(t, 3)) + " s" for t in self.x]
                colors = [self.colors[k] for k, _ in enumerate(self.x)]
                for s in coverage_names:
                    ym = [self.coverage[k][:, self.cls.surf.species_index(s)] for k, _ in enumerate(self.x)]
                    self.plot_x_vector_and_y_matrix(plt,
                                                    self.length,
                                                    ym,
                                                    'Length [m]',
                                                    'Site fraction',
                                                    legend,
                                                    title=s,
                                                    colors=colors)

    def plot_temperature(self):
        """
        Plotting temperature
        :return:
        """
        if self.cls.solution_parser.reactor_type == ReactorType.BATCH or self.cls.solution_parser.reactor_type == ReactorType.CSTR:
            self.plot_x_vector_and_y_vector(plt,
                                            self.x,
                                            self.temperature,
                                            'Time [s]',
                                            'Temperature [K]')

        if self.cls.solution_parser.reactor_type == ReactorType.PSEUDOHOMOGENEOUSPFR:
            if self.cls.solution_parser.resolution_method == ResolutionMethod.STEADYSTATE:
                self.plot_x_vector_and_y_vector(plt,
                                                self.x,
                                                self.temperature,
                                                'Length [m]',
                                                'Temperature [K]')

            if self.cls.solution_parser.resolution_method == ResolutionMethod.TRANSIENT:
                legend = ["Time: " + str(round(t, 3)) + " s" for t in self.x]
                colors = [self.colors[k] for k, _ in enumerate(self.x)]
                ym = [self.temperature[k] for k, _ in enumerate(self.x)]
                self.plot_x_vector_and_y_matrix(plt,
                                                self.length,
                                                ym,
                                                'Length [m]',
                                                'Temperature [K]',
                                                legend,
                                                colors=colors)

        if self.cls.solution_parser.reactor_type == ReactorType.HETEROGENEOUSPRF:
            if self.cls.solution_parser.resolution_method == ResolutionMethod.STEADYSTATE:
                self.plot_x_vector_and_y_vector(plt,
                                                self.x,
                                                self.temperature,
                                                'Length [m]',
                                                'Gas temperature [K]')

                self.plot_x_vector_and_y_vector(plt,
                                                self.x,
                                                self.temperature_wall,
                                                'Length [m]',
                                                'Solid temperature [K]')

            if self.cls.solution_parser.resolution_method == ResolutionMethod.TRANSIENT:
                legend = ["Time: " + str(round(t, 3)) + " s" for t in self.x]
                colors = [self.colors[k] for k, _ in enumerate(self.x)]
                ym = [self.temperature[k] for k, _ in enumerate(self.x)]
                self.plot_x_vector_and_y_matrix(plt,
                                                self.length,
                                                ym,
                                                'Length [m]',
                                                'Gas temperature [K]',
                                                legend,
                                                colors=colors)

                ym = [self.temperature_wall[k] for k, _ in enumerate(self.x)]
                self.plot_x_vector_and_y_matrix(plt,
                                                self.length,
                                                ym,
                                                'Length [m]',
                                                'Solid temperature [K]',
                                                legend,
                                                colors=colors)

    @staticmethod
    def show():
        """
        Show plots
        :return:
        """
        plt.show()
