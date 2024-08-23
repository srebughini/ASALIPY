from asali.plotters.basic import BasicPlotter
from asali.utils.input_parser import ResolutionMethod


class Heterogeneous1DReactorPlotter(BasicPlotter):
    def __init__(self, cls, colormap):
        """
        Class to plot BATCH reactor results
        :param cls: Object of class to be plotted
        :param colormap: String representing the color map
        """
        super().__init__(cls, colormap)
        self.mass_fraction = self.cls.solution_parser.get_mass_fraction()
        self.mass_fraction_wall = self.cls.solution_parser.get_mass_fraction_wall()
        self.mole_fraction = self.cls.solution_parser.get_mole_fraction()
        self.mole_fraction_wall = self.cls.solution_parser.get_mole_fraction_wall()
        self.coverage = self.cls.solution_parser.get_coverage()
        self.temperature = self.cls.solution_parser.get_temperature()
        self.temperature_wall = self.cls.solution_parser.get_temperature_wall()
        self.x = self.cls.solution_parser.x
        self.length = self.cls.solution_parser.length

    def plot_species_mass_fraction(self, plt, species_names):
        """
        Plotting mass fraction
        :param plt: Matplotlib object
        :param species_names: List of species to be plotted
        :return:
        """
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

    def plot_species_mole_fraction(self, plt, species_names):
        """
        Plotting mole fraction
        :param plt: Matplotlib object
        :param species_names: List of species to be plotted
        :return:
        """
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

    def plot_coverage(self, plt, coverage_names):
        """
        Plotting coverage
        :param plt: Matplotlib object
        :param coverage_names: List of coverage species to be plotted
        :return:
        """
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

    def plot_temperature(self, plt):
        """
        Plotting temperature
        :param plt: Matplotlib object
        :return:
        """
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
