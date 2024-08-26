from asali.plotters.basic import BasicPlotter
from asali.utils.input_parser import ResolutionMethod


class PseudoHomogeneous1DReactorPlotter(BasicPlotter):
    def __init__(self, cls, colormap):
        """
        Class to plot PSEUDO HOMOGENOUS 1D REACTOR results
        :param cls: Object of class to be plotted
        :param colormap: String representing the color map
        """
        super().__init__(cls, colormap)
        self.mass_fraction = self.cls.solution_parser.get_mass_fraction()
        self.mole_fraction = self.cls.solution_parser.get_mole_fraction()
        self.coverage = self.cls.solution_parser.get_coverage()
        self.temperature = self.cls.solution_parser.get_temperature()
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
                                            self._variables_name.length,
                                            self._variables_name.gas_mass_fraction,
                                            species_names)

        if self.cls.solution_parser.resolution_method == ResolutionMethod.TRANSIENT:
            legend = [f"{self._variables_name.time}:{str(round(t, 3))}" for t in self.x]
            colors = [self.colors[k] for k, _ in enumerate(self.x)]
            for s in species_names:
                ym = [self.mass_fraction[k][:, self.cls.gas.species_index(s)] for k, _ in enumerate(self.x)]
                self.plot_x_vector_and_y_matrix(plt,
                                                self.length,
                                                ym,
                                                self._variables_name.length,
                                                self._variables_name.gas_mass_fraction,
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
                                            self._variables_name.length,
                                            self._variables_name.gas_mole_fraction,
                                            species_names)

        if self.cls.solution_parser.resolution_method == ResolutionMethod.TRANSIENT:
            legend = [f"{self._variables_name.time}:{str(round(t, 3))}" for t in self.x]
            colors = [self.colors[k] for k, _ in enumerate(self.x)]
            for s in species_names:
                ym = [self.mole_fraction[k][:, self.cls.gas.species_index(s)] for k, _ in enumerate(self.x)]
                self.plot_x_vector_and_y_matrix(plt,
                                                self.length,
                                                ym,
                                                self._variables_name.length,
                                                self._variables_name.gas_mole_fraction,
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
                                            self._variables_name.length,
                                            self._variables_name.coverage,
                                            coverage_names)

        if self.cls.solution_parser.resolution_method == ResolutionMethod.TRANSIENT:
            legend = [f"{self._variables_name.time}:{str(round(t, 3))}" for t in self.x]
            colors = [self.colors[k] for k, _ in enumerate(self.x)]
            for s in coverage_names:
                ym = [self.coverage[k][:, self.cls.surf.species_index(s)] for k, _ in enumerate(self.x)]
                self.plot_x_vector_and_y_matrix(plt,
                                                self.length,
                                                ym,
                                                self._variables_name.length,
                                                self._variables_name.coverage,
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
                                            self._variables_name.length,
                                            self._variables_name.gas_temperature)

        if self.cls.solution_parser.resolution_method == ResolutionMethod.TRANSIENT:
            legend = [f"{self._variables_name.time}:{str(round(t, 3))}" for t in self.x]
            colors = [self.colors[k] for k, _ in enumerate(self.x)]
            ym = [self.temperature[k] for k, _ in enumerate(self.x)]
            self.plot_x_vector_and_y_matrix(plt,
                                            self.length,
                                            ym,
                                            self._variables_name.length,
                                            self._variables_name.gas_temperature,
                                            legend,
                                            colors=colors)
