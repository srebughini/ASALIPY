from asali.plotters.basic import BasicPlotter


class BatchAndCstrPlotter(BasicPlotter):
    def __init__(self, cls, colormap):
        """
        Class to plot BATCH and CSTR reactor results
        :param cls: Object of class to be plotted
        :param colormap: String representing the color map
        """
        super().__init__(cls, colormap)
        self.mass_fraction = self.cls.solution_parser.get_mass_fraction()
        self.mole_fraction = self.cls.solution_parser.get_mole_fraction()
        self.coverage = self.cls.solution_parser.get_coverage()
        self.temperature = self.cls.solution_parser.get_temperature()
        self.x = self.cls.solution_parser.x

    def plot_species_mass_fraction(self, plt, species_names):
        """
        Plotting mass fraction
        :param plt: Matplotlib object
        :param species_names: List of species to be plotted
        :return:
        """
        ym = [self.mass_fraction[:, self.cls.gas.species_index(s)] for s in species_names]
        self.plot_x_vector_and_y_matrix(plt,
                                        self.x,
                                        ym,
                                        self._variables_name.time,
                                        self._variables_name.gas_mass_fraction,
                                        species_names)

    def plot_species_mole_fraction(self, plt, species_names):
        """
        Plotting mole fraction
        :param plt: Matplotlib object
        :param species_names: List of species to be plotted
        :return:
        """
        ym = [self.mole_fraction[:, self.cls.gas.species_index(s)] for s in species_names]
        self.plot_x_vector_and_y_matrix(plt,
                                        self.x,
                                        ym,
                                        self._variables_name.time,
                                        self._variables_name.gas_mole_fraction,
                                        species_names)

    def plot_coverage(self, plt, coverage_names):
        """
        Plotting coverage
        :param plt: Matplotlib object
        :param coverage_names: List of coverage species to be plotted
        :return:
        """
        ym = [self.coverage[:, self.cls.surf.species_index(s)] for s in coverage_names]
        self.plot_x_vector_and_y_matrix(plt,
                                        self.x,
                                        ym,
                                        self._variables_name.time,
                                        self._variables_name.coverage,
                                        coverage_names)

    def plot_temperature(self, plt):
        """
        Plotting temperature
        :param plt: Matplotlib object
        :return:
        """
        self.plot_x_vector_and_y_vector(plt,
                                        self.x,
                                        self.temperature,
                                        self._variables_name.time,
                                        self._variables_name.gas_temperature)
