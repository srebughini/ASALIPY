from asali.savers.basic import BasicSaver


class BatchAndCstrSaver(BasicSaver):
    def __init__(self, cls):
        """
        Class to plot BATCH reactor results
        :param cls: Object of class to be plotted
        """
        super().__init__(cls)
        self.mass_fraction = self.cls.solution_parser.get_mass_fraction()
        self.mole_fraction = self.cls.solution_parser.get_mole_fraction()
        self.coverage = self.cls.solution_parser.get_coverage()
        self.temperature = self.cls.solution_parser.get_temperature()
        self.x = self.cls.solution_parser.x

    def species_mass_fraction_to_df(self, species_names):
        """
        Saving mass fraction
        :param species_names: List of species to be plotted
        :return: Dictionary with sheet name and pd.DataFrame
        """
        ym = [self.mass_fraction[:, self.cls.gas.species_index(s)] for s in species_names]
        return self.create_df_from_x_vector_and_y_matrix(self.x,
                                                         ym,
                                                         self._variables_name.time,
                                                         self._variables_name.gas_mass_fraction,
                                                         column_names=species_names)

    def species_mole_fraction_to_df(self, species_names):
        """
        Saving mass fraction
        :param species_names: List of species to be plotted
        :return: Dictionary with sheet name and pd.DataFrame
        """
        ym = [self.mole_fraction[:, self.cls.gas.species_index(s)] for s in species_names]
        return self.create_df_from_x_vector_and_y_matrix(self.x,
                                                         ym,
                                                         self._variables_name.time,
                                                         self._variables_name.gas_mole_fraction,
                                                         column_names=species_names)

    def coverage_to_df(self, coverage_names):
        """
        Saving coverage
        :param coverage_names: List of coverage species to be plotted
        :return: Dictionary with sheet name and pd.DataFrame
        """
        ym = [self.coverage[:, self.cls.surf.species_index(s)] for s in coverage_names]
        return self.create_df_from_x_vector_and_y_matrix(self.x,
                                                         ym,
                                                         self._variables_name.time,
                                                         self._variables_name.coverage,
                                                         column_names=coverage_names)

    def temperature_to_df(self):
        """
        Saving temperature
        :return: Dictionary with sheet name and pd.DataFrame
        """
        return self.create_df_from_x_vector_and_y_vector(self.x,
                                                         self.temperature,
                                                         self._variables_name.time,
                                                         self._variables_name.gas_temperature)
