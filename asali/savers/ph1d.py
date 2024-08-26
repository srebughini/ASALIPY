from asali.savers.basic import BasicSaver
from asali.utils.input_parser import ResolutionMethod


class PseudoHomogeneous1DReactorSaver(BasicSaver):
    def __init__(self, cls):
        """
        Class to save PSEUDO HOMOGENOUS 1D REACTOR results
        :param cls: Object of class to be plotted
        """
        super().__init__(cls)
        self.mass_fraction = self.cls.solution_parser.get_mass_fraction()
        self.mole_fraction = self.cls.solution_parser.get_mole_fraction()
        self.coverage = self.cls.solution_parser.get_coverage()
        self.temperature = self.cls.solution_parser.get_temperature()
        self.x = self.cls.solution_parser.x
        self.length = self.cls.solution_parser.length

    def species_mass_fraction_to_df(self, species_names):
        """
        Saving mass fraction
        :param species_names: List of species to be plotted
        :return: Dictionary with sheet name and pd.DataFrame
        """
        if self.cls.solution_parser.resolution_method == ResolutionMethod.STEADYSTATE:
            ym = [self.mass_fraction[:, self.cls.gas.species_index(s)] for s in species_names]
            return self.create_df_from_x_vector_and_y_matrix(self.x,
                                                             ym,
                                                             self._variables_name.length,
                                                             self._variables_name.gas_mass_fraction,
                                                             species_names)

        if self.cls.solution_parser.resolution_method == ResolutionMethod.TRANSIENT:
            output_dict = {}
            column_names = [f"{self._variables_name.time}:{str(round(t, 3))}" for t in self.x]
            for s in species_names:
                ym = [self.mass_fraction[k][:, self.cls.gas.species_index(s)] for k, _ in enumerate(self.x)]
                sheet_name = f"{self._variables_name.gas_mass_fraction} {s}"
                output_dict.update(self.create_df_from_x_vector_and_y_matrix(self.length,
                                                                             ym,
                                                                             self._variables_name.length,
                                                                             sheet_name,
                                                                             column_names))
            return output_dict

    def species_mole_fraction_to_df(self, species_names):
        """
        Saving mole fraction
        :param species_names: List of species to be plotted
        :return: Dictionary with sheet name and pd.DataFrame
        """
        if self.cls.solution_parser.resolution_method == ResolutionMethod.STEADYSTATE:
            ym = [self.mole_fraction[:, self.cls.gas.species_index(s)] for s in species_names]
            return self.create_df_from_x_vector_and_y_matrix(self.x,
                                                             ym,
                                                             self._variables_name.length,
                                                             self._variables_name.gas_mole_fraction,
                                                             species_names)

        if self.cls.solution_parser.resolution_method == ResolutionMethod.TRANSIENT:
            output_dict = {}
            column_names = [f"{self._variables_name.time}:{str(round(t, 3))}" for t in self.x]
            for s in species_names:
                ym = [self.mole_fraction[k][:, self.cls.gas.species_index(s)] for k, _ in enumerate(self.x)]
                sheet_name = f"{self._variables_name.gas_mole_fraction} {s}"
                output_dict.update(self.create_df_from_x_vector_and_y_matrix(self.length,
                                                                             ym,
                                                                             self._variables_name.length,
                                                                             sheet_name,
                                                                             column_names))
            return output_dict

    def coverage_to_df(self, coverage_names):
        """
        Saving coverage
        :param coverage_names: List of coverage species to be plotted
        :return: Dictionary with sheet name and pd.DataFrame
        """
        if self.cls.solution_parser.resolution_method == ResolutionMethod.STEADYSTATE:
            ym = [self.coverage[:, self.cls.surf.species_index(s)] for s in coverage_names]
            return self.create_df_from_x_vector_and_y_matrix(self.x,
                                                             ym,
                                                             self._variables_name.length,
                                                             self._variables_name.coverage,
                                                             coverage_names)

        if self.cls.solution_parser.resolution_method == ResolutionMethod.TRANSIENT:
            output_dict = {}
            column_names = [f"{self._variables_name.time}:{str(round(t, 3))}" for t in self.x]
            for s in coverage_names:
                ym = [self.coverage[k][:, self.cls.surf.species_index(s)] for k, _ in enumerate(self.x)]
                sheet_name = f"{self._variables_name.coverage} {s}"
                output_dict.update(self.create_df_from_x_vector_and_y_matrix(self.length,
                                                                             ym,
                                                                             self._variables_name.length,
                                                                             sheet_name,
                                                                             column_names))
            return output_dict

    def temperature_to_df(self):
        """
        Saving temperature
        :return: Dictionary with sheet name and pd.DataFrame
        """
        if self.cls.solution_parser.resolution_method == ResolutionMethod.STEADYSTATE:
            return self.create_df_from_x_vector_and_y_vector(self.x,
                                                             self.temperature,
                                                             self._variables_name.length,
                                                             self._variables_name.gas_temperature)

        if self.cls.solution_parser.resolution_method == ResolutionMethod.TRANSIENT:
            column_names = [f"{self._variables_name.time}:{str(round(t, 3))}" for t in self.x]
            ym = [self.temperature[k] for k, _ in enumerate(self.x)]
            return self.create_df_from_x_vector_and_y_matrix(self.length,
                                                             ym,
                                                             self._variables_name.length,
                                                             self._variables_name.gas_temperature,
                                                             column_names)
