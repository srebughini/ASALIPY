from abc import abstractmethod

import numpy as np

from asali.utils.basic_supporter import BasicSupporter

import pandas as pd


class BasicSaver(BasicSupporter):
    def __init__(self, cls):
        """
        Abstract class that can be used to save results
        :param cls: Object of class to be plotted
        """
        super().__init__(cls)

    @abstractmethod
    def species_mass_fraction_to_df(self, species_names):
        """
        Saving mass fraction
        :param species_names: List of species to be plotted
        :return: Dictionary with sheet name and pd.DataFrame
        """
        pass

    @abstractmethod
    def species_mole_fraction_to_df(self, species_names):
        """
        Saving mass fraction
        :param species_names: List of species to be plotted
        :return: Dictionary with sheet name and pd.DataFrame
        """
        pass

    @abstractmethod
    def coverage_to_df(self, coverage_names):
        """
        Saving coverage
        :param coverage_names: List of coverage species to be plotted
        :return: Dictionary with sheet name and pd.DataFrame
        """
        pass

    @abstractmethod
    def temperature_to_df(self):
        """
        Saving temperature
        :return: Dictionary with sheet name and pd.DataFrame
        """
        pass

    @staticmethod
    def save_to_excel(input_dict, file_path):
        """
        Save results into Excel format
        :param input_dict: Input dictionary with {sheetName: pd.Dataframe} format
        :param file_path: File path where to save the results
        :return:
        """
        with pd.ExcelWriter(file_path) as writer:
            for sheet_name, df in input_dict.items():
                df.to_excel(writer, sheet_name=sheet_name, index=False)

    @staticmethod
    def create_df_from_x_vector_and_y_vector(xv, yv, x_name, y_name):
        """
        Create pd.DataFrame from x and array of y by adding columns names
        :param xv: Independent variable vector
        :param yv: Dependent variable matrix
        :param x_name: Independent variable name
        :param y_name: Dependent variable label
        :return: Dictionary with sheet name and pd.DataFrame
        """

        df = pd.DataFrame({x_name: xv,
                           y_name: yv})

        sheet_name = y_name

        return {sheet_name: df}

    @staticmethod
    def create_df_from_x_vector_and_y_matrix(xv, ym, x_name, y_name, column_names):
        """
        Create pd.DataFrame from x and array of y by adding columns names
        :param xv: Independent variable vector
        :param ym: Dependent variable matrix
        :param x_name: Independent variable name
        :param y_name: Dependent variable label
        :param column_names: Legend list
        :return: Dictionary with sheet name and pd.DataFrame
        """

        columns = [x_name]
        columns.extend(column_names)

        data = np.asarray(ym)
        data = data.T
        data = np.insert(data, 0, xv, 1)

        df = pd.DataFrame(data=data,
                          columns=columns)

        sheet_name = y_name

        return {sheet_name: df}

    def parse_specie_and_coverage_names(self, species_names, coverage_names):
        """
        Parse specie and coverage names
        :param species_names: List of species to be saved
        :param coverage_names: List of coverage to be saved
        :return: species_names, coverage_names
        """
        if species_names is None:
            species_names = self.cls.gas.species_names

        if coverage_names is None:
            coverage_names = self.cls.surf.species_names

        return species_names, coverage_names

    def save_using_mass_fraction(self, file_path, species_names=None, coverage_names=None):
        """
        Saving output to file
        :param file_path: File path where to save the results
        :param species_names: List of species to be saved
        :param coverage_names: List of coverage to be saved
        :return:
        """
        species_names, coverage_names = self.parse_specie_and_coverage_names(species_names, coverage_names)

        output_dict = self.species_mass_fraction_to_df(species_names)
        output_dict.update(self.coverage_to_df(coverage_names))
        output_dict.update(self.temperature_to_df())

        self.save_to_excel(output_dict, file_path)

    def save_using_mole_fraction(self, file_path, species_names=None, coverage_names=None):
        """
        Saving output to file
        :param file_path: File path where to save the results
        :param species_names: List of species to be saved
        :param coverage_names: List of coverage to be saved
        :return:
        """
        species_names, coverage_names = self.parse_specie_and_coverage_names(species_names, coverage_names)

        output_dict = self.species_mole_fraction_to_df(species_names)
        output_dict.update(self.coverage_to_df(coverage_names))
        output_dict.update(self.temperature_to_df())

        self.save_to_excel(output_dict, file_path)
