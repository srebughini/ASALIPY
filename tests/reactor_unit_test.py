from asali.utils.input_parser import ResolutionMethod, ReactorType
from tests.basic_unit_test import BasicUnitTest
from importlib import import_module

import numpy as np

import os


class ReactorUnitTest(BasicUnitTest):
    def __init__(self, cls, file_name, is_local, folder_path="tests/"):
        """
        Derived class for unit test of ASALIPY reactor models
        :param cls: Class to be tested
        :param file_name: Json file with the tests to be performed
        :param folder_path: Folder path where input can be found
        """
        super().__init__(cls=cls, file_name=file_name, folder_path=folder_path)
        self.is_local = is_local

    def _convert_path(self, f, results, results_format):
        """
        Function to convert a .asali file into an array
        :param f: Function to be tested
        :param results: Expected results
        :param results_format: Expected results format
        :return: Expected results, Expected results format
        """
        if results_format == "path":
            if self.cls.solution_parser.resolution_method == ResolutionMethod.STEADYSTATE:
                file_name = "{}.asali".format(f.__name__)
                results_new = np.loadtxt(os.path.join(self.folder_path, self.cls.__class__.__name__, file_name))
                results_format_new = "array"
                return results_new, results_format_new

            if self.cls.solution_parser.resolution_method == ResolutionMethod.TRANSIENT:
                file_name = "{}.asali".format(f.__name__)
                results_new = np.loadtxt(os.path.join(self.folder_path, self.cls.__class__.__name__, file_name))
                results_format_new = "array"
                return results_new, results_format_new

            file_name = "{}.asali".format(f.__name__)
            results_new = np.loadtxt(os.path.join(self.folder_path, self.cls.__class__.__name__, file_name))
            results_format_new = "array"
            return results_new, results_format_new

        return results, results_format

    def check_all(self):
        """
        Function to test all functions of a class
        :return:
        """
        for function in self.function_to_check:
            if function == "example":
                self.run_example()
            else:
                f = getattr(self.cls, function)
                self.check_function(f)

    def run_example(self,
                    atol=1.e-04,
                    rtol=1.e-04):
        """
        Function to run the examples as Unit Test
        :param atol: Comparison absolute tolerance
        :param rtol: Comparison relative tolerance
        :return:
        """

        if self.is_local:
            example_module_format = 'tests.{}.example'

            mod = import_module(example_module_format.format(self.cls.__class__.__name__))
            outputs = mod.main(self.cls.cantera_input_file, self.cls.gas_phase_name, self.cls.surface_phase_name)
            results = np.loadtxt(os.path.join(self.folder_path, self.cls.__class__.__name__, "example.asali"))

            if outputs.shape == results.shape:
                check = np.allclose(outputs, results, atol=atol, rtol=rtol)
            else:
                check = False

            if check:
                self._print_on_screen("example", "OK", "green")
            else:
                self._print_on_screen("example", "ASALI::ERROR", "red")
                self._print_comparison_on_screen(outputs, results, "array")
        else:
            self._print_on_screen("example", "ASALI::WARNING::NOT PERFORMED", "yellow")
