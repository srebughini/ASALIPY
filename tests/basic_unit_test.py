from termcolor import colored

import json
import os
import filecmp
import yaml

import numpy as np


class BasicUnitTest:
    def __init__(self, cls, file_name, folder_path="tests/"):
        """
        Basic class for unit test of ASALIPY
        :param cls: Class to be tested
        :param file_name: Json file with the tests to be performed
        :param folder_path: Folder path where input can be found
        """
        with open(os.path.join(folder_path, file_name)) as json_file:
            solution_dict = json.load(json_file)

        self.folder_path = folder_path
        self.cls = cls
        self.setting_dict = solution_dict[cls.__class__.__name__]
        self.function_to_check = self.setting_dict.keys()
        self.keys_to_ignore_for_yaml = ["input-files", "date"]

    @staticmethod
    def _get_output(f, args, args_format):
        """
        Function to run a general function
        :param f: Function to be run
        :param args: Args of the function
        :param args_format: Args type
        :return: Output of the run function
        """
        if args_format == "none":
            return f()

        if args_format == "tuple":
            return f(*args)

        return f(args)

    @staticmethod
    def _print_comparison_on_screen(outputs, results, results_format):
        """
        Print on screen the comparison between expected and obtained results
        :param outputs: Obtained results
        :param results: Expected results
        :param results_format: Format of the results
        :return:
        """
        msg = 'Output value: {}\nExpected value: {}'
        if results_format == "enum":
            print(msg.format(int(outputs), int(results)))
        else:
            print(msg.format(outputs, results))

    def _print_on_screen(self, f, output, color):
        """
        Print on screen the test function and its results
        :param f: Tested function
        :param output: Test output (OK/NOT OK)
        :param color: Color to be used
        :return:
        """
        if not isinstance(f, str):
            function_name = f.__name__
        else:
            function_name = f

        tested_function = '{}::{}'.format(self.cls.__class__.__name__, function_name)
        msg = 'ASALI::{}{} --> {}'.format(tested_function, ' ' * (85 - len(tested_function)), output)
        print(colored(msg, color))

    def _convert_path(self, f, results, results_format):
        """
        Function to convert a .asali file into an array
        :param f: Function to be tested
        :param results: Expected results
        :param results_format: Expected results format
        :return: Expected results, Expected results format
        """
        if results_format == "path":
            file_name = "{}.asali".format(f.__name__)
            results_new = np.loadtxt(os.path.join(self.folder_path, self.cls.__class__.__name__, file_name))
            results_format_new = "array"
            return results_new, results_format_new

        return results, results_format

    def _check_float(self, f, results, args, args_format):
        """
        Function to check float outputs
        :param f: Function to be tested
        :param results: Expected results
        :param args: Args of the function
        :param args_format: Args type
        :return: True/False output of the test, Function results
        """
        outputs = self._get_output(f, args, args_format)
        return np.fabs(outputs - results) < 1.e-12, outputs, results

    def _check_array(self, f, results, args, args_format, atol=1.e-02, rtol=1.e-02):
        """
        Function to check array outputs
        :param f: Function to be tested
        :param results: Expected results
        :param args: Args of the function
        :param args_format: Args type
        :param atol: Comparison absolute tolerance
        :param rtol: Comparison relative tolerance
        :return: True/False output of the test, Function results
        """
        outputs = self._get_output(f, args, args_format)

        outputs_as_array = np.asarray(outputs)
        results_as_array = np.asarray(results)

        if outputs_as_array.shape == results_as_array.shape:
            return np.allclose(outputs, results, atol=atol, rtol=rtol), outputs, results

        if outputs_as_array.size == results_as_array.size:
            return np.allclose(outputs.flatten(), results.flatten(), atol=atol,
                               rtol=rtol), outputs.flatten(), results.flatten()

        return False, outputs, results

    def _check_enum(self, f, results, args, args_format):
        """
        Function to check enum outputs
        :param f: Function to be tested
        :param results: Expected results
        :param args: Args of the function
        :param args_format: Args type
        :return: True/False output of the test, Function results
        """
        outputs = self._get_output(f, args, args_format)
        return int(outputs) == results, outputs, results

    def _check_others(self, f, results, args, args_format):
        """
        Function to check not float, array or enum outputs
        :param f: Function to be tested
        :param results: Expected results
        :param args: Args of the function
        :param args_format: Args type
        :return: True/False output of the test, Function results
        """
        outputs = self._get_output(f, args, args_format)
        return outputs == results, outputs, results

    def _check_files(self, f, results, args, args_format):
        """
        Function to check file outputs
        :param f: Function to be tested
        :param results: Expected results
        :param args: Args of the function
        :param args_format: Args type
        :return: True/False output of the test, Function results
        """
        outputs = self._get_output(f, args, args_format)

        output_file_as_list = [line.strip() for line in open(outputs)]
        results_file_as_list = [line.strip() for line in open(results)]

        if len(output_file_as_list) == len(results_file_as_list):
            for i, output_line in enumerate(output_file_as_list):
                if output_line != results_file_as_list[i]:
                    return False, '\n'.join(output_file_as_list), '\n'.join(results_file_as_list)

            return True, '\n'.join(output_file_as_list), '\n'.join(results_file_as_list)

        return False, '\n'.join(output_file_as_list), '\n'.join(results_file_as_list)

    def _check_yaml(self, f, results, args, args_format):
        """
        Function to check YAML file outputs
        :param f: Function to be tested
        :param results: Expected results
        :param args: Args of the function
        :param args_format: Args type
        :return: True/False output of the test, Function results
        """
        outputs = self._get_output(f, args, args_format)

        with open(outputs, 'r') as stream:
            outputs_as_yaml = yaml.safe_load(stream)

        with open(results, 'r') as stream:
            results_as_yaml = yaml.safe_load(stream)

        for key_to_remove in self.keys_to_ignore_for_yaml:
            if key_to_remove in results_as_yaml.keys():
                del results_as_yaml[key_to_remove]

            if key_to_remove in outputs_as_yaml.keys():
                del outputs_as_yaml[key_to_remove]

        return outputs_as_yaml == results_as_yaml, outputs, results

    def check_function(self, f):
        """
        Function to test a single function
        :param f: Function to be tested
        :return:
        """
        args = self.setting_dict[f.__name__]["input"]["value"]
        results = self.setting_dict[f.__name__]["output"]["value"]
        args_format = self.setting_dict[f.__name__]["input"]["format"]
        results_format = self.setting_dict[f.__name__]["output"]["format"]

        results, results_format = self._convert_path(f, results, results_format)

        if results_format == "float":
            check, outputs, results = self._check_float(f, results, args, args_format)
        elif results_format == "array":
            check, outputs, results = self._check_array(f, results, args, args_format)
        elif results_format == "enum":
            check, outputs, results = self._check_enum(f, results, args, args_format)
        elif results_format == "file":
            check, outputs, results = self._check_files(f, results, args, args_format)
        elif results_format == "yaml":
            check, outputs, results = self._check_yaml(f, results, args, args_format)
        else:
            check, outputs, results = self._check_others(f, results, args, args_format)

        if check:
            self._print_on_screen(f, "OK", "green")
        else:
            self._print_on_screen(f, "ASALI::ERROR", "red")
            self._print_comparison_on_screen(outputs, results, results_format)

    def check_function_from_string(self, f_as_string):
        """
        Function to test a single function from its name
        :param f_as_string: Function name as str
        :return:
        """
        f = getattr(self.cls, f_as_string)
        self.check_function(f)

    def check_all(self):
        """
        Function to test all functions of a class
        :return:
        """
        for function in self.function_to_check:
            f = getattr(self.cls, function)
            self.check_function(f)
