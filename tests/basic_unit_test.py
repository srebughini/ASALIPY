from termcolor import colored

import json
import os

import numpy as np


class BasicUnitTest:
    def __init__(self, cls, file_name, folder_path="tests/"):
        with open(os.path.join(folder_path, file_name)) as json_file:
            solution_dict = json.load(json_file)

        self.folder_path = folder_path
        self.cls = cls
        self.setting_dict = solution_dict[cls.__class__.__name__]
        self.function_to_check = self.setting_dict.keys()

    @staticmethod
    def _get_output(f, args, args_format):
        if args_format == "none":
            return f()

        if args_format == "tuple":
            return f(*args)

        return f(args)

    @staticmethod
    def _print_comparison_on_screen(outputs, results, results_format):
        msg = 'Output value: {}\nExpected value: {}'
        if results_format == "enum":
            print(msg.format(int(outputs), int(results)))
        else:
            print(msg.format(outputs, results))

    def _print_on_screen(self, f, output, color):
        tested_function = '{}::{}'.format(self.cls.__class__.__name__, f.__name__)
        msg = 'ASALI::{}{} --> {}'.format(tested_function, ' ' * (60 - len(tested_function)), output)
        print(colored(msg, color))

    def _convert_path(self, f, results, results_format):
        if results_format == "path":
            file_name = "{}.asali".format(f.__name__)
            results_new = np.loadtxt(os.path.join(self.folder_path, self.cls.__class__.__name__, file_name))
            results_format_new = "array"
            return results_new, results_format_new

        return results, results_format

    def _check_float(self, f, results, args, args_format):
        outputs = self._get_output(f, args, args_format)
        return np.fabs(outputs - results) < 1.e-12, outputs

    def _check_array(self, f, results, args, args_format, atol=1.e-04, rtol=1.e-04):
        outputs = self._get_output(f, args, args_format)
        return np.allclose(outputs, results, atol=atol, rtol=rtol), outputs

    def _check_enum(self, f, results, args, args_format):
        outputs = self._get_output(f, args, args_format)
        return int(outputs) == results, outputs

    def _check_others(self, f, results, args, args_format):
        outputs = self._get_output(f, args, args_format)
        return outputs == results, outputs

    def check(self, f):
        args = self.setting_dict[f.__name__]["input"]["value"]
        results = self.setting_dict[f.__name__]["output"]["value"]
        args_format = self.setting_dict[f.__name__]["input"]["format"]
        results_format = self.setting_dict[f.__name__]["output"]["format"]

        results, results_format = self._convert_path(f, results, results_format)

        if results_format == "float":
            check, outputs = self._check_float(f, results, args, args_format)
        elif results_format == "array":
            check, outputs = self._check_array(f, results, args, args_format)
        elif results_format == "enum":
            check, outputs = self._check_enum(f, results, args, args_format)
        else:
            check, outputs = self._check_others(f, results, args, args_format)

        if check:
            self._print_on_screen(f, "OK", "green")
        else:
            self._print_on_screen(f, "ASALI::ERROR", "red")
            self._print_comparison_on_screen(outputs, results, results_format)

    def check_function(self, f_as_string):
        f = getattr(self.cls, f_as_string)
        self.check(f)

    def check_all(self):
        for function in self.function_to_check:
            f = getattr(self.cls, function)
            self.check(f)
