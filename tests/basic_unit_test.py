import json
import os
import sys
import traceback

import numpy as np


class BasicUnitTestException(Exception):
    def __init__(self, message):
        self.message = message
        super().__init__(self.message)

    @staticmethod
    def output_value(f, args, args_format, results_format):
        msg = 'Output value: {}'
        if args_format == "none":
            if results_format == "enum":
                return msg.format(int(f()))
            return msg.format(f())

        if args_format == "tuple":
            if results_format == "enum":
                return msg.format(int(f(*args)))
            return msg.format(f(*args))

        if results_format == "enum":
            return msg.format(int(f(args)))

        return msg.format(f(args))

    @staticmethod
    def expected_value(results, results_format):
        msg = 'Expected value: {}'
        if results_format == "enum":
            return msg.format(int(results))

        return msg.format(results)


class BasicUnitTest:
    def __init__(self, cls, file_name, folder_path="tests/"):
        with open(os.path.join(folder_path, file_name)) as json_file:
            solution_dict = json.load(json_file)

        self.folder_path = folder_path
        self.cls = cls
        self.setting_dict = solution_dict[cls.__class__.__name__]
        self.function_to_check = self.setting_dict.keys()

    def _print_on_screen(self, f, output):
        tested_function = '{}::{}'.format(self.cls.__class__.__name__, f.__name__)
        msg = 'ASALI::{}{} --> {}'.format(tested_function, ' ' * (60 - len(tested_function)), output)
        print(msg)

    def _convert_path(self, f, results, results_format):
        if results_format == "path":
            file_name = "{}.asali".format(f.__name__)
            results_new = np.loadtxt(os.path.join(self.folder_path, self.cls.__class__.__name__, file_name))
            results_format_new = "array"
            return results_new, results_format_new

        return results, results_format

    @staticmethod
    def _check_float(f, results, args, args_format):
        if args_format == "none":
            if not np.fabs(f() - results) < 1.e-12:
                raise BasicUnitTestException(str(f))
        elif args_format == "tuple":
            if not np.fabs(f(*args) - results) < 1.e-12:
                raise BasicUnitTestException(str(f))
        else:
            if not np.fabs(f(args) - results) < 1.e-12:
                raise BasicUnitTestException(str(f))

    @staticmethod
    def _check_array(f, results, args, args_format, rounded=3):
        if args_format == "none":
            if not np.array_equal(np.round(f(), rounded), np.round(results, rounded)):
                raise BasicUnitTestException(str(f))
        elif args_format == "tuple":
            if not np.array_equal(np.round(f(*args), rounded), np.round(results, rounded)):
                raise BasicUnitTestException(str(f))
        else:
            if not np.array_equal(np.round(f(args), rounded), np.round(results, rounded)):
                raise BasicUnitTestException(str(f))

    @staticmethod
    def _check_enum(f, results, args, args_format):
        if args_format == "none":
            if not int(f()) == results:
                raise BasicUnitTestException(str(f))
        elif args_format == "tuple":
            if not int(f(*args)) == results:
                raise BasicUnitTestException(str(f))
        else:
            if not int(f(args)) == results:
                raise BasicUnitTestException(str(f))

    @staticmethod
    def _check_others(f, results, args, args_format):
        if args_format == "none":
            if not f() == results:
                raise BasicUnitTestException(str(f))
        elif args_format == "tuple":
            if not f(*args) == results:
                raise BasicUnitTestException(str(f))
        else:
            if not f(args) == results:
                raise BasicUnitTestException(str(f))

    def check(self, f):
        args = self.setting_dict[f.__name__]["input"]["value"]
        results = self.setting_dict[f.__name__]["output"]["value"]
        args_format = self.setting_dict[f.__name__]["input"]["format"]
        results_format = self.setting_dict[f.__name__]["output"]["format"]
        try:
            results, results_format = self._convert_path(f, results, results_format)

            if results_format == "float":
                self._check_float(f, results, args, args_format)
            elif results_format == "array":
                self._check_array(f, results, args, args_format)
            elif results_format == "enum":
                self._check_enum(f, results, args, args_format)
            else:
                self._check_others(f, results, args, args_format)

            self._print_on_screen(f, "OK")
        except Exception as e:
            self._print_on_screen(f, "ASALI::ERROR")
            if isinstance(e, BasicUnitTestException):
                print(e.output_value(f, args, args_format, results_format))
                print(e.expected_value(results, results_format))
            else:
                print(str(traceback.format_exc()))

            sys.exit()

    def check_all(self):
        for function in self.function_to_check:
            f = getattr(self.cls, function)
            self.check(f)
