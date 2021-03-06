from asali.reactors.basic import ResolutionMethod
from tests.basic_unit_test import BasicUnitTest

import numpy as np

import os


class ReactorUnitTest(BasicUnitTest):
    def __init__(self, cls, file_name, folder_path="tests/"):
        super().__init__(cls=cls, file_name=file_name, folder_path=folder_path)

    def _convert_path(self, f, results, results_format):
        if results_format == "path":
            if self.cls.resolution_method == ResolutionMethod.STEADYSTATE:
                file_name = "{}.asali".format(f.__name__)
                results_new = np.loadtxt(os.path.join(self.folder_path, self.cls.__class__.__name__, file_name))
                results_format_new = "array"
                return results_new, results_format_new

            if self.cls.resolution_method == ResolutionMethod.TRANSIENT:
                file_name = "{}.asali".format(f.__name__)
                results_new = np.loadtxt(os.path.join(self.folder_path, self.cls.__class__.__name__, file_name))
                results_format_new = "array"
                return results_new, results_format_new

            file_name = "{}.asali".format(f.__name__)
            results_new = np.loadtxt(os.path.join(self.folder_path, self.cls.__class__.__name__, file_name))
            results_format_new = "array"
            return results_new, results_format_new

        return results, results_format
