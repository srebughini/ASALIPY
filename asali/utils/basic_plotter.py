from abc import ABC

import json


class BasicPlotter(ABC):
    def __init__(self, cls, config_file=None):
        self.cls = cls
        self.nFig = 0

        if config_file is not None:
            with open(config_file) as json_file:
                self.config_dict = json.load(json_file)
        else:
            self.config_dict[self.cls.__class__.__name__]["colormap"] = None
