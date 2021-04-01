from abc import ABC


class BasicPlotter(ABC):
    def __init__(self, cls, colormap):
        self.cls = cls
        self.nFig = 0
        self.colormap = colormap
