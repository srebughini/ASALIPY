from abc import ABC


class BasicPlotter(ABC):
    def __init__(self, cls, colormap):
        """
        Abstract class that can be used to plot results
        :param cls: Object of class to be plotted
        :param colormap: String representing the color map
        """
        self.cls = cls
        self.nFig = 0
        self.colormap = colormap
