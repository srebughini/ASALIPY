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

    def plot_x_vector_and_y_vector(self, plt, xv, yv, xlabel, ylabel):
        """
        Plotting x and y
        :param plt: Matplotlib object
        :param xv: Independent variable vector
        :param yv: Dependent variable vector
        :param xlabel: Independent variable label
        :param ylabel: Dependent variable label
        :return:
        """
        self.nFig = self.nFig + 1
        plt.figure(self.nFig)
        plt.plot(xv, yv)

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)

    def plot_x_vector_and_y_matrix(self, plt, xv, ym, xlabel, ylabel, legend, title=None, colors=None):
        """
        Plotting x and array of y by adding legend
        :param plt: Matplotlib object
        :param xv: Independent variable vector
        :param ym: Dependent variable matrix
        :param xlabel: Independent variable label
        :param ylabel: Dependent variable label
        :param legend: Legend list
        :param title: Plot title
        :param colors: List of colors
        :return:
        """
        self.nFig = self.nFig + 1
        plt.figure(self.nFig)
        if colors is None:
            for yv in ym:
                plt.plot(xv, yv)
        else:
            for i, yv in enumerate(ym):
                plt.plot(xv, yv, color=colors[i])

        plt.xlabel(xlabel)
        plt.ylabel(ylabel)
        plt.legend(legend, loc='best').get_frame().set_linewidth(0.0)

        if title is not None:
            plt.title(title, loc='center')
