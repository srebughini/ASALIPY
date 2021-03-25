from asali.reactors.basic import ReactorType, ResolutionMethod
from asali.plotters.basic import BasicPlotter

import matplotlib.pyplot as plt
import matplotlib
import numpy as np


class ReactorPlotter(BasicPlotter):
    def __init__(self, reactor_cls, colormap="Blues"):
        super().__init__(cls=reactor_cls, colormap=colormap)

        if self.cls.reactor_type not in [r for r in ReactorType]:
            raise Exception("ASALI::ERROR::Unknown class type: ", str(type(self.cls)))

        if not self.cls.is_solved:
            raise Exception("ASALI::ERROR::Nothing to plot, no solution found")

        cmap = matplotlib.cm.get_cmap(self.colormap)
        self.colors = cmap(np.linspace(0.2, 1., num=self.cls.tspan.size))

    def _plot_species_mass_fraction_for_batch_and_cstr(self, species_names):
        self.nFig = self.nFig + 1
        plt.figure(self.nFig)
        for s in species_names:
            plt.plot(self.cls.tspan, self.cls.y_sol[:, self.cls.gas.species_index(s)])

        plt.ylabel('Mass fraction')
        plt.xlabel('Time [s]')
        plt.legend(species_names, loc='best').get_frame().set_linewidth(0.0)

    def _plot_species_mole_fraction_for_batch_and_cstr(self, species_names):
        self.nFig = self.nFig + 1
        plt.figure(self.nFig)

        for s in species_names:
            plt.plot(self.cls.tspan, self.cls.x_sol[:, self.cls.gas.species_index(s)])

        plt.ylabel('Mole fraction')
        plt.xlabel('Time [s]')
        plt.legend(species_names, loc='best').get_frame().set_linewidth(0.0)

    def _plot_coverage_for_batch_and_cstr(self, coverage_names):
        self.nFig = self.nFig + 1
        plt.figure(self.nFig)

        for s in coverage_names:
            plt.plot(self.cls.tspan, self.cls.coverage_sol[:, self.cls.surf.species_index(s)])

        plt.ylabel('Site fraction')
        plt.xlabel('Time [s]')
        plt.legend(coverage_names, loc='best').get_frame().set_linewidth(0.0)

    def _plot_temperature_for_batch_and_cstr(self):
        self.nFig = self.nFig + 1
        plt.figure(self.nFig)
        plt.plot(self.cls.tspan, self.cls.temperature_sol)

        plt.ylabel('Temperature [K]')
        plt.xlabel('Time [s]')

    def _plot_species_mass_fraction_pseudo_homogeneous_steady_state(self, species_names):
        self.nFig = self.nFig + 1
        plt.figure(self.nFig)
        for s in species_names:
            plt.plot(self.cls.length, self.cls.y_sol[:, self.cls.gas.species_index(s)])

        plt.ylabel('Mass fraction')
        plt.xlabel('Length [m]')
        plt.legend(species_names, loc='best').get_frame().set_linewidth(0.0)

    def _plot_species_mole_fraction_pseudo_homogeneous_steady_state(self, species_names):
        self.nFig = self.nFig + 1
        plt.figure(self.nFig)

        for s in species_names:
            plt.plot(self.cls.length, self.cls.x_sol[:, self.cls.gas.species_index(s)])

        plt.ylabel('Mole fraction')
        plt.xlabel('Length [m]')
        plt.legend(species_names, loc='best').get_frame().set_linewidth(0.0)

    def _plot_coverage_pseudo_homogeneous_steady_state(self, coverage_names):
        self.nFig = self.nFig + 1
        plt.figure(self.nFig)
        for s in coverage_names:
            plt.plot(self.cls.length, self.cls.coverage_sol[:, self.cls.surf.species_index(s)])

        plt.ylabel('Site fraction')
        plt.xlabel('Length [m]')
        plt.legend(coverage_names, loc='best').get_frame().set_linewidth(0.0)

    def _plot_temperature_pseudo_homogeneous_steady_state(self):
        self.nFig = self.nFig + 1
        plt.figure(self.nFig)
        plt.plot(self.cls.length, self.cls.temperature_sol)

        plt.ylabel('Temperature [K]')
        plt.xlabel('Length [m]')

    def _plot_species_mass_fraction_pseudo_homogeneous_transient(self, species_names):
        for s in species_names:
            self.nFig = self.nFig + 1
            plt.figure(self.nFig)
            legend = list()
            for k, t in enumerate(self.cls.tspan):
                plt.plot(self.cls.length, self.cls.y_sol[k][:, self.cls.gas.species_index(s)], color=self.colors[k])
                legend.append("Time: " + str(round(t, 3)) + " s")

            plt.ylabel('Mass fraction')
            plt.xlabel('Length [m]')
            plt.legend(legend, loc='best').get_frame().set_linewidth(0.0)
            plt.title(s, loc='center')

    def _plot_species_mole_fraction_pseudo_homogeneous_transient(self, species_names):
        for s in species_names:
            self.nFig = self.nFig + 1
            plt.figure(self.nFig)
            legend = list()
            for k, t in enumerate(self.cls.tspan):
                plt.plot(self.cls.length, self.cls.x_sol[k][:, self.cls.gas.species_index(s)], color=self.colors[k])
                legend.append("Time: " + str(round(t, 3)) + " s")

            plt.ylabel('Mole fraction')
            plt.xlabel('Length [m]')
            plt.legend(legend, loc='best').get_frame().set_linewidth(0.0)
            plt.title(s, loc='center')

    def _plot_coverage_pseudo_homogeneous_transient(self, coverage_names):
        for s in coverage_names:
            self.nFig = self.nFig + 1
            plt.figure(self.nFig)
            legend = list()
            for k, t in enumerate(self.cls.tspan):
                plt.plot(self.cls.length, self.cls.coverage_sol[k][:, self.cls.surf.species_index(s)], color=self.colors[k])
                legend.append("Time: " + str(round(t, 3)) + " s")

            plt.ylabel('Site fraction')
            plt.xlabel('Length [m]')
            plt.legend(legend, loc='best').get_frame().set_linewidth(0.0)
            plt.title(s, loc='center')

    def _plot_temperature_pseudo_homogeneous_transient(self):
        self.nFig = self.nFig + 1
        plt.figure(self.nFig)
        legend = list()
        for k, t in enumerate(self.cls.tspan):
            plt.plot(self.cls.length, self.cls.temperature_sol[k], color=self.colors[k])
            legend.append("Time: " + str(round(t, 3)) + " s")

        plt.ylabel('Temperature [K]')
        plt.xlabel('Length [m]')
        plt.legend(legend, loc='best').get_frame().set_linewidth(0.0)

    def _plot_species_mass_fraction_heterogeneous(self, species_names):
        for s in species_names:
            self.nFig = self.nFig + 1
            plt.figure(self.nFig)
            legend = list()
            for k, t in enumerate(self.cls.tspan):
                plt.plot(self.cls.length, self.cls.y_sol[k][:, self.cls.gas.species_index(s)], color=self.colors[k])
                legend.append("Time: " + str(round(t, 3)) + " s")

            plt.ylabel('Gas mass fraction')
            plt.xlabel('Length [m]')
            plt.legend(legend, loc='best').get_frame().set_linewidth(0.0)
            plt.title(s, loc='center')

            self.nFig = self.nFig + 1
            plt.figure(self.nFig)
            legend = list()
            for k, t in enumerate(self.cls.tspan):
                plt.plot(self.cls.length, self.cls.y_sol_wall[k][:, self.cls.gas.species_index(s)], color=self.colors[k])
                legend.append("Time: " + str(round(t, 3)) + " s")

            plt.ylabel('Solid mass fraction')
            plt.xlabel('Length [m]')
            plt.legend(legend, loc='best').get_frame().set_linewidth(0.0)
            plt.title(s, loc='center')

    def _plot_species_mole_fraction_heterogeneous(self, species_names):
        for s in species_names:
            self.nFig = self.nFig + 1
            plt.figure(self.nFig)
            legend = list()
            for k, t in enumerate(self.cls.tspan):
                plt.plot(self.cls.length, self.cls.x_sol[k][:, self.cls.gas.species_index(s)], color=self.colors[k])
                legend.append("Time: " + str(round(t, 3)) + " s")

            plt.ylabel('Gas mole fraction')
            plt.xlabel('Length [m]')
            plt.legend(legend, loc='best').get_frame().set_linewidth(0.0)
            plt.title(s, loc='center')

            self.nFig = self.nFig + 1
            plt.figure(self.nFig)
            legend = list()
            for k, t in enumerate(self.cls.tspan):
                plt.plot(self.cls.length, self.cls.x_sol_wall[k][:, self.cls.gas.species_index(s)], color=self.colors[k])
                legend.append("Time: " + str(round(t, 3)) + " s")

            plt.ylabel('Solid mole fraction')
            plt.xlabel('Length [m]')
            plt.legend(legend, loc='best').get_frame().set_linewidth(0.0)
            plt.title(s, loc='center')

    def _plot_temperature_heterogeneous(self):
        self.nFig = self.nFig + 1
        plt.figure(self.nFig)
        legend = list()
        for k, t in enumerate(self.cls.tspan):
            plt.plot(self.cls.length, self.cls.temperature_sol[k], color=self.colors[k])
            legend.append("Time: " + str(round(t, 3)) + " s")

        plt.ylabel('Gas temperature [K]')
        plt.xlabel('Length [m]')
        plt.legend(legend, loc='best').get_frame().set_linewidth(0.0)

        self.nFig = self.nFig + 1
        plt.figure(self.nFig)
        legend = list()
        for k, t in enumerate(self.cls.tspan):
            plt.plot(self.cls.length, self.cls.temperature_sol_wall[k], color=self.colors[k])
            legend.append("Time: " + str(round(t, 3)) + " s")

        plt.ylabel('Solid temperature [K]')
        plt.xlabel('Length [m]')
        plt.legend(legend, loc='best').get_frame().set_linewidth(0.0)

    def plot_species_mass_fraction(self, species_names):
        if self.cls.reactor_type == ReactorType.BATCH or self.cls.reactor_type == ReactorType.CSTR:
            self._plot_species_mass_fraction_for_batch_and_cstr(species_names)
        elif self.cls.reactor_type == ReactorType.PSEUDOHOMOGENEOUSPFR:
            if self.cls.resolution_method == ResolutionMethod.STEADYSTATE:
                self._plot_species_mass_fraction_pseudo_homogeneous_steady_state(species_names)
            elif self.cls.resolution_method == ResolutionMethod.TRANSIENT:
                self._plot_species_mass_fraction_pseudo_homogeneous_transient(species_names)
        elif self.cls.reactor_type == ReactorType.HETEROGENEOUSPRF:
            self._plot_species_mass_fraction_heterogeneous(species_names)

    def plot_species_mole_fraction(self, species_names):
        if self.cls.reactor_type == ReactorType.BATCH or self.cls.reactor_type == ReactorType.CSTR:
            self._plot_species_mole_fraction_for_batch_and_cstr(species_names)
        elif self.cls.reactor_type == ReactorType.PSEUDOHOMOGENEOUSPFR:
            if self.cls.resolution_method == ResolutionMethod.STEADYSTATE:
                self._plot_species_mole_fraction_pseudo_homogeneous_steady_state(species_names)
            elif self.cls.resolution_method == ResolutionMethod.TRANSIENT:
                self._plot_species_mole_fraction_pseudo_homogeneous_transient(species_names)
        elif self.cls.reactor_type == ReactorType.HETEROGENEOUSPRF:
            self._plot_species_mole_fraction_heterogeneous(species_names)

    def plot_coverage(self, coverage_names):
        if self.cls.reactor_type == ReactorType.BATCH or self.cls.reactor_type == ReactorType.CSTR:
            self._plot_coverage_for_batch_and_cstr(coverage_names)
        elif self.cls.reactor_type == ReactorType.PSEUDOHOMOGENEOUSPFR:
            if self.cls.resolution_method == ResolutionMethod.STEADYSTATE:
                self._plot_coverage_pseudo_homogeneous_steady_state(coverage_names)
            elif self.cls.resolution_method == ResolutionMethod.TRANSIENT:
                self._plot_coverage_pseudo_homogeneous_transient(coverage_names)
        elif self.cls.reactor_type == ReactorType.HETEROGENEOUSPRF:
            self._plot_coverage_pseudo_homogeneous_transient(coverage_names)

    def plot_temperature(self):
        if self.cls.reactor_type == ReactorType.BATCH or self.cls.reactor_type == ReactorType.CSTR:
            self._plot_temperature_for_batch_and_cstr()
        elif self.cls.reactor_type == ReactorType.PSEUDOHOMOGENEOUSPFR:
            if self.cls.resolution_method == ResolutionMethod.STEADYSTATE:
                self._plot_temperature_pseudo_homogeneous_steady_state()
            elif self.cls.resolution_method == ResolutionMethod.TRANSIENT:
                self._plot_temperature_pseudo_homogeneous_transient()
        elif self.cls.reactor_type == ReactorType.HETEROGENEOUSPRF:
            self._plot_temperature_heterogeneous()

    @staticmethod
    def show():
        plt.show()
