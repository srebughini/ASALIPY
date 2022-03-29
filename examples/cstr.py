import os

from asali.reactors.cstr import CstrReactor
from asali.plotters.reactor import ReactorPlotter


def main(kinetic_schem_file_path, plot=True):
    c = CstrReactor(kinetic_schem_file_path, 'gas', 'Rh_surface')
    c.set_volume(5., 'm3')
    c.set_pressure(5, 'bar')
    c.set_catalytic_load(150, '1/m')
    c.set_volumetric_flow_rate(1, 'm3/h')
    c.set_inlet_temperature(120, 'degC')
    c.set_inlet_mass_fraction({'O2': 0.4, 'AR': 0.5, 'H2': 0.1})
    c.set_initial_mass_fraction({'AR': 1})
    c.set_initial_temperature(25, 'degC')
    c.set_initial_coverage({'Rh(s)': 1})
    c.set_energy(1)
    solution = c.solve([0, 0.1, 0.5, 1, 2, 5], 'min')

    if plot:
        plt = ReactorPlotter(c)
        plt.plot_species_mass_fraction(['H2', 'H2O', 'O2'])
        plt.plot_species_mole_fraction(['H2', 'H2O', 'O2'])
        plt.plot_species_mole_fraction(['AR'])
        plt.plot_coverage(['Rh(s)', 'H(s)', 'O(s)', 'OH(s)'])
        plt.plot_temperature()
        plt.show()

    return solution


if __name__ == "__main__":
    main(os.path.join('files', 'H2-O2-Rh.xml'))
