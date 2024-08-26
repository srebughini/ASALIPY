import os

from asali.reactors.cstr import CstrReactor
from asali.plotters.reactor import ReactorPlotter
from asali.savers.reactor import ReactorSaver

if __name__ == "__main__":
    c = CstrReactor(os.path.join('files', 'H2-O2-Rh.yaml'), 'gas', 'Rh_surface')
    c.set_volume(5., 'dm3')
    c.set_pressure(5, 'bar')
    c.set_catalytic_load(150, '1/m')
    c.set_volumetric_flow_rate(1, 'm3/h')
    c.set_inlet_temperature(120, 'degC')
    c.set_inlet_mass_fraction({'O2': 0.4, 'AR': 0.5, 'H2': 0.1})
    c.set_initial_temperature(50, 'degC')
    c.set_initial_mass_fraction({'AR': 1.0})
    c.set_initial_coverage({'Rh(s)': 1})
    c.set_energy(1)
    c.solve(list(range(0, 30, 1)), 's')

    svr = ReactorSaver(c)
    svr.save_using_mole_fraction(os.path.join('files', 'output_cstr_mole_fraction.xlsx'),
                                 species_names=['H2', 'H2O', 'O2', 'AR'],
                                 coverage_names=['Rh(s)', 'H(s)', 'O(s)', 'OH(s)'])

    plt = ReactorPlotter(c)
    plt.plot_species_mass_fraction(['H2', 'H2O', 'O2'])
    plt.plot_species_mole_fraction(['H2', 'H2O', 'O2'])
    plt.plot_species_mole_fraction(['AR'])
    plt.plot_coverage(['Rh(s)', 'H(s)', 'O(s)', 'OH(s)'])
    plt.plot_temperature()
    plt.show()
