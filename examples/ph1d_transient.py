import os

from asali.reactors.ph1d_transient import TransientPseudoHomogeneous1DReactor
from asali.plotters.reactor import ReactorPlotter
from asali.savers.reactor import ReactorSaver

if __name__ == "__main__":
    p = TransientPseudoHomogeneous1DReactor(os.path.join('files', 'H2-O2-Rh.yaml'), 'gas', 'Rh_surface')
    p.set_length(2.5, 'm')
    p.set_diameter(10., 'mm')
    p.set_pressure(20, 'bar')
    p.set_catalytic_load(75, '1/m')
    p.set_volumetric_flow_rate(10, 'm3/h')
    p.set_inlet_temperature(240, 'degC')
    p.set_inlet_mass_fraction({'O2': 0.4, 'AR': 0.5, 'H2': 0.1})
    p.set_initial_coverage({'Rh(s)': 1})
    p.set_energy(True)
    p.set_initial_mass_fraction({'O2': 0.4, 'AR': 0.5, 'H2': 0.1})
    p.set_inert_specie('AR')
    p.set_inert_coverage('Rh(s)')
    p.set_initial_temperature(240, 'degC')
    p.set_gas_diffusion(True)
    p.set_verbosity(False)
    p.set_relative_tolerance(1.e-04)
    p.set_absolute_tolerance(1.e-04)
    p.solve([0, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06], 's')

    svr = ReactorSaver(p)
    svr.save_using_mole_fraction(os.path.join('files', 'output_ph1d_transient_mole_fraction.xlsx'),
                                 species_names=['H2', 'H2O', 'O2', 'AR'],
                                 coverage_names=['Rh(s)', 'H(s)', 'O(s)', 'OH(s)'])

    plt = ReactorPlotter(p, colormap="Greens")
    plt.plot_species_mass_fraction(['H2', 'H2O', 'O2'])
    plt.plot_species_mole_fraction(['H2', 'H2O', 'O2'])
    plt.plot_species_mole_fraction(['AR'])
    plt.plot_coverage(['Rh(s)', 'H(s)', 'O(s)', 'OH(s)'])
    plt.plot_temperature()
    plt.show()