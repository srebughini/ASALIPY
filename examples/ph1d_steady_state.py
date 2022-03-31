import os

from asali.reactors.ph1d import PseudoHomogeneous1DReactor
from asali.plotters.reactor import ReactorPlotter

if __name__ == "__main__":
    p = PseudoHomogeneous1DReactor(os.path.join('files', 'H2-O2-Rh.xml'), 'gas', 'Rh_surface')
    p.set_length(2.5, 'm')
    p.set_diameter(10., 'mm')
    p.set_pressure(20, 'bar')
    p.set_catalytic_load(75, '1/m')
    p.set_volumetric_flow_rate(10, 'm3/h')
    p.set_inlet_temperature(240, 'degC')
    p.set_inlet_mass_fraction({'O2': 0.4, 'AR': 0.5, 'H2': 0.1})
    p.set_initial_coverage({'Rh(s)': 1})
    p.set_energy(True)
    p.set_inert_specie('AR')
    p.set_inert_coverage('Rh(s)')
    p.set_resolution_method("STEADYSTATE")
    p.set_gas_diffusion(True)
    p.set_verbosity(False)
    p.solve()

    plt = ReactorPlotter(p)
    plt.plot_species_mass_fraction(['H2', 'H2O', 'O2', 'AR'])
    plt.plot_species_mole_fraction(['H2', 'H2O', 'O2', 'AR'])
    plt.plot_coverage(['Rh(s)', 'H(s)', 'O(s)', 'OH(s)'])
    plt.plot_temperature()
    plt.show()
