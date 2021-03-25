from asali.reactors.ph1d import PseudoHomogeneous1DReactor
from asali.plotters.reactor import ReactorPlotter

if __name__ == "__main__":
    p = PseudoHomogeneous1DReactor('H2-O2-Rh.xml', 'gas', 'Rh_surface')
    p.set_length(0.5, 'cm')
    p.set_diameter(10., 'mm')
    p.set_pressure(5, 'bar')
    p.set_catalytic_load(150, '1/m')
    p.set_volumetric_flow_rate(10, 'm3/h')
    p.set_inlet_temperature(240, 'degC')
    p.set_inlet_mass_fraction({'O2': 0.4, 'AR': 0.5, 'H2': 0.1})
    p.set_initial_coverage({'Rh(s)': 1})
    p.set_energy(1)
    p.set_initial_mass_fraction({'AR': 1})
    p.set_initial_temperature(25, 'degC')
    p.set_resolution_method("TRANSIENT")
    p.solve([0, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1], 's')

    plt = ReactorPlotter(p, colormap="Greens")
    plt.plot_species_mass_fraction(['H2', 'H2O', 'O2'])
    plt.plot_species_mole_fraction(['H2', 'H2O', 'O2'])
    plt.plot_species_mole_fraction(['AR'])
    plt.plot_coverage(['Rh(s)', 'H(s)', 'O(s)', 'OH(s)'])
    plt.plot_temperature()
    plt.show()