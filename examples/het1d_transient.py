import os

from asali.plotters.reactor import ReactorPlotter
from asali.reactors.het1d_transient import TransientHeterogeneous1DReactor
from asali.savers.reactor import ReactorSaver

if __name__ == "__main__":
    h = TransientHeterogeneous1DReactor(os.path.join('files', 'H2-O2-Rh.yaml'), 'gas', 'Rh_surface')
    h.set_length([0, 0.5, 1.0, 1.5, 2.0], 'm')
    h.set_pressure(5, 'bar')
    h.set_catalytic_load(35, '1/m')
    h.set_volumetric_flow_rate(15., 'm3/h')
    h.set_inlet_temperature(250, 'degC')
    h.set_inlet_mass_fraction({'O2': 0.4, 'AR': 0.5, 'H2': 0.1})
    h.set_initial_coverage({'Rh(s)': 1})
    h.set_solid_density(2300, 'kg/m3')
    h.set_solid_specific_heat(750, 'J/kg/degK')
    h.set_solid_thermal_conductivity(2.5, 'W/m/degK')
    h.set_initial_solid_temperature(250, 'degC')
    h.set_energy(True)
    h.set_gas_diffusion(True)
    h.set_verbosity(False)
    h.set_initial_mass_fraction({'O2': 0.4, 'AR': 0.5, 'H2': 0.1})
    h.set_initial_temperature(250, 'degC')
    h.set_packed_bed_reactor(0.3, 'mm', 1, 'cm', 0.75)
    h.solve([0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 2.0, 4.0, 10., 20.], 's')

    svr = ReactorSaver(h)
    svr.save_using_mole_fraction(os.path.join('files', 'output_het1d_steady_state_mole_fraction.xlsx'),
                                 species_names=['H2', 'H2O', 'O2'],
                                 coverage_names=['H(s)', 'O(s)', 'Rh(s)'])

    plt = ReactorPlotter(h)
    plt.plot_species_mass_fraction(['H2', 'H2O', 'O2'])
    plt.plot_species_mole_fraction(['H2', 'H2O', 'O2'])
    plt.plot_coverage(['Rh(s)', 'H(s)', 'O(s)', 'OH(s)'])
    plt.plot_temperature()
    plt.show()
