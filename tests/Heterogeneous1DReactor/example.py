from asali.reactors.het1d import Heterogeneous1DReactor


def main(cantera_input_file, gas_phase_name, surface_phase_name):
    h = Heterogeneous1DReactor(cantera_input_file, gas_phase_name, surface_phase_name)
    h.set_length([0, 0.05, 0.1, 0.15, 0.2, 0.6, 0.65], 'm')
    h.set_pressure(5, 'bar')
    h.set_catalytic_load(10, '1/m')
    h.set_volumetric_flow_rate(1., 'm3/h')
    h.set_inlet_temperature(250, 'degC')
    h.set_inlet_mass_fraction({'O2': 0.4, 'AR': 0.5, 'H2': 0.1})
    h.set_initial_coverage({'Rh(s)': 1})
    h.set_solid_density(2300, 'kg/m3')
    h.set_solid_specific_heat(750, 'J/kg/degK')
    h.set_solid_thermal_conductivity(2.5, 'W/m/degK')
    h.set_initial_solid_temperature(250, 'degC')
    h.set_energy(True)
    h.set_verbosity(False)
    h.set_initial_mass_fraction({'O2': 0.4, 'AR': 0.5, 'H2': 0.1})
    h.set_initial_temperature(250, 'degC')
    h.set_packed_bed_reactor(0.3, 'mm', 1, 'cm', 0.75)
    return h.solve([0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 2.0, 4.0, 10., 20.], 's')
