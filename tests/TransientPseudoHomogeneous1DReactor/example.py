from asali.reactors.ph1d_transient import TransientPseudoHomogeneous1DReactor


def main(cantera_input_file, gas_phase_name, surface_phase_name):
    p = TransientPseudoHomogeneous1DReactor(cantera_input_file, gas_phase_name, surface_phase_name)
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
    return p.solve([0, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06], 's')
