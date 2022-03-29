from asali.reactors.batch import BatchReactor


def main(cantera_input_file, gas_phase_name, surface_phase_name):
    b = BatchReactor(cantera_input_file, gas_phase_name, surface_phase_name)
    b.set_volume(10., 'mm3')
    b.set_pressure(5, 'bar')
    b.set_catalytic_load(15, '1/m')
    b.set_initial_mole_fraction({'O2': 0.4, 'AR': 0.5, 'H2': 0.1})
    b.set_initial_temperature(120, 'degC')
    b.set_initial_coverage({'Rh(s)': 1})
    b.set_energy(1)
    return b.solve([0, 0.1, 0.5, 5], 's')
