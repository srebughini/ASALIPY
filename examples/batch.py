import os

from asali.reactors.batch import BatchReactor
from asali.plotters.reactor import ReactorPlotter
from asali.savers.reactor import ReactorSaver

if __name__ == "__main__":
    b = BatchReactor(os.path.join('files', 'H2-O2-Rh.yaml'), 'gas', 'Rh_surface')
    b.set_user_defined_kinetic_model(os.path.join('files', 'H2-O2.json'))
    b.set_volume(10., 'm3')
    b.set_pressure(5, 'bar')
    b.set_catalytic_load(15, '1/m')
    b.set_initial_mole_fraction({'O2': 0.4, 'AR': 0.5, 'H2': 0.1})
    b.set_initial_temperature(120, 'degC')
    b.set_energy(1)
    b.solve([0, 0.1, 0.5, 5], 's')

    c = BatchReactor(os.path.join('files', 'H2-O2-Rh.yaml'), 'gas', 'Rh_surface')
    c.set_results(b.get_time("s"), b.get_results())

    svr = ReactorSaver(b)
    svr.save_using_mass_fraction(os.path.join('files', 'output_batch_mass_fraction.xlsx'),
                                 species_names=['H2', 'H2O', 'O2'])

    plt = ReactorPlotter(c, style='classic')
    plt.set_rc_params({'toolbar': 'None'})
    plt.plot_species_mass_fraction(['H2', 'H2O', 'O2'])
    plt.plot_species_mole_fraction(['H2', 'H2O', 'O2'])
    plt.plot_temperature()
    plt.show()
