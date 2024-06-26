import os

from asali.reactors.batch import BatchReactor
from asali.plotters.reactor import ReactorPlotter

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

    plt = ReactorPlotter(b)
    plt.plot_species_mass_fraction(['H2', 'H2O', 'O2'])
    plt.plot_species_mole_fraction(['H2', 'H2O', 'O2'])
    plt.plot_temperature()
    plt.show()