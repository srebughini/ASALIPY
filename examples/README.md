## **Batch Reactor**

This [example](https://github.com/srebughini/ASALIPY/blob/main/examples/batch.py) show how to solve a **batch reactor** for the [catalytic combustion of hydrogen over rhodium](https://www.detchem.com/mechanisms).

```python
from asali.reactors.batch import BatchReactor

if __name__ == "__main__":
    b = BatchReactor('H2-O2-Rh.xml', 'gas', 'Rh_surface')  # Initialize reactor class
    b.set_volume(10., 'mm3')  # Set reactor volume in [mm3]
    b.set_pressure(5, 'bar')  # Set reactor pressure in [bar]
    b.set_catalytic_load(15, '1/m')  # Set catalytic load in [1/m]
    b.set_initial_mole_fraction({'O2': 0.4, 'AR': 0.5, 'H2': 0.1})  # Set reactor initial composition using mole fraction
    b.set_initial_temperature(120, 'degC')  # Set reactor initial temperature in [°C]
    b.set_initial_coverage({'Rh(s)': 1})  # Set reactor initial coverage
    b.set_energy(1)  # Enable energy balance
    b.solve([0, 0.1, 0.5, 5], 's')  # Solve for different time steps
```

## **Continuous Stirred Tank Reactor**

This [example](https://github.com/srebughini/ASALIPY/blob/main/examples/cstr.py) show how to solve a **continuous stirred tank reactor** for the [catalytic combustion of hydrogen over rhodium](https://www.detchem.com/mechanisms).

```python
from asali.reactors.cstr import CstrReactor

if __name__ == "__main__":
    c = CstrReactor('H2-O2-Rh.xml', 'gas', 'Rh_surface')  # Initialize reactor class
    c.set_volume(5., 'm3')  # Set reactor volume in [m3]
    c.set_pressure(5, 'bar')  # Set reactor pressure in [bar]
    c.set_catalytic_load(150, '1/m')  # Set catalytic load in [1/m]
    c.set_volumetric_flow_rate(1, 'm3/h')  # Set volumetric flow rate in [m3/h]
    c.set_inlet_temperature(120, 'degC')  # Set inlet gas temperature in [°C]
    c.set_inlet_mass_fraction({'O2': 0.4, 'AR': 0.5, 'H2': 0.1})  # Set inlet gas composition using mass fraction
    c.set_initial_mass_fraction({'AR': 1})  # Set reactor initial composition using mass fraction
    c.set_initial_temperature(25, 'degC')  # Set reactor initial temperature in [°C]
    c.set_initial_coverage({'Rh(s)': 1})  # Set reactor initial coverage
    c.set_energy(1)  # Enable energy balance
    c.solve([0, 0.1, 0.5, 1, 2, 5], 'min')  # Solve for different time steps in [min]
```

## **1-D Pseudo-Homogeneous Plug Flow Reactor**

This [example](https://github.com/srebughini/ASALIPY/blob/main/examples/ph1d.py) show how to solve a **1-D pseudo-homogeneous plug flow reactor** for the [catalytic combustion of hydrogen over rhodium](https://www.detchem.com/mechanisms).

```python
from asali.reactors.ph1d import PseudoHomogeneous1DReactor

if __name__ == "__main__":
    p = PseudoHomogeneous1DReactor('H2-O2-Rh.xml', 'gas', 'Rh_surface')  # Initialize reactor class
    p.set_length(0.5, 'cm')  # Set reactor length in [cm]
    p.set_diameter(10., 'mm')  # Set reactor diameter in [mm]
    p.set_pressure(5, 'bar')  # Set reactor pressure in [bar]
    p.set_catalytic_load(150, '1/m')  # Set catalytic load in [1/m]
    p.set_volumetric_flow_rate(10, 'm3/h')  # Set volumetric flow rate in [m3/h]
    p.set_inlet_temperature(240, 'degC')  # Set inlet gas temperature in [°C]
    p.set_inlet_mass_fraction({'O2': 0.4, 'AR': 0.5, 'H2': 0.1})  # Set inlet gas composition using mass fraction
    p.set_initial_coverage({'Rh(s)': 1})  # Set reactor initial coverage
    p.set_energy(1)  # Enable energy balance
    p.set_initial_mass_fraction({'AR': 1})  # Set reactor initial composition using mass fraction
    p.set_initial_temperature(25, 'degC')  # Set reactor initial temperature in [°C]
    p.set_resolution_method("TRANSIENT")  # Set resolution method
    p.solve([0, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1], 's')  # Solve for different time steps in [s]
```

## **1-D Heterogeneous Plug Flow Reactor**

This [example](https://github.com/srebughini/ASALIPY/blob/main/examples/het1d.py) show how to solve a **1-D heterogeneous plug flow reactor** for the [catalytic combustion of hydrogen over rhodium](https://www.detchem.com/mechanisms).

```python
from asali.reactors.het1d import Heterogeneous1DReactor

if __name__ == "__main__":
    h = Heterogeneous1DReactor('H2-O2-Rh.xml', 'gas', 'Rh_surface')  # Initialize reactor class
    h.set_length([0, 0.05, 0.1, 0.15, 0.2, 0.6, 0.65], 'm')  # Set reactor length in [m]
    h.set_pressure(5, 'bar')  # Set reactor pressure in [bar]
    h.set_catalytic_load(10, '1/m')  # Set catalytic load in [1/m]
    h.set_volumetric_flow_rate(1., 'm3/h')  # Set volumetric flow rate in [m3/h]
    h.set_inlet_temperature(250, 'degC')  # Set inlet gas temperature in [°C]
    h.set_inlet_mass_fraction({'O2': 0.4, 'AR': 0.5, 'H2': 0.1})  # Set inlet gas composition using mass fraction
    h.set_initial_coverage({'Rh(s)': 1})  # Set reactor initial coverage
    h.set_solid_density(2300, 'kg/m3')  # Set catalyst density in [kg/m3]
    h.set_solid_specific_heat(750, 'J/kg/degK')  # Set catalyst specific heat in [J/kg/K]
    h.set_solid_thermal_conductivity(2.5, 'W/m/degK')  # Set catalyst thermal conductivity in [W/m/K]
    h.set_initial_solid_temperature(250, 'degC')  # Set initial catalyst temperature in [°C]
    h.set_energy(True)  # Enable energy balance
    h.set_initial_mass_fraction({'O2': 0.4, 'AR': 0.5, 'H2': 0.1})  # Set reactor initial composition using mass fraction
    h.set_initial_temperature(250, 'degC')  # Set reactor initial temperature in [°C]
    h.set_packed_bed_reactor(0.3, 'mm', 1, 'cm', 0.75)  # Set packed bed reactor properties
    h.solve([0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 2.0, 4.0, 10., 20.], 's')  # Solve for different time steps in [s]
```

## **Reactor Plotter**

This [example](https://github.com/srebughini/ASALIPY/blob/main/examples/het1d.py) show how to **solve and plot** 1-D heterogeneous plug flow reactor for the [catalytic combustion of hydrogen over rhodium](https://www.detchem.com/mechanisms).

```python
from asali.reactors.het1d import Heterogeneous1DReactor
from asali.plotters.reactor import ReactorPlotter

if __name__ == "__main__":
    h = Heterogeneous1DReactor('H2-O2-Rh.xml', 'gas', 'Rh_surface')
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
    h.set_initial_mass_fraction({'O2': 0.4, 'AR': 0.5, 'H2': 0.1})
    h.set_initial_temperature(250, 'degC')
    h.set_packed_bed_reactor(0.3, 'mm', 1, 'cm', 0.75)
    h.solve([0, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 2.0, 4.0, 10., 20.], 's')

    plt = ReactorPlotter(h)  # Initialize plotting object
    plt.plot_species_mass_fraction(['H2', 'H2O', 'O2'])  # Plot mass fraction from species names
    plt.plot_species_mole_fraction(['H2', 'H2O', 'O2'])  # Plot mole fraction from species names
    plt.plot_coverage(['Rh(s)', 'H(s)', 'O(s)', 'OH(s)'])  # Plot coverage from coverage names
    plt.plot_temperature()  # Plot temperature
    plt.show()  # Show plots
```
