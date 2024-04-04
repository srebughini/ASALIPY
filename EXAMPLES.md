## **Batch Reactor**

This [example](https://github.com/srebughini/ASALIPY/blob/main/examples/batch.py) show how to solve a **batch reactor**
for the [catalytic combustion of hydrogen over rhodium](https://www.detchem.com/mechanisms).

```python
import os
from asali.reactors.batch import BatchReactor

if __name__ == "__main__":
    b = BatchReactor(os.path.join('examples/files', 'H2-O2-Rh.xml'), 'gas', 'Rh_surface')  # Initialize reactor class
    b.set_volume(10., 'mm3')  # Set reactor volume in [mm3]
    b.set_pressure(5, 'bar')  # Set reactor pressure in [bar]
    b.set_catalytic_load(15, '1/m')  # Set catalytic load in [1/m]
    b.set_initial_mole_fraction(
        {'O2': 0.4, 'AR': 0.5, 'H2': 0.1})  # Set reactor initial composition using mole fraction
    b.set_initial_temperature(120, 'degC')  # Set reactor initial temperature in [°C]
    b.set_initial_coverage({'Rh(s)': 1})  # Set reactor initial coverage
    b.set_energy(1)  # Enable energy balance
    b.solve([0, 0.1, 0.5, 5], 's')  # Solve for different time steps
```

## **Continuous Stirred Tank Reactor**

This [example](https://github.com/srebughini/ASALIPY/blob/main/examples/cstr.py) show how to solve a **continuous
stirred tank reactor** for the [catalytic combustion of hydrogen over rhodium](https://www.detchem.com/mechanisms).

```python
import os
from asali.reactors.cstr import CstrReactor

if __name__ == "__main__":
    c = CstrReactor(os.path.join('files', 'H2-O2-Rh.xml'), 'gas', 'Rh_surface')  # Initialize reactor class
    c.set_volume(5., 'dm3')  # Set reactor volume in [dm3]
    c.set_pressure(5, 'bar')  # Set reactor pressure in [bar]
    c.set_catalytic_load(150, '1/m')  # Set catalytic load in [1/m]
    c.set_volumetric_flow_rate(1, 'm3/h')  # Set volumetric flow rate in [m3/h]
    c.set_inlet_temperature(120, 'degC')  # Set inlet gas temperature in [°C]
    c.set_inlet_mass_fraction({'O2': 0.4, 'AR': 0.5, 'H2': 0.1})  # Set inlet gas composition using mass fraction
    c.set_initial_mass_fraction({'AR': 1})  # Set reactor initial composition using mass fraction
    c.set_initial_temperature(50, 'degC')  # Set reactor initial temperature in [°C]
    c.set_initial_coverage({'Rh(s)': 1})  # Set reactor initial coverage
    c.set_energy(1)  # Enable energy balance
    c.solve(list(range(0, 30, 1)), 's')  # Solve for different time steps in [s]
```

## **1-D Pseudo-Homogeneous Plug Flow Reactor**
### **Transient**

This [example](https://github.com/srebughini/ASALIPY/blob/main/examples/ph1d_transient.py) show how to solve a **transient 1-D
pseudo-homogeneous plug flow reactor** for
the [catalytic combustion of hydrogen over rhodium](https://www.detchem.com/mechanisms).

```python
import os
from asali.reactors.ph1d import PseudoHomogeneous1DReactor

if __name__ == "__main__":
    p = PseudoHomogeneous1DReactor(os.path.join('files', 'H2-O2-Rh.xml'), 'gas', 'Rh_surface')  # Initialize reactor class
    p.set_length(2.5, 'm')  # Set reactor length in [m]
    p.set_diameter(10., 'mm')  # Set reactor diameter in [mm]
    p.set_pressure(20, 'bar')  # Set reactor pressure in [bar]
    p.set_catalytic_load(75, '1/m')  # Set catalytic load in [1/m]
    p.set_volumetric_flow_rate(10, 'm3/h')  # Set volumetric flow rate in [m3/h]
    p.set_inlet_temperature(240, 'degC')  # Set inlet gas temperature in [°C]
    p.set_inlet_mass_fraction({'O2': 0.4, 'AR': 0.5, 'H2': 0.1})  # Set inlet gas composition using mass fraction
    p.set_initial_coverage({'Rh(s)': 1})  # Set reactor initial coverage
    p.set_energy(True)  # Enable energy balance
    p.set_initial_mass_fraction( {'O2': 0.4, 'AR': 0.5, 'H2': 0.1})  # Set reactor initial composition using mass fraction
    p.set_inert_specie('AR')  # Set inert specie
    p.set_inert_coverage('Rh(s)')  # Set inert coverage
    p.set_initial_temperature(240, 'degC')  # Set reactor initial temperature in [°C]
    p.set_resolution_method("TRANSIENT")  # Set resolution method
    p.set_gas_diffusion(False)  # Disable gas diffusion
    p.set_verbosity(False)  # Disable solver verbosity
    p.set_relative_tolerance(1.e-04)  # Set solver relative tolerance
    p.set_absolute_tolerance(1.e-04)  # Set solver absolute tolerance
    p.solve([0, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06], 's')  # Solve for different time steps in [s]
```

### **Steady State**

This [example](https://github.com/srebughini/ASALIPY/blob/main/examples/ph1d_steady_state.py) show how to solve a **steady state 1-D
pseudo-homogeneous plug flow reactor** for
the [catalytic combustion of hydrogen over rhodium](https://www.detchem.com/mechanisms).

```python
import os
from asali.reactors.ph1d_steady_state import SteadyStatePseudoHomogeneous1DReactor

if __name__ == "__main__":
    p = SteadyStatePseudoHomogeneous1DReactor(os.path.join('files', 'H2-O2-Rh.yaml'), 'gas', 'Rh_surface')
    p.set_length(2.5, 'm')  # Set reactor length in [m]
    p.set_diameter(10., 'mm')  # Set reactor diameter in [mm]
    p.set_pressure(20, 'bar')  # Set reactor pressure in [bar]
    p.set_catalytic_load(75, '1/m')  # Set catalytic load in [1/m]
    p.set_volumetric_flow_rate(10, 'm3/h')  # Set volumetric flow rate in [m3/h]
    p.set_inlet_temperature(240, 'degC')  # Set inlet gas temperature in [°C]
    p.set_inlet_mass_fraction({'O2': 0.4, 'AR': 0.5, 'H2': 0.1})  # Set inlet gas composition using mass fraction
    p.set_initial_coverage({'Rh(s)': 1})  # Set reactor initial coverage
    p.set_energy(True)  # Enable energy balance
    p.set_inert_specie('AR')  # Set inert specie
    p.set_inert_coverage('Rh(s)')  # Set inert coverage
    p.set_gas_diffusion(False)  # Disable gas diffusion
    p.set_verbosity(False)  # Disable solver verbosity
    p.solve()  # Solve
```

## **1-D Heterogeneous Plug Flow Reactor**
### **Transient**
This [example](https://github.com/srebughini/ASALIPY/blob/main/examples/het1d_transient.py) show how to solve a **transient 1-D
heterogeneous plug flow reactor** for
the [catalytic combustion of hydrogen over rhodium](https://www.detchem.com/mechanisms).

```python
import os
from asali.reactors.het1d import Heterogeneous1DReactor

if __name__ == "__main__":
    h = Heterogeneous1DReactor(os.path.join('files', 'H2-O2-Rh.xml'), 'gas', 'Rh_surface')  # Initialize reactor class
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
### **Steady State**
This [example](https://github.com/srebughini/ASALIPY/blob/main/examples/het1d_steady_state.py) show how to solve a **steady state 1-D
heterogeneous plug flow reactor** for the [catalytic combustion of hydrogen over rhodium](https://www.detchem.com/mechanisms).

```python
import os
from asali.reactors.het1d import Heterogeneous1DReactor

if __name__ == "__main__":
    h = Heterogeneous1DReactor(os.path.join('files', 'H2-O2-Rh.xml'), 'gas', 'Rh_surface')  # Initialize reactor class
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
    h.set_packed_bed_reactor(0.3, 'mm', 1, 'cm', 0.75)  # Set packed bed reactor properties
    h.solve()  # Solve
```

## **Reactor Plotter**

This [example](https://github.com/srebughini/ASALIPY/blob/main/examples/het1d_transiet.py) show how to **solve and plot** 1-D
heterogeneous plug flow reactor for the [catalytic combustion of hydrogen over rhodium](https://www.detchem.com/mechanisms).

```python
import os
from asali.reactors.het1d import Heterogeneous1DReactor
from asali.plotters.reactor import ReactorPlotter

if __name__ == "__main__":
    h = Heterogeneous1DReactor(os.path.join('files', 'H2-O2-Rh.xml'), 'gas', 'Rh_surface')
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
    h.set_resolution_method("TRANSIENT")
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

## **Cantera file converter**

This [example](https://github.com/srebughini/ASALIPY/blob/main/examples/canterafiles.py) show how to **convert** [Cantera](https://cantera.org/) file formats.

```python
import os
from asali.utils.cantera_file_converter import CanteraFileConverter

if __name__ == "__main__":
    # Convert from CHEMKIN format to YAML format
    CanteraFileConverter.from_chemkin_to_yaml(kinetic_file_path=os.path.join("files", "kinetic.kin"),
                                              thermodynamic_file_path=os.path.join("files", "thermo.tdc"),
                                              transport_file_path=os.path.join("files", "transport.tra"),
                                              surface_file_path=os.path.join("files", "surface.sur"),
                                              output_file_path=os.path.join("files", "output_v3.yaml"))
```
