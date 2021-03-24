The models.

## **1. Batch Reactor**

This [example](https://github.com/srebughini/ASALIPY/blob/main/examples/batch.py) show how to solve a **batch reactor** for the [catalytic combustion of hydrogen over rhodium](https://www.detchem.com/mechanisms).

```python
from asali.reactors.batch import BatchReactor

if __name__ == "__main__":
    b = BatchReactor('H2-O2-Rh.xml', 'gas', 'Rh_surface')  # Initialize reactor class
    b.set_volume(10., 'mm3')  # Set reactor volume in mm3
    b.set_pressure(5, 'bar')  # Set reactor pressure in bar
    b.set_catalytic_load(15, '1/m')  # Set catalytic load in 1/m
    b.set_initial_mole_fraction({'O2': 0.4, 'AR': 0.5, 'H2': 0.1})  # Set reactor initial composition using mole fraction
    b.set_initial_temperature(120, 'degC')  # Set reactor initial temperature in °C
    b.set_initial_coverage({'Rh(s)': 1})  # Set reactor initial coverage
    b.set_energy(1)  # Enable energy balance
    b.solve([0, 0.1, 0.5, 5], 's')  # Solve for different time steps
```

## **2. Continuous Stirred Tank Reactor**

This [example](https://github.com/srebughini/ASALIPY/blob/main/examples/cstr.py) show how to solve a **continuous stirred tank reactor** for the [catalytic combustion of hydrogen over rhodium](https://www.detchem.com/mechanisms).

```python
from asali.reactors.cstr import CstrReactor

if __name__ == "__main__":
    c = CstrReactor('H2-O2-Rh.xml', 'gas', 'Rh_surface')
    c.set_volume(5., 'm3')
    c.set_pressure(5, 'bar')
    c.set_catalytic_load(150, '1/m')
    c.set_volumetric_flow_rate(1, 'm3/h')
    c.set_inlet_temperature(120, 'degC')
    c.set_inlet_mass_fraction({'O2': 0.4, 'AR': 0.5, 'H2': 0.1})
    c.set_initial_mass_fraction({'AR': 1})
    c.set_initial_temperature(25, 'degC')
    c.set_initial_coverage({'Rh(s)': 1})
    c.set_energy(1)
    c.solve([0, 0.1, 0.5, 1, 2, 5], 'min')
```