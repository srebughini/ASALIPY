

The models.
## **1. BATCH reactor**
This example show how to solve a BATCH reactor for the [catalytic combustion of hydrogen over rhodium](https://www.detchem.com/mechanisms). 

```python
from asali.reactors.batch import BatchReactor

if __name__ == "__main__":
    # Initialize reactor class
    b = BatchReactor('H2-O2-Rh.xml', 'gas', 'Rh_surface')
    # Set reactor volume in mm3
    b.set_volume(10., 'mm3')
    # Set reactor pressure in bar
    b.set_pressure(5, 'bar')
    # Set catalytic load in 1/m
    b.set_catalytic_load(15, '1/m')
    # Set reactor initial composition using mole fraction
    b.set_initial_mole_fraction({'O2': 0.4, 'AR': 0.5, 'H2': 0.1})
    # Set reactor initial temperature in °C
    b.set_initial_temperature(120, 'degC')
    # Set reactor initial coverage
    b.set_initial_coverage({'Rh(s)': 1})
    # Enable energy balance
    b.set_energy(1)
    # Solve for different time steps
    b.solve([0, 0.1, 0.5, 5], 's')
```
