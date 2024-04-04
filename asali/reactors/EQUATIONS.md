# Equations
In the following the list of the equations solved for each reactor model and a table explaining the symbols:
|Symbol|Meaning|Unit dimension|
|-|-|-|
|$i$|Gas specie index|n.a.|
|$j$|Coverage specie index|n.a.|
|$&alpha;$|Catalytic load|$\frac{1}{m}$|
|$&omega;$|Gas mass fraction|-|
|$&omega;^0$|Gas mass fraction at initial conditions|-|
|$&Gamma;$|Site density|$\frac{kmol}{m^2}$|
|$\dot{m}$|Inlet mass flow rate|$\frac{kg}{m^3s}$|
|$m$|Total mass|$kg$|
|$c_{p}$|Specific heat|$\frac{J}{kgK}$|
|$t$|Time|$s$|
|$T$|Temperature|$K$|
|$V$|Reactor volume|$m^3$|
|$MW_{i}$|Gas specie molecular weight|$\frac{kg}{kmol}$|
|$Q^{hom}$|Heat of reaction from homogeneous reactions|$\frac{J}{m^3s}$|
|$Q^{het}$|Heat of reaction from heterogeneous reactions|$\frac{J}{m^2s}$|
|${R_{i}}^{hom}$|Gas specie reaction rate from homogeneous reactions|$\frac{kmol}{m^3s}$|
|${R_{i}}^{het}$|Gas specie reaction rate from heterogeneous reactions|$\frac{kmol}{m^2s}$|



## **Batch Reactor**
### Mass balance
```math
\frac{\partial &omega;_{i}}{\partial t} = \frac{MW_{i}{R_{i}}^{hom}}{&rho;} + \frac{&alpha;MW_{i}{R_{i}}^{het}}{&rho;} - \frac{&omega;_i}{m}\frac{\partial m}{\partial t}
```
### Total mass balance
```math
\frac{\partial m}{\partial t} = V&alfa;\sum{MW_{i}R_{i}}^{het}}
```
### Coverage balance
```math
\frac{\partial &theta;_{j}}{\partial t} = \frac{{R_{j}}^{het}}{&Gamma;}
```
### Energy balance
```math
\frac{\partial T}{\partial t} = \\frac{Q^{hom}}{&rho;c_{p}} + \frac{&alpha;Q^{het}}{&rho;c_{p}}
```


## **CSTR Reactor**

### Mass balance
```math
\frac{\partial &omega;_{i}}{\partial t} = \frac{\dot{m}(&omega;^0 - &omega;_i )}{V&rho;} + \frac{MW_{i}{R_{i}}^{hom}}{&rho;} + \frac{&alpha;MW_{i}{R_{i}}^{het}}{&rho;}
```
### Coverage balance
```math
\frac{\partial &theta;_{j}}{\partial t} = \frac{{R_{j}}^{het}}{&Gamma;}
```
### Energy balance
```math
\frac{\partial T}{\partial t} = \frac{\dot{m}(T^0 - T)}{V&rho;} + \frac{Q^{hom}}{&rho;c_{p}} + \frac{&alpha;Q^{het}}{&rho;c_{p}}
```
