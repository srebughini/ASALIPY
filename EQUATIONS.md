# Equations
In the following the list of the equations solved for each reactor model



## **Batch Reactor**

#### Mass balance
```math
\frac{\partial &omega;_{i}}{\partial t} = \frac{MW_{i}{R_{i}}^{hom}}{&rho;} + \frac{&alpha;MW_{i}{R_{i}}^{het}}{&rho;} - \frac{&omega;_i}{m}\frac{\partial m}{\partial t}
```

#### Total mass balance
```math
\frac{\partial m}{\partial t} = &alpha;V\sum_{i}^{NS}{MW_{i}R_{i}^{het}}
```
#### Coverage balance
```math
\frac{\partial &theta;_{j}}{\partial t} = \frac{{R_{j}}^{het}}{&Gamma;}
```
#### Energy balance
```math
\frac{\partial T}{\partial t} = \frac{Q^{hom}}{&rho;c_{p}} + \frac{&alpha;Q^{het}}{&rho;c_{p}}
```


## **CSTR Reactor**

#### Mass balance
```math
\frac{\partial &omega;_{i}}{\partial t} = \frac{\dot{m}(&omega;^0 - &omega;_i )}{V&rho;} + \frac{MW_{i}{R_{i}}^{hom}}{&rho;} + \frac{&alpha;MW_{i}{R_{i}}^{het}}{&rho;}
```
#### Coverage balance
```math
\frac{\partial &theta;_{j}}{\partial t} = \frac{{R_{j}}^{het}}{&Gamma;}
```
#### Energy balance
```math
\frac{\partial T}{\partial t} = \frac{\dot{m}(T^0 - T)}{V&rho;} + \frac{Q^{hom}}{&rho;c_{p}} + \frac{&alpha;Q^{het}}{&rho;c_{p}}
```

## **1-D Pseudo-Homogeneous Plug Flow Reactor: *Steady-State***
#### Mass balance
```math
\frac{\partial &omega;_{i}}{\partial z} = \frac{MW_{i}A({R_{i}}^{hom} + &alpha;{R_{i}}^{het})}{\dot{m}}
```
#### Coverage balance
```math
0 = \frac{{R_{j}}^{het}}{&Gamma;}
```
#### Energy balance
```math
\frac{\partial T}{\partial z} = \frac{A(Q^{hom} + &alpha;Q^{het})}{\dot{m}c_{p}}
```

## **1-D Pseudo-Homogeneous Plug Flow Reactor: *Transient***
#### Mass balance
```math
\frac{\partial &omega;_{i}}{\partial t} = -\frac{\dot{m}}{A}\frac{\partial &omega;_{i}}{\partial z} + D^{mix}_{i}\frac{\partial^2 &omega;_{i}}{\partial^2 z} + \frac{MW_{i}{R_{i}}^{hom}}{&rho;} + \frac{&alpha;MW_{i}{R_{i}}^{het}}{&rho;}
```
#### Coverage balance
```math
\frac{\partial &theta;_{j}}{\partial t} = \frac{{R_{j}}^{het}}{&Gamma;}
```
#### Energy balance
```math
\frac{\partial T}{\partial t} = -\frac{\dot{m}}{A&rho;}\frac{\partial T}{\partial z} +  \frac{k^{gas}_{mix}}{&rho;c_{p}}\frac{\partial^2 T}{\partial^2 z} + \frac{Q^{hom}}{&rho;c_{p}} + \frac{&alpha;Q^{het}}{&rho;c_{p}}
```

## **1-D Heterogeneous Plug Flow Reactor: *Steady-State***
#### Mass balance
```math
\frac{\dot{m}}{A&rho;}\frac{\partial &omega;_{i}}{\partial z} = + D^{mix}_{i}\frac{\partial^2 &omega;_{i}}{\partial^2 z} -\frac{A_{s}K_{mat}}{&epsi;}(&omega;_{i} - &omega;^S_{i})+ \frac{MW_{i}{R_{i}}^{hom}}{&rho;}
```
#### Solid mass balance
```math
0 = A_{s}K_{mat}&rho;&epsi;(&omega;_{i} - &omega;^S_{i}) + &epsi;&alpha;MW_{i}{R_{i}}^{het}
```
#### Coverage balance
```math
0 = \frac{{R_{j}}^{het}}{&Gamma;}
```
#### Energy balance

\frac{\partial T}{\partial t} =  + \frac{Q^{hom}}{&rho;c_{p}} + \frac{&alpha;Q^{het}}{&rho;c_{p}}


#### Solid energy balance

## **1-D Heterogeneous Plug Flow Reactor: *Transient***
#### Mass balance
```math
\frac{1partial &omega;_{i}}{dt} = -\frac{\dot{m}}{A&rho;}\frac{\partial &omega;_{i}}{\partial z}  + D^{mix}_{i}\frac{\partial^2 &omega;_{i}}{\partial^2 z} -\frac{A_{s}K_{mat}}{&epsi;}(&omega;_{i} - &omega;^S_{i})+ \frac{MW_{i}{R_{i}}^{hom}}{&rho;}
```
#### Solid mass balance
```math
0 = A_{s}K_{mat}&rho;&epsi;(&omega;_{i} - &omega;^S_{i}) + &epsi;&alpha;MW_{i}{R_{i}}^{het}
```
#### Coverage balance
```math
\frac{\partial &theta;_{j}}{\partial t} = \frac{{R_{j}}^{het}}{&Gamma;}
```
#### Energy balance
```math
\frac{\partial T}{\partial t} = -\frac{\dot{m}}{A&rho;}\frac{\partial T}{\partial z} +  \frac{k^{gas}_{mix}}{&rho;c_{p}}\frac{\partial^2 T}{\partial^2 z} + \frac{Q^{hom}}{&rho;c_{p}} 
```

#### Solid energy balance
```math
\frac{\partial T^S}{\partial t} = \frac{k^S}{&rho;^Sc_{p}^S}\frac{\partial^2 T^S}{\partial^2 z} + \frac{&alpha;Q^{het}}{&rho;^Sc_{p}^S(1-&epsi;)} + \frac{A_{s}h(T - T^s)}{&rho;^Sc_{p}^S(1-&epsi;)}
```


## Symbols
Here is the symbols meaning:
|Symbol|Meaning|Unit dimension|
|-|-|-|
|$i$|Gas specie index|$-$|
|$j$|Coverage specie index|$-$|
|$&alpha;$|Catalytic load|$\frac{1}{m}$|
|$&epsi;$|Reactor void fraction|$-$|
|$&omega;$|Gas mass fraction|$-$|
|$&omega;^0$|Gas mass fraction at initial conditions|$-$|
|$&omega;^S$|Gas mass fraction in the solid phase|$-$|
|&theta;|Coverage fraction|$-$|
|$&Gamma;$|Site density|$\frac{kmol}{m^2}$|
|$A$|Reactor cross section area|$m^2$|
|$A_{s}$|Reactor specific area|$\frac{1}{m}$|
|$c_{p}$|Specific heat|$\frac{J}{kgK}$|
|$D^{mix}$|Mixture diffusion coefficient|$\frac{m^2}{s}$|
|$k^{gas}_{mix}$|Mixture thermal conductivity|$\frac{W}{mK}$|
|$K_{mat}$|Gas-to-solid mass transfer coefficient|$\frac{m}{s}$|
|$m$|Total mass|$kg$|
|$\dot{m}$|Inlet mass flow rate|$\frac{kg}{m^3s}$|
|$MW_{i}$|Gas specie molecular weight|$\frac{kg}{kmol}$|
|$Q^{hom}$|Heat of reaction from homogeneous reactions|$\frac{J}{m^3s}$|
|$Q^{het}$|Heat of reaction from heterogeneous reactions|$\frac{J}{m^2s}$|
|${R_{i}}^{hom}$|Gas specie reaction rate from homogeneous reactions|$\frac{kmol}{m^3s}$|
|${R_{i}}^{het}$|Gas specie reaction rate from heterogeneous reactions|$\frac{kmol}{m^2s}$|
|${R_{j}}^{het}$|Coverage reaction rate|$\frac{kmol}{m^2s}$|
|$t$|Time|$s$|
|$T$|Temperature|$K$|
|$V$|Reactor volume|$m^3$|
|$z$|Reactor lenght|$m$|
