<p align="center">
  <a href="https://srebughini.github.io/ASALI/"><img src="https://github.com/srebughini/ASALI/blob/8bab03eebf389c0a6f30ececd553eed8cedf2f73/GUI/src/resources/images/BigLogo.png"/></a>
</p>
<p align="center">
  <a href="https://anaconda.org/ASALIcode/asali"><img alt="Conda" src="https://img.shields.io/conda/pn/asalicode/asali?color=orange&style=popout-square"></a>
  <a href="https://anaconda.org/asalicode/asali"><img src="https://anaconda.org/asalicode/asali/badges/version.svg" /></a>
  <a href="https://anaconda.org/ASALIcode/asali"><img alt="Conda - License" src="https://img.shields.io/conda/l/asalicode/asali?style=popout-square"></a>
  <a href="https://anaconda.org/ASALIcode/asali"><img alt="Conda" src="https://img.shields.io/conda/dn/asalicode/asali?style=popout-square"></a>
  <img alt="CodeFactor Grade" src="https://img.shields.io/codefactor/grade/github/srebughini/ASALIPY?style=flat-square">
  <a href="https://github.com/srebughini/ASALIPY/stargazers"><img src="https://img.shields.io/github/stars/srebughini/ASALIPY.svg?style=popout-square"></a>
</p>

## 1. Introduction
**ASALIPY** is part of the [ASALI](https://github.com/srebughini/ASALI) project and it is a **python** library to model **chemical reactors** based on [Cantera](https://cantera.org/). Here a list of the available reactor models:
* Batch Reactor
* Continuous Stirred Tank Reactor
* 1-Dimensional Pseudo-Homogeneous Plug Flow Reactor
* 1-Dimensional Heterogeneous Plug Flow Reactor

## 2. Installation
**ASALIPY** requires [Anaconda](https://www.anaconda.com/) as package manager, since stable versions of [Cantera](https://cantera.org/) and [Assimulo](https://jmodelica.org/assimulo/) 
are not available for [Pypi](https://pypi.org/). [Here](https://www.anaconda.com/products/individual) you can find how to install Anaconda on your operating system.  
### 2.1 Using [Anaconda](https://www.anaconda.com/)
**ASALIPY** [conda](https://www.anaconda.com/) package can be installed as follow:  
```bash
conda install -c conda-forge asalicode::asali
```  
### 2.2 Using [Github](https://github.com/srebughini/ASALIPY.git)
If you want to use **ASALIPY** locally, without installing its conda package, it can be installed as follow:  
```bash
git clone https://github.com/srebughini/ASALIPY.git
cd ASALIPY
conda env create -f environment.yml
conda activate asali
```  
## 3. Equations
Equations solved for each reactor model can be found [here](EQUATIONS.md)

## 4. Examples
Examples on how to use reactor models can be found [here](EXAMPLES.md).

## 5. Contacts
If you want to contribute, ask questions, report bugs compile the form [here](https://srebughini.github.io/ASALI/pages/contacts/) or [open an issue](https://github.com/srebughini/ASALIPY/issues).
