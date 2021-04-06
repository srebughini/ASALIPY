<p align="center">
  <a href="https://imgur.com/T5W1buF"><img src="https://i.imgur.com/T5W1buFm.png" title="source: imgur.com" /></a>
</p>
<p align="center">
  <a href="https://github.com/srebughini/JASALI/blob/main/LICENSE"><img alt="GitHub" src="https://img.shields.io/github/license/srebughini/ASALIPY"></a>
  <a href="https://www.codefactor.io/repository/github/srebughini/asalipy"><img src="https://www.codefactor.io/repository/github/srebughini/asalipy/badge" alt="CodeFactor" /></a>
</p>

## 1. Introduction
**ASALIPY** is part of the [ASALI](https://github.com/srebughini/ASALI) project and it is a **python** library to model **chemical reactors** based on [Cantera](https://cantera.org/). Here a list of the available reactor models:
* Batch Reactor
* Continuous Stirred Tank Reactor
* 1-Dimensional Pseudo-Homogeneous Plug Flow Reactor
* 1-Dimensional Heterogeneous Plug Flow Reactor

## 2. Installation
**ASALIPY** requires [Anaconda](https://www.anaconda.com/) as package manager, since stable versions of [Cantera](https://cantera.org/) and [Assimulo](https://jmodelica.org/assimulo/) 
are not for [Pypi](https://pypi.org/). [Here](https://www.anaconda.com/products/individual) you can find how to install Anaconda on your operating system.  
### 2.1 Using [Anaconda](https://www.anaconda.com/)
**ASALIPY** [conda](https://www.anaconda.com/) package can be installed as follow:  
```bash
conda install -c conda-forge asali #STILL WORKING ON IT
```  
### 2.2 Using [Github](https://github.com/srebughini/ASALIPY.git)
If you want to use **ASALIPY** locally, without installing its conda package, it can be installed as follow:  
```bash
git clone https://github.com/srebughini/ASALIPY.git
cd ASALIPY
conda env create -f environment.yml
conda activate asali
```  

## 2. Examples
Examples on how to use reactor models can be found [here](https://github.com/srebughini/ASALIPY/tree/main/examples).

## 3. Contacts
If you want to contribute, ask questions, report bugs compile the form [here](https://srebughini.github.io/ASALI/pages/contacts/) or [open an issue](https://github.com/srebughini/ASALIPY/issues).
