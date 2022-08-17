# TPE-Source
MatLab functions and scripts for displacement, strain and stress induced by a disk-shaped Thermo-Poro-Elastic source in a poro-elastic half-space.

This software is for academic use. Do not use this in a commercial environment. 

The functions collected in this repository implement the semi-analytical solutions for displacement, strain and stress field due to a disk-shaped Thermo-Poro-Elastic source (TPE) included in an elastic half-space, and subject to a sudden increase in pore pressure and temperature, as described in:

## Citation

[Mantiloni, L., Nespoli, M., Belardinelli, M. E., & Bonafede, M. (2020). Deformation and stress in hydrothermal regions: The case of a disk-shaped inclusion in a half-space. Journal of Volcanology and Geothermal Research, 403, 107011] (MSc. thesis, Università di Bologna) doi: https://doi.org/10.1016/j.jvolgeores.2020.107011"

The above reference should be cited whenever the software is employed.

This code was created in Windows MATLAB versions:

2017b:2019b

and tested on 2020b version.

## Index
- [Introduction](#introduction)
- [1. How to run](#1-how-to-run)
- [2. Code highlights and applications](#6-code-highlights-and-applications)  
- [3. Acknowledgements](#7-acknowledgements)  

## Introduction
The solutions provided in this software calculate the displacement, strain and stress field due to the TPE described in the Introduction, at given observation points within a half-space. The inputs and outputs are expressed in a cartesian reference frame and with the following units of measures: 

Length/Displacement - m

Pore pressure/Rigidity modulus/Stress - Pa

Temperture - K

The functions are currently optimized to work with horizontal grids of observation points (X,Y) at fixed depth (e.g. computing displacement at the free surface in an inversion of ground deformation data). 

The TPE is meant to represent a permeable rock layer stressed and strained by hot and pressurized volatiles released upward by an underlying magmatic
reservoir or, alternatively, a mushy magma storage region where an input of fresh magma is injected. The source is modeled as a thin horizontal disk inside which a sudden change of temperature (ΔT) and pore pressure (Δp) occurs. Semi-analytical solutions for the displacement and stress fields both within and outside
the source are provided. 

## How to run

[Index ^](#tpe-source)

To run:
* Download the file from GitHub
* Unzip the folder and move it to the directory you would like to work from
* Open MATLAB
* Add the folder to the MatLab search path
* Run using and editing 'MainFrame.m'

### Functions





