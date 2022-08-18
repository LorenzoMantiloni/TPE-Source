# TPE-Source
MatLab functions and scripts for displacement, strain and stress induced by a disk-shaped Thermo-Poro-Elastic source in a poro-elastic half-space.

This software is for academic use. Do not use this in a commercial environment. 

The functions collected in this repository implement the semi-analytical solutions for displacement, strain and stress field due to a disk-shaped Thermo-Poro-Elastic source (TPE) included in an elastic half-space, and subject to a sudden increase in pore pressure and temperature, as described in:

## Citation

[Mantiloni, L., Nespoli, M., Belardinelli, M. E., & Bonafede, M. (2020). Deformation and stress in hydrothermal regions: The case of a disk-shaped inclusion in a half-space. Journal of Volcanology and Geothermal Research, 403, 107011] (MSc. thesis, Università di Bologna) doi: https://doi.org/10.1016/j.jvolgeores.2020.107011"

The above reference should be cited whenever the software is employed. It is referred to as "Mantiloni et al., 2020" in the following and in the scripts.

This code was created in Windows MATLAB versions:

2017b:2019b

and tested on 2020b version.

## Index
- [Introduction](#introduction)
- [1. How to run](#1-how-to-run)
- [2. Acknowledgements](#7-acknowledgements)  

## Introduction
The solutions provided in this software calculate the displacement, strain and stress field due to the TPE described in the Introduction, at given observation points within a half-space. The inputs and outputs are expressed in a cartesian reference frame and with the following units of measures: 

Length/Displacement - m

Pore pressure/Rigidity modulus/Stress - Pa

Temperture - K

Please note that the z-axis is positive downward in Mantiloni et al., 2020; nevertheless, the vertical component of the displacement field provided by the software is positive for uplift, and negative for subsidence.

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
* Use the functions "TPE_Displacement" and "TPE_Stresses" as illustrated in the usage examples provided in their respective descriptions.

### Functions

#### TPE_Displacement

Thermo-Poro-Elastic (Source) Displacement. Computes the displaceent components at points (x_obs, y_obs, z_obs) in a half-space, due to a TPE centred at point (tpe_x, tpe_y, c).          

Arguments: (input)

     x_obs, y_obs, z_obs  - The observation point location in cartesian coordinates. 
     They should be provided either as a NxN grid, or as a
     one-dimensional array (e.g. x_obs=[x_obs(1)...x_obs(N)],y_obs=0)

     tpe_x, tpe_y - The horizontal cartesian coordinates of the TPE
     centre.

     c - The depth of the TPE centre with respect to the free surface
     (must be positive).

     a - The radius of the TPE

     d - The thickness of the TPE

     Delta_p - The increase in pore pressure within the TPE (must be
     positive)

     Delta_T - The increase in temperature within the TPE (must be
     positive)

     H - Biot's constant for the medium 
     (see Mantiloni et al., 2020, Equation 3)
        
     alpha - coefficient of thermal expansion for the medium 
     (see Mantiloni et al., 2020, Equation 3)
        
     nu - Poisson's ratio for the medium 
        
     mu - Modulus of rigidity for the medium 

     varargin: % n - The number of coefficients in the Legendre's 
     polynomials expansion for singular displacement components (default is 
     200), see Mantiloni et al., 2020, Equation 12.
    
Arguments: (output)

     U,V,W - The displacement components in a cartesian reference frame at point 
     (x_obs,y_obs,z_obs). U is the x-component, V the y-component, W the
     vertical component. For NxN obs points, there will be three Nx1
     arrays.

#### TPE_Stresses

Thermo-Poro-Elastic (Source) Stresses. Computes the components of the stress tensor at points (x_obs, y_obs, z_obs) in a half-space, due to a TPE centred at point (tpe_x, tpe_y, c).

Arguments: (input)

     x_obs, y_obs, z_obs  - The observation point location in cartesian coordinates. 
     They should be provided either as a NxN grid, or as a
     one-dimensional array (e.g. x_obs=[x_obs(1)...x_obs(N)],y_obs=0)

     tpe_x, tpe_y - The horizontal cartesian coordinates of the TPE
     centre.

     c - The depth of the TPE centre with respect to the free surface
     (must be positive).

     a - The radius of the TPE

     d - The thickness of the TPE

     Delta_p - The increase in pore pressure within the TPE (must be
     positive)

     Delta_T - The increase in temperature within the TPE (must be
     positive)

     H - Biot's constant for the medium 
     (see Mantiloni et al., 2020, Equation 3)
        
     alpha - coefficient of thermal expansion for the medium 
     (see Mantiloni et al., 2020, Equation 3)
        
     nu - Poisson's ratio for the medium 
        
     mu - Modulus of rigidity for the medium 

     varargin: % n - The number of coefficients in the Legendre's 
     polynomials expansion for singular displacement components (default is 
     200), see Mantiloni et al., 2020, Equation 12.
    
Arguments: (output)

     Stress - The stress tensor in a cartesian reference frame at point 
     (x_obs,y_obs,z_obs), given as a one-row array for each obs point
     [T_xx,T_yy,T_zz,T_xy,T_xz,T_yz] (for N obs points, it will populate a
     N x 6 matrix)

     Strain - The strain tensor in a cartesian reference frame at point 
     (x_obs,y_obs,z_obs) given as a one-row array for each os point
     [E_xx,E_yy,E_zz,E_xy,E_xz,E_yz] (for N obs points, it will populate a
     N x 6 matrix)
     
# Acknowledgments

This software was developed under the supervision of Dr. Massimo Nespoli, Prof. Dr. Maria Elina Belardinelli and Prof. Dr. Maurizio Bonafede at Università di Bologna, as part of a M.Sc. thesis, in the years 2018-2019. Valuable comments and suggestions have also been provided by Tim Davis, Mehdi Nikkhoo, Eleonora Rivalta and Francesca Silverii.
