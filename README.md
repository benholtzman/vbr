# The Very Broadband Rheology (VBR) Calculator

**Licensing**: LICENSE?

**Bug Reporting and Support**: Join our slack channel at vbr-calc.slack.com to let us know if you find bugs or if you need help getting started!

**Citing**: Manuscript in prep. Until publication, email Holtzman and Havlin before using for publishable research.

## Overview

The Very Broadband Rheology (VBR) Calculator provides a useful framework for calculating material properties from thermodynamic state variables (e.g., temperature, pressure, melt fraction) using a wide range of experimental scalings. The main goal is to allow easy comparison between methods for calculating anelastic-dependent seismic properties, but the VBR Calculator can also be used for calculating steady state viscosity, pure elastic (anharmonic) seismic properties and more. The VBR Calculator is nominally for olivine, but may be applied to other compositions (at your own risk).

At present, the code is in Matlab, but it is functional in [GNU Octave](https://www.gnu.org/software/octave/). There are plans for a python release, pending funding.

# Basic Usage

The following outlines the basic usage for the VBR calculator. Additionally, there is a growing number of examples in  Projects/ to illustrate more complex usage, particularly in developing a statistical framework for comparing predicted mechanical properties to observed properties.  

### Initialize VBR

To start, add the top level directory to your Matlab path (relative or absolute path) and run vbr_init to add all the required directories to your path:
```
vbr_path='~/src/vbr/';
addpath(vbr_path)
vbr_init
```

### Initialize Methods List

The VBR Calculator is built around Matlab structures. All direction and data is stored in the ```VBR``` structure, which gets passed around to where it needs to go. ```VBR.in``` contains the user's input. ```VBR.out``` contains the results of any calculations.

**First**, the user must supply a cell array called ```methods_list``` for each property for which they want to calculate:
```Matlab
VBR.in.elastic.methods_list={'anharmonic';'anh_poro';};
VBR.in.viscous.methods_list={'HK2003','HZK2011'};
VBR.in.anelastic.methods_list={'eburgers_psp';'andrade_psp';'xfit_mxw'};
```

Each method will have a field in ```VBR.out```  beneath the property, e.g.,

```Matlab
VBR.out.elastic.anharmonic
VBR.out.viscous.HK2003
VBR.out.anelastic.eburgers_psp
VBR.out.anelastic.andrade_psp
```
beneath which there will be fields for the output for the calculations, e.g., ```VBR.out.anelastic.andrade_psp.Q``` for quality factor Q (attenuation=Q<sup>-1</sup>).

After VBR is initialized, a list of available methods can be printed by running `vbrListMethods()`. For theoretical background on the different methods, see the accompanying VBR Calculator Manual.

### Initialize the State Variables

The input structure ```VBR.in.SV``` contains the state variables that define the conditions at which you want to apply the methods. The following fields **MUST** be defined:

```Matlab
%  frequencies to calculate at
   VBR.in.SV.f = logspace(-2.2,-1.3,4); % [Hz]

%  size of the state variable arrays. arrays can be any shape
%  but all arays must be the same shape.
   n1 = 1;

   VBR.in.SV.P_GPa = 2 * ones(n1,1); % pressure [GPa]
   VBR.in.SV.T_K = 1473 * ones(n1,1); % temperature [K]
   VBR.in.SV.rho = 3300 * ones(n1,1); % density [kg m^-3]
   VBR.in.SV.sig_MPa = 10 * ones(n1,1); % differential stress [MPa]

   VBR.in.SV.phi = 0.0 * ones(n1,1); % melt fraction
   VBR.in.SV.dg_um = 0.01 * 1e6 * ones(n1,1); % grain size [um]

%  optional state variables
   VBR.in.SV.chi=1*ones(n1,1); % composition fraction: 1 for olivine, 0 for crust (OPTIONAL, DEFAULT 1)
   VBR.in.SV.Ch2o = 0 * ones(n1,1) ; % water concentration  (OPTIONAL, DEFAULT 0)

```

All SV arrays must be the same size and shape, except for the frequency ```VBR.in.SV.f```. They can be any length and shape as long as they are the same. Frequency dependent variables store the frequency dependencein an extra dimension of the output. If ```shape(VBR.in.SV.T)``` is (50,100) and ```numel(VBR.in.SV.f)``` is 3, then  ```shape(VBR.out.anelastic.eburgers_psp.V)``` will be (50,100,3).

### Adjust parameters (optional)

The VBR calculator allows the user to change any parameter they see fit. Parameters are stored in the VBR.in.(property).(method) structure, e.g.:

```Matlab
VBR.in.elastic.anharmonic.Gu_0_ol = 75.5; % olivine reference shear modulus [GPa]
VBR.in.viscous.HZK2011.diff.Q=350e3; % diffusion creep activation energy
```

The default parameters are stored in ```vbr/4_VBR/VBR_v0p95/params/``` and can be loaded and explored with

```Matlab
VBR.in.elastic.anharmonic=Params_Elastic('anharmonic'); % unrelaxed elasticity
VBR.in.viscous.HZK2011=Params_Viscous('HZK2011'); % HZK2011 params
```

### Run the VBR Calculator

The VBR Calculator begins calculations by passing the ```VBR``` structure to the ``VBR_spine()```:

```Matlab
[VBR] = VBR_spine(VBR) ;
```

### Pulling out results

Results are stored in ```VBR.out``` for each property type and method:

```Matlab
VBR.out.elastic.anharmonic.Vsu % unrelaxed seismic shear wave velocity
VBR.out.anelastic.eburgers_psp.V % anelastic-dependent seismic shear wave velocity
VBR.out.viscous.HZK2011.eta_total % composite steady state creep viscosity
```

### Notes for GNU Octave users

The VBR Calculator nominally works in GNU Octave, but you may find that you need to install some packages.

https://octave.sourceforge.io/io/index.html
https://octave.sourceforge.io/statistics/index.html
https://octave.org/doc/interpreter/Installing-and-Removing-Packages.html
