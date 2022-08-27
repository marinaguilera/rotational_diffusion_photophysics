# rotational_diffusion_photophysics
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7009614.svg)](https://doi.org/10.5281/zenodo.7009614)  
Tools for computing time-dependent fluorescence signals with rotational diffusion and complex photophysics.  
An arbitrary kinetic scheme can be used and the program analytically solves the diffusion-kinetics problem. All the angular probability distribution functions are expanded with spherical harmonics.  

## Content
- `rotational_diffusion_photophysics.py`  
Main engine of the kinetics and diffusion solver for an arbitrary kinetics scheme and isotorpic free rotatioal diffusion.

- `rotational_diffusion_photophysics_models.py`  
Modeling of STARSS experiments using rsEGFP2 protein.

- `plot` folder  
Some useful plotting tools for orientational probabilities and STARSS pulse schemes.

- `starss` folder  
Notebooks with computation and plotting of example STARSS simulations.

- `notes` folder  
Mathematical notes about the rotational diffusion and kinetics model.
