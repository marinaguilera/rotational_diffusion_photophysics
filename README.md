# rotational_diffusion_photophysics
Tools for computing time-dependent fluorescence signals with rotational diffusion and complex photophysics.
An arbitrary kinetic scheme can be used and the program analytically solves the diffusion-kinetics problem. All the angular probability distribution functions are expanded with spherical harmonics. Supporting the method implementation and findings in the manuscript "Extending fluorescence anisotropy to large complexes
using reversibly switchable proteins", Andrea Volpato, Dirk Ollech, Jonatan Alvelid, Martina Damenti, Barbara Müller, Andrew York, Maria Ingaramo, Ilaria Testa (manuscript just-accepted, Nature Biotechnology, 2022).

## Content
- `rotational_diffusion_photophysics.py`  
Main engine of the kinetics and diffusion solver for an arbitrary kinetics scheme and isotorpic free rotatioal diffusion.

- `rotational_diffusion_photophysics_models.py`  
Modeling of STARSS experiments using rsEGFP2 protein.

- `plot` folder  
Some useful plotting tools for orientational probabilities and STARSS pulse schemes.

- `starss` folder  
Notebooks with computation and plotting of example STARSS simulations discussed in the notes.

- `notes` folder  
Mathematical notes about the rotational diffusion and kinetics model.