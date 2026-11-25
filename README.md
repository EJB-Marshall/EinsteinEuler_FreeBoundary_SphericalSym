# EinsteinEuler_FreeBoundary_SphericalSym
A code to numerically evolve the Einstein-Euler equations with a free boundary in spherical symmetry.

The code requires the HDF5 library and the following python modules:
- NumPy
- SciPy
- h5py
- matplotlib

## Usage
After compiling the fortran code, a simulation can be run using run_sim.py which takes the following command line arguments: 

- N - Number of grid points/cells
- K - The sound speed parameter
- f - The filename to be saved (must include .hdf5 extension)

For example, a simulation with 1000 grid points and K=1 can be run using 
```
python3 -O run_sim.py -N 1000 -K 1.0 -f output.hdf5 
```
The results of the simulation are saved to the file output.hdf5 in the folder /examples/einstein_euler/HDF_Files.


Further parameters, e.g. perturbation size, domain size, CFL number etc. can be modified in run_sim.py and evolve_einstein_euler.f90


The boundaries of the fluid domain are explicitly tracked and boundary conditions are enforced via the inverse Lax-Wendroff procedure of Tan, Shu, and collaborators, see for example: https://www.sciencedirect.com/science/article/pii/S002199911100235X, https://arxiv.org/abs/2010.05126

