# EinsteinEuler_FreeBoundary_SphericalSym
A code to numerically evolve the Einstein-Euler equations with a free boundary in spherical symmetry.

The code requires HDF5 and the following python modules:
- NumPy
- SciPy
- h5py
- matplotlib

After compiling the fortran code, a simulation can be run using run_sim.py which takes the following command line arguments:
-N - Number of grid points/cells
-K - The sound speed parameter
-f - The filename to be saved (must include .hdf5 extension)
Further parameters, e.g. perturbation size, domain size, CFL number etc. can be modified in run_sim.py and evolve_einstein_euler.f90

