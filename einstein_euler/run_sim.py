import subprocess
import numpy as np
import argparse
import sys
from TOV_ID import Construct_Initial_Data, Construct_Perturbed_Initial_Data



#######################################
#Parser Settings
#######################################

# Initialise parser
parser = argparse.ArgumentParser(description=\
"""This program is used to initialise parameters and run simulations. """)

# Parse files
parser.add_argument('-N', type=int,help=\
"""The number of grid points.""")
parser.add_argument('-f','-file', help=\
"""The name of the hdf file to be produced.""")
parser.add_argument('-K',type=float, help=\
"""The EOS parameter K""")
args = parser.parse_args()


#Check Inputs
if args.N is None:
    print("Error: argument -N is required.")
    sys.exit(1)
if args.f is None:
    print("Error: argument -f is required.")
    sys.exit(1)

### Define Grid/Simulation/Aux Parameters

nml_1 = "&Grid_Params"
Nx = args.N
Ngz = 3


nml_2 = "&Sim_Params"
t0 = 0.0
tend = 100
filename = args.f

nml_3 = "&Aux_Params"
aux1 = args.K
aux2 = 0.0
aux3 = 0.0


### TOV Initial Data
# Construct_Initial_Data(0,5,Nx,Ngz,1.1,args.K)

### Perturbed TOV Initial Data
star_bdy_r = Construct_Perturbed_Initial_Data(0,5.0,Nx,Ngz,1.1,args.K,5e-3)



### Create namelist file
f = open("params.nml","w")

f.write(nml_1+"\n")
f.write("Nx = " + str(Nx) + "\n")
f.write("Ngz = " + str(Ngz) + "\n")
f.write("/\n")
f.write("\n")
f.write(nml_2+"\n")
f.write(f"t0 = {t0:.15e}\n")
f.write(f"tend = {tend:.15e}\n")
f.write("filename = '" + str(filename) + "'" + "\n")
f.write("/\n")
f.write(nml_3+"\n")
f.write(f"aux1 = {aux1:.15e}\n")
f.write(f"aux2 = {aux2:.15e}\n")
f.write(f"aux3 = {aux3:.15e}\n")
f.write("/\n")
f.close() # Need to close namelist file before running Fortran!


### Run our Fortran script
subprocess.run("./main")
