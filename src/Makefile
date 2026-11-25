# Makefile for src

FC = gfortran
FFLAGS = -Wall -Wextra -O3 -fbounds-check -I/opt/homebrew/include -fopenmp  -g \
		-finit-real=snan -ffpe-trap=invalid,zero,overflow \
         -fbacktrace
# LDFLAGS = -fopenmp -L/opt/homebrew/lib \
# 		-lhdf5_fortran -lhdf5_f90cstub -lhdf5
AR = ar
ARFLAGS = rcs

LIB = libPDE_modules.a

#-------------------------
# Source Files
#-------------------------


# Source files (in dependency order)
SRC = real_type_mod.f90 \
read_parameters_mod.f90 \
grid_mod.f90 \
diff_op_mod.f90 \
timelevel_mod.f90 \
pde_system_mod.f90 \
patch_mod.f90 \
timestepping_mod.f90 \
compute_norms_mod.f90 \
hdf5_IO_mod.f90 \
simulation_mod.f90 


# Object files
OBJ = $(SRC:.f90=.o)

# Default target:
all: $(LIB)

# Build the library
$(LIB): $(OBJ)
	$(AR) $(ARFLAGS) $@ $^
	@echo "Common library $(LIB) created successfully"

%.o : %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

.PHONY: clean

clean:
	@rm -f *.o *.mod *.a


#-------------------------
# Dependencies 
#-------------------------

read_parameters_mod.o: real_type_mod.o
grid_mod.o: 
diff_op_mod.o: real_type_mod.o grid_mod.o
timelevel_mod.o: real_type_mod.o
pde_system_mod.o: real_type_mod.o grid_mod.o
patch_mod.o: real_type_mod.o grid_mod.o timelevel_mod.o pde_system_mod.o 
timestepping_mod.o: real_type_mod.o read_parameters_mod.o patch_mod.o 
compute_norms_mod.90 : real_type_mod.o read_parameters_mod.o patch_mod.o
compute_norms_mod.o: real_type_mod.o patch_mod.o diff_op_mod.o 
hdf5_IO_mod.o: real_type_mod.o patch_mod.o 
simulation_mod.o: real_type_mod.o patch_mod.o timestepping_mod.o hdf5_IO_mod.o 









