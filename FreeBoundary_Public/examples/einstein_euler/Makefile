FC = gfortran
FFLAGS = -Wall -Wextra -Wuninitialized -O0 -fbounds-check -I/opt/homebrew/include -fopenmp -g \
		 -finit-real=snan -ffpe-trap=invalid,zero,overflow \
         -fbacktrace
# FFLAGS = -Wall -Wextra -O3 -fbounds-check -fopenmp 
SRC_DIR = ../../src
LIB_DIR = ../../src
LIB_NAME = libPDE_modules.a
LDFLAGS = -fopenmp -L/opt/homebrew/lib \
		-lhdf5_fortran -lhdf5_f90cstub -lhdf5 

# Include path for modules
INCLUDE_FLAGS = -I$(SRC_DIR)

# Library linking
LIB_FLAGS = -L$(LIB_DIR) -lPDE_modules



#-------------------------
# Local SRC and Objects 
#-------------------------
SRC = einstein_euler_initial_data_mod.f90 \
einstein_euler_mod.f90 \
evolve_einstein_euler.f90 \

# Object files
OBJ = $(SRC:.f90=.o)


#-------------------------
# Compilation
#-------------------------

.PHONY: clean


# Build executable
main: setup $(OBJ) $(LIB_DIR)/$(LIB_NAME)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $(LIB_FLAGS) $(LDFLAGS)
	@echo "Executable created successfully"

# Check AMR Library is (re)built
$(LIB_DIR)/$(LIB_NAME): force-src
	$(MAKE) -C $(SRC_DIR)

# Force check of SRC (make will decide if rebuild is needed)
force-src:

%.o : %.f90
	$(FC) $(FFLAGS) $(INCLUDE_FLAGS)  -o $@ -c $<


# Create HDF_Files directory if it doesn't exist
setup:
	@mkdir -p HDF_Files

# Clean up object/mod/executable files etc. 

clean:
	@rm -f *.o *.mod 

#-------------------------
# Dependencies 
#-------------------------
einstein_euler_initial_data_mod.o: $(LIB_DIR)/$(LIB_NAME)
einstein_euler_mod.o: $(LIB_DIR)/$(LIB_NAME)
evolve_einstein_euler.o: einstein_euler_mod.o einstein_euler_initial_data_mod.o $(LIB_DIR)/$(LIB_NAME)


