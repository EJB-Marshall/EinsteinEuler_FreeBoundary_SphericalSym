program evolve_einstein_euler
use real_type_mod, only: dp
use einstein_euler_mod
use grid_mod
use patch_mod
use read_parameters_mod
use einstein_euler_initial_data_mod
use simulation_mod
use diff_op_mod


implicit none
type(grid) :: base_grid
type(patch) :: base_patch
type(einstein_euler) :: einstein_euler_sys
type(simulation) :: euler_simulation
integer :: beginning, end
real(dp) :: rate
! real(dp), parameter :: pi = acos(-1.0_dp)

!Timing code
call system_clock(beginning, rate)


!#######################################
!# Load Simulation Parameters
!#######################################
call read_params()

call einstein_euler_sys%system_init()

!#######################################
!# Create Base Patch
!#######################################

call create_patch_grid(0.0_dp,5.0_dp,Nx,Ngz,base_grid) 

call create_base_patch(base_grid,einstein_euler_sys,t0,base_patch)

call initialise_euler_data(base_patch,t0)


! !#######################################
! !# Initialise Simulation
! !#######################################

! !Initialise simulation
call initialise_simulation(base_patch,0.3_dp,t0,tend,filename,euler_simulation)

! Evolve system
call euler_simulation%evolve()

call system_clock(end)
print *, "elapsed time: ", real(end - beginning) / real(rate)


end program evolve_einstein_euler




