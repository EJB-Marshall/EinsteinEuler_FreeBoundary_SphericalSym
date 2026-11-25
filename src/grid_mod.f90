module grid_mod
    use real_type_mod, only: dp
    implicit none

    !!################################ GRID MOD ######################################################
    !!                       
    !! A grid contains the array of cell centres
    !! 
    !! 
    !! x_min - The location of the left boundary \
    !! x_max - The location of the right boundary \
    !! dx - The spatial step size \
    !! Nx - Number of cell centres in the x-direction \
    !! Ngz - Number of ghost points at each boundary (same in x and y directions) \
    !! coordinates - The coordinates of the grid *including* ghost points, indexing (1-Ngz:Nx+Ngz) \
    !! 
    !!
    !!##############################################################################################


    type :: grid

        real(dp) :: x_min, x_max, dx
        integer ::  Nx, Ngz
        real(dp), allocatable, dimension(:) :: coordinates_x ! Computational coordinates

    end type grid 


    contains

    !!##############################################
    !! This subroutine sets up a single patch's 
    !! grid in *computational* coordinates
    !!
    !!##############################################

    subroutine create_patch_grid(x_min,x_max,Nx,Ngz,patch_grid)
        implicit none
        real(dp), intent(in) :: x_min, x_max 
        integer, intent(in) :: Nx, Ngz
        type(grid), intent(out) :: patch_grid
        integer ::  i

        patch_grid%x_min = x_min
        patch_grid%x_max = x_max
        patch_grid%Nx = Nx
        patch_grid%Ngz = Ngz
        
        ! Compute spatial step size 
        patch_grid%dx = (x_max-x_min)/(real(Nx,dp))

        ! Create grid meshes
        allocate(patch_grid%coordinates_x(1-Ngz:Nx+Ngz))

        do i = 1-Ngz, Nx + Ngz
            patch_grid%coordinates_x(i) = x_min + (real(i,dp) - 0.5_dp)*patch_grid%dx ! Compute cell centres in x direction
        end do

        
    end subroutine create_patch_grid

end module grid_mod