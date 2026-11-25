module patch_mod
    use real_type_mod, only: dp
    use grid_mod
    use timelevel_mod
    use pde_system_mod
    implicit none

    !!################################ PATCH MOD ###############################################
    !!                       
    !! A patch links the grid, solution data, and PDE system together.
    !! It is essentially a grid with solution data. 
    !!
    !! cgrid - The computational grid \
    !! p_num - The patch number (w.r.t Ronchi et al. patch numbering) \
    !! tl_current - The current timelevel (i.e. solution data) \
    !! tl_previous - The previous timelevel (i.e. solution data) \
    !! system - The PDE system which is being solved \
    !!#########################################################################################


    type :: patch

        type(grid) ::  cgrid
        real(dp) :: t 
        real(dp) :: dt
        real(dp) :: x_boundary(2) ! Locations of the free boundary
        type(timelevel), allocatable :: tl_current 
        real(dp), allocatable :: local_error(:)
        class(pde_system), allocatable :: system
        integer :: iterations
    

    end type patch

    contains 
    
    !#----------------------------------------------
    !# Create Base Patch
    !#---------------------------------------------- 

    subroutine create_base_patch(new_grid,sys,t0,new_patch)
        implicit none
        type(grid), intent(in) :: new_grid
        class(pde_system), intent(in) :: sys
        real(dp), intent(in) :: t0
        integer :: Nx, Ngz, Nvar
        type(patch), intent(inout) :: new_patch

        ! Grid size and number of variables
        Nx = new_grid%Nx
        Ngz = new_grid%Ngz
        Nvar = sys%Nvar

        ! Set initial time + timestep
        new_patch%t = t0
        new_patch%dt = 0.0_dp

        ! Set the grid
        new_patch%cgrid = new_grid

        ! Set iterations
        new_patch%iterations = 0

        ! Set initial error
        allocate(new_patch%local_error(1-Ngz:Nx+Ngz))
        new_patch%local_error = 0.0_dp
        
        ! Allocate the PDE system
        if (allocated(new_patch%system)) deallocate(new_patch%system)
        allocate(new_patch%system, source=sys)
        
        
        ! Allocate time level objects
        allocate(new_patch%tl_current)
        
        ! Allocate solution arrays with ghost zones
        allocate(new_patch%tl_current%soln(1-Ngz:Nx+Ngz, Nvar))

        
        ! Initialise arrays
        new_patch%tl_current%soln = 0.0_dp
        new_patch%tl_current%t = t0
        new_patch%x_boundary = 0.0_dp


    end subroutine create_base_patch
    
        

end module patch_mod
