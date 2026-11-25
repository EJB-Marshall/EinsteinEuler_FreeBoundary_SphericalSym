module pde_system_mod
    use real_type_mod, only: dp
    use grid_mod
    implicit none


    !!################################ PDE_SYSTEM MOD ########################################
    !!                       
    !! This provides an abstract base class for PDE systems. \
    !! When creating a specific PDE system, we will *extend* this type.
    !! 
    !! 
    !! 
    !!#########################################################################################

    type, abstract :: pde_system 

    integer :: Nvar 
    real(dp) :: CS_x 
    real(dp), allocatable :: flux(:,:)
    real(dp), allocatable :: lapse_tmp(:)

    contains

    ! Deferred Initialisation Procedure - This is how we will set system dependent parameters, e.g. Nvar etc.
    procedure(system_init_interface), deferred :: system_init

    ! Deferred procedures for computing RHS, boundary conditions etc.
    procedure(compute_rhs_interface), deferred :: compute_rhs
    procedure(compute_boundary_rhs_interface), deferred :: compute_boundary_rhs
    procedure(apply_BCs_interface), deferred :: apply_BCs
    procedure(apply_boundary_BCs_interface), deferred :: apply_boundary_BCs

    end type pde_system

    ! Interfaces for PDE functions
    abstract interface 

    subroutine system_init_interface(this)
        import pde_system
        implicit none
        class(pde_system), intent(inout) ::  this
    
    end subroutine system_init_interface


    subroutine compute_rhs_interface(this,t,u,cgrid,dudt)
        import :: dp, pde_system, grid
        implicit none
        class(pde_system), intent(inout) :: this
        type(grid), intent(in) :: cgrid
        real(dp), intent(in) :: u(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,this%Nvar)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: dudt(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,this%Nvar)

    end subroutine compute_rhs_interface


    subroutine compute_boundary_rhs_interface(this,t,u,x_boundary,cgrid,dudt)
        import :: dp, pde_system, grid
        implicit none
        class(pde_system), intent(inout) :: this
        type(grid), intent(in) :: cgrid
        real(dp), intent(in) :: u(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,this%Nvar)
        real(dp), intent(in) :: x_boundary(2)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: dudt(2)

    end subroutine compute_boundary_rhs_interface


    subroutine apply_BCs_interface(this,cgrid,u)
        import :: dp, pde_system, grid
        implicit none
        class(pde_system), intent(in) :: this
        type(grid), intent(in) :: cgrid
        real(dp), intent(inout) :: u(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,this%Nvar)

    end subroutine apply_BCs_interface


    subroutine apply_boundary_BCs_interface(this,cgrid,t,u)
        import :: dp, pde_system, grid
        implicit none
        class(pde_system), intent(inout) :: this
        type(grid), intent(in) :: cgrid
        real(dp), intent(in) :: t
        real(dp), intent(inout) :: u(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,this%Nvar)

    end subroutine apply_boundary_BCs_interface
    

    end interface


end module pde_system_mod