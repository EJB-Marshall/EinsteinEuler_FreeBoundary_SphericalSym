module einstein_euler_mod
    use real_type_mod, only: dp
    use pde_system_mod
    use grid_mod
    use read_parameters_mod, only: Nx, Ngz, aux1, aux2, aux3
    use diff_op_mod
    implicit none
    real(dp), parameter :: pi = acos(-1.0_dp)


    type, extends(pde_system) :: einstein_euler 

    real(dp) :: K
    real(dp), allocatable :: pressure(:)
    real(dp) :: x_left, x_right ! Values of the left and right boundaries
    integer :: i_left, i_right ! Indices of the left and right boundaries


    contains 

    procedure :: system_init => einstein_euler_init
    procedure :: compute_rhs => einstein_euler_compute_rhs    
    procedure :: compute_boundary_rhs => einstein_euler_compute_boundary_rhs   
    procedure :: apply_BCs => einstein_euler_apply_BCs
    procedure :: apply_boundary_BCs => einstein_euler_apply_boundary_BCs
    procedure :: conserved_to_primitive
    procedure :: primitive_to_conserved
    procedure :: compute_characteristic_speeds
    procedure :: compute_flux_F1
    procedure :: compute_pressure_flux
    procedure :: compute_sources
    procedure :: compute_lapse

    end type einstein_euler


    contains

    !####################################################
    !# Initialise System
    !####################################################

    subroutine einstein_euler_init(this)
        implicit none
        class(einstein_euler), intent(inout) :: this

        this%Nvar = 4 ! Number of variables
        this%CS_x = 1.0_dp
        
        this%K = aux1

        allocate(this%pressure(1-Ngz:Nx+Ngz))
        this%pressure = 1.0_dp

        allocate(this%lapse_tmp(1-Ngz:Nx+Ngz))
        this%lapse_tmp = 0.0_dp

    
    end subroutine einstein_euler_init


    !####################################################
    !# Solve equation for the Lapse
    !####################################################

    subroutine compute_lapse(this,cgrid,t,u,prims,lapse)
        implicit none
        class(einstein_euler), intent(inout) :: this
        type(grid), intent(in) :: cgrid
        real(dp), intent(in) :: t
        real(dp), intent(in) :: u(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,this%Nvar)
        real(dp), intent(in) :: prims(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp), intent(out) :: lapse(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz)
        real(dp) :: a(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz), mu(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz), v(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz), Gamma2(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz)
        real(dp) :: dx, f_i, f_i1
        integer :: i, Nx, Ngz

        dx = cgrid%dx 
        Nx = cgrid%Nx
        Ngz = cgrid%Ngz

        a = u(:,3)
        mu = prims(:,1)
        v = prims(:,2)
        Gamma2 = 1.0_dp/(1.0_dp - v**2)

        lapse(Nx) = 1.0_dp/u(Nx,3)
        ! print *, mu
        do i = Nx-1,1,-1 ! Need to include step parameter for decreasing loop index

            ! ODE RHS at index i
            if (i <= this%i_right) then
                f_i = (a(i)**2 - 1.0_dp)/(2.0_dp*cgrid%coordinates_x(i)) &
                    + (cgrid%coordinates_x(i)/2.0_dp)*a(i)**2*(-1.0_dp + this%K*mu(i) + v(i)**2*Gamma2(i)*((this%K + 1.0_dp)*mu(i) - 1.0_dp))

            else 
                f_i = (a(i)**2 - 1.0_dp)/(2.0_dp*cgrid%coordinates_x(i)) 
            end if
            
            ! ODE RHS at index i+1
            if (i+1 <= this%i_right) then
                f_i1 = (a(i+1)**2 - 1.0_dp)/(2.0_dp*cgrid%coordinates_x(i+1)) &
                     + (cgrid%coordinates_x(i+1)/2.0_dp)*a(i+1)**2*(-1.0_dp + this%K*mu(i+1) + v(i+1)**2*Gamma2(i+1)*((this%K + 1.0_dp)*mu(i+1) - 1.0_dp))
            else 
                f_i1 = (a(i+1)**2 - 1.0_dp)/(2.0_dp*cgrid%coordinates_x(i+1)) 
            end if

            ! Solve for lapse using the trapezoid rule for ODEs 
            ! The formula below is obtained by solving for \del_{r}log(alpha) = f(r) and then re-expressing in terms of alpha
            lapse(i) = lapse(i+1)*exp(-0.5_dp*dx*(f_i + f_i1))

        end do

        ! Fill in left ghost points using symmetries and right ghost points using first order extrapolation
        do i = 1, Ngz
            lapse(Nx+i) = lapse(Nx)
            lapse(1-i) = lapse(i)
        end do

    end subroutine compute_lapse




    !####################################################
    !# Primitive and Conserved Variable Conversions
    !####################################################


    subroutine conserved_to_primitive(this,cgrid,t,interior,cons,prims)
        implicit none
        class(einstein_euler), intent(inout) :: this
        type(grid), intent(in) :: cgrid
        real(dp), intent(in) :: t
        logical :: interior
        real(dp), intent(in) :: cons(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp), intent(out) :: prims(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp) :: v(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz), Gamma2(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz), &
                    mu(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz)
        real(dp) :: K
        integer :: left, right, Nx, Ngz, i


        K = this%K

        Nx = cgrid%Nx
        Ngz = cgrid%Ngz


        ! First if case if the primitives should only be recovered at interior points
        ! second case when recovering primitives at ghost points too
        ! Recall that in spherical symmetry, the left boundary is fixed at r=0!
        if (interior .eqv. .True.) then
            left = this%i_left
            right = this%i_right
        else 
            left = this%i_left - Ngz
            right = this%i_right + Ngz +1
        end if

        !--------- Use Explicit Formula ------------

        ! Compute primitives
        v = 0.0_dp
        Gamma2 = 0.0_dp
        mu = 0.0_dp

        associate(T00 => cons(left:right,1), T01 => cons(left:right,2))

        mu(left:right) = (1.0_dp + (this%K - 1.0_dp)*T00 + sqrt((-1.0_dp + T00 + this%K*T00)**2 - 4.0_dp*this%K*T01**2))/(2.0_dp*this%K)

        v(left:right) = (2.0_dp*T01)/(-1.0_dp + (this%K + 1.0_dp)*T00 + sqrt((-1.0_dp + (this%K+1.0_dp)*T00)**2 - 4.0_dp*this%K*T01**2))


        !### End Evolution if Negative Pressures detected
        if (any(mu(this%i_left:this%i_right)<=1.0_dp/this%K)) then
            print *, "WARNING: Negative pressure in the interior of star!"
            print *, "t = ", t
            print *, minval(mu(this%i_left:this%i_right)-1.0_dp/this%K)
            print *, mu(this%i_left:this%i_right)<1.0_dp/this%K
            print *, abs(cgrid%coordinates_x(this%i_right)-this%x_right)
            stop
        end if

        end associate



        ! Update primitives array
        prims(:,1) = mu
        prims(:,2) = v

        


    end subroutine conserved_to_primitive


    subroutine primitive_to_conserved(this,cgrid,t,prims,cons)
        implicit none
        class(einstein_euler), intent(inout) :: this
        type(grid), intent(in) :: cgrid
        real(dp), intent(in) :: t
        real(dp), intent(in) :: prims(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp), intent(out) :: cons(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp) :: Gamma2, v, mu
        integer :: i

        do i = 1-cgrid%Ngz, cgrid%Nx+cgrid%Ngz

            mu = prims(i,1)
            v = prims(i,2)
            Gamma2 = 1.0_dp/(1.0_dp - v**(2))

            cons(i,1) = Gamma2*((this%K + 1.0_dp)*mu - 1.0_dp) - (this%K*mu - 1.0_dp)

            cons(i,2) = Gamma2*((this%K + 1.0_dp)*mu - 1.0_dp)*v


        end do

    end subroutine primitive_to_conserved


    !####################################################
    !# Characteristic Speeds, Fluxes, and Sources
    !####################################################

    subroutine compute_characteristic_speeds(this, cgrid, t, grav_plus, grav_minus,prims_plus, prims_minus, char_speed)
        implicit none
        class(einstein_euler), intent(inout) :: this
        type(grid), intent(in) :: cgrid
        real(dp), intent(in) :: t
        real(dp), intent(in) :: grav_plus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp), intent(in) :: grav_minus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp), intent(in) :: prims_plus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp), intent(in) :: prims_minus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp), intent(out) :: char_speed(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz)
        integer :: i, Nx, Ngz
        real(dp) :: v_plus, v_minus, a_plus, a_minus, alpha_plus, alpha_minus
        real(dp) :: K
        real(dp) :: lam1_plus, lam1_minus, lam2_plus, lam2_minus

        Nx = cgrid%Nx
        Ngz = cgrid%Ngz

        K = this%K

        do i = 0, Nx

            v_plus = prims_plus(i,2)
            v_minus = prims_minus(i,2)

            a_plus = grav_plus(i,1)
            a_minus = grav_minus(i,1)

            alpha_plus = grav_plus(i,2)
            alpha_minus = grav_minus(i,2)

            lam1_plus = (alpha_plus/a_plus)*((-1.0_dp+K)*v_plus + sqrt(K)*(1.0_dp - v_plus**2))/(-1.0_dp + K*v_plus**2)
            lam1_minus = (alpha_minus/a_minus)*((-1.0_dp+K)*v_minus + sqrt(K)*(1.0_dp - v_minus**2))/(-1.0_dp + K*v_minus**2)

            lam2_plus = (alpha_plus/a_plus)*((-1.0_dp+K)*v_plus - sqrt(K)*(1.0_dp - v_plus**2))/(-1.0_dp + K*v_plus**2)
            lam2_minus = (alpha_minus/a_minus)*((-1.0_dp+K)*v_minus - sqrt(K)*(1.0_dp - v_minus**2))/(-1.0_dp + K*v_minus**2)

            char_speed(i) = max(abs(lam1_plus), abs(lam2_plus), abs(lam1_minus), abs(lam2_minus))

        end do


    end subroutine

    subroutine compute_flux_F1(this,cgrid,t,grav_plus,grav_minus,prims_plus,prims_minus,flux_F1_plus, flux_F1_minus)
        implicit none
        class(einstein_euler), intent(inout) :: this
        type(grid), intent(in) :: cgrid
        real(dp), intent(in) :: t
        real(dp), intent(in) :: grav_plus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp), intent(in) :: grav_minus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp), intent(in) :: prims_plus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp), intent(in) :: prims_minus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp), intent(out) :: flux_F1_plus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp), intent(out) :: flux_F1_minus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        integer :: i, Nx, Ngz
        real(dp) :: mu_plus, mu_minus, v_plus, v_minus, a_plus, a_minus, alpha_plus, alpha_minus
        real(dp) :: Gamma2_plus, Gamma2_minus

        Nx = cgrid%Nx
        Ngz = cgrid%Ngz

        do i = 0, Nx


            mu_plus = prims_plus(i,1)
            mu_minus = prims_minus(i,1)

            v_plus = prims_plus(i,2)
            v_minus = prims_minus(i,2)

            Gamma2_plus = 1.0_dp/(1.0_dp - v_plus**(2))
            Gamma2_minus = 1.0_dp/(1.0_dp - v_minus**(2))

            a_plus = grav_plus(i,1)
            a_minus = grav_minus(i,1)

            alpha_plus = grav_plus(i,2)
            alpha_minus = grav_minus(i,2)
            
            ! Construct flux vector
            flux_F1_plus(i,1) = (alpha_plus/a_plus)*(Gamma2_plus*(mu_plus*(this%K + 1.0_dp) - 1.0_dp)*v_plus)
            flux_F1_minus(i,1) = (alpha_minus/a_minus)*(Gamma2_minus*(mu_minus*(this%K + 1.0_dp) - 1.0_dp)*v_minus)

            flux_F1_plus(i,2) = (alpha_plus/a_plus)*(Gamma2_plus*(mu_plus*(this%K + 1.0_dp) - 1.0_dp)*v_plus**(2))
            flux_F1_minus(i,2) = (alpha_minus/a_minus)*(Gamma2_minus*(mu_minus*(this%K + 1.0_dp) - 1.0_dp)*v_minus**(2))


        end do

    end subroutine compute_flux_F1


    subroutine compute_pressure_flux(this,cgrid,t,grav_plus,grav_minus,prims_plus,prims_minus,flux_pressure)
        implicit none
        class(einstein_euler), intent(inout) :: this
        type(grid), intent(in) :: cgrid
        real(dp), intent(in) :: t
        real(dp), intent(in) :: grav_plus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp), intent(in) :: grav_minus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp), intent(in) :: prims_plus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp), intent(in) :: prims_minus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp), intent(out) :: flux_pressure(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz)
        real(dp) :: flux_pressure_minus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz), flux_pressure_plus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz)
        integer :: i, Nx, Ngz
        real(dp) :: mu_plus, mu_minus, v_plus, v_minus, a_plus, a_minus, alpha_plus, alpha_minus
        real(dp) :: Gamma2_plus, Gamma2_minus


        Nx = cgrid%Nx
        Ngz = cgrid%Ngz

        do i = 0, Nx


            mu_plus = prims_plus(i,1)
            mu_minus = prims_minus(i,1)

            a_plus = grav_plus(i,1)
            a_minus = grav_minus(i,1)

            alpha_plus = grav_plus(i,2)
            alpha_minus = grav_minus(i,2)
            
            ! Construct flux vector
            flux_pressure_plus(i) = (alpha_plus/a_plus)*(this%K*mu_plus-1.0_dp)
            flux_pressure_minus(i) = (alpha_minus/a_minus)*(this%K*mu_minus-1.0_dp)

            flux_pressure(i) = 0.5_dp*(flux_pressure_minus(i) + flux_pressure_plus(i))


        end do

    end subroutine compute_pressure_flux


    subroutine compute_sources(this,cgrid,t,grav_vars,prims,sources)
        implicit none
        class(einstein_euler), intent(inout) :: this
        type(grid), intent(in) :: cgrid
        real(dp), intent(in) :: t
        real(dp), intent(in) :: grav_vars(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp), intent(in) :: prims(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp), intent(out) :: sources(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,3)
        real(dp) :: mu, v, alpha, a, Gamma2
        integer :: i, Nx, Ngz

        Nx = cgrid%Nx
        Ngz = cgrid%Ngz

        do i = 1, Nx

            a = grav_vars(i,1)
            alpha = grav_vars(i,2)
            mu = prims(i,1)
            v = prims(i,2)
            Gamma2 = 1.0_dp/(1.0_dp - v**2)


            sources(i,1) = 0.0_dp

            sources(i,2) = ((a**2 - 1.0_dp)*alpha*((this%K - 1.0_dp)*mu - 1.0_dp))/(2.0_dp*a*cgrid%coordinates_x(i)) &
                            + a*alpha*cgrid%coordinates_x(i)*mu*(1.0_dp - this%K*mu)

            ! Evolution equation for a is identically zero outside of the star (vacuum solutions are static in spherical symmetry)
            if (cgrid%coordinates_x(i) < this%x_right ) then
                sources(i,3) = -(cgrid%coordinates_x(i)/2.0_dp)*a**2*alpha*Gamma2*v*((this%K+1.0_dp)*mu - 1.0_dp)
            else 
                sources(i,3) = 0.0_dp
            end if

        end do



    end subroutine compute_sources


    !####################################################
    !# Compute RHS of boundary ODEs
    !####################################################

    subroutine einstein_euler_compute_boundary_rhs(this,t,u,x_boundary,cgrid,dudt)
        implicit none
        class(einstein_euler), intent(inout) :: this
        type(grid), intent(in) :: cgrid
        real(dp), intent(in) :: u(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,this%Nvar)
        real(dp), intent(in) :: x_boundary(2)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: dudt(2)
        real(dp) :: prims(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp) :: v(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz), v_left, v_right
        real(dp) :: l0, l1, l2, l3
        real(dp) :: mu(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz)
        real(dp) :: alpha(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz), a(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz)
        real(dp) :: alpha_right, a_right, theta, theta_min
        integer :: i, Nx

        Nx = cgrid%Nx

        ! Set the boundary values and indices for this timestep
        this%x_left = x_boundary(1)
        this%x_right = x_boundary(2)

        this%i_left = -1
        this%i_right = -1

        do i = 1, Nx
            if ((this%i_left == -1) .and. (cgrid%coordinates_x(i) > this%x_left)) then
                this%i_left = i
            end if

            if ((this%i_right == -1) .and. (cgrid%coordinates_x(i) >= this%x_right)) then
                this%i_right = i-1
            end if

            if ((this%i_left /= -1) .and. (this%i_right /= -1)) then
                exit
            end if 
        end do

        if ((this%i_left == -1) .or. (this%i_right == -1)) then
                print *, "WARNING: Fluid boundary index not assigned!"
                print *, "x_L = ", x_boundary(1)
                print *, "x_R = ",  x_boundary(2)
        end if 

        ! Compute primitive fluid variables
        call conserved_to_primitive(this,cgrid,t,.True.,u,prims)
        call compute_lapse(this,cgrid,t,u,prims,alpha)
        this%lapse_tmp = alpha
        a = u(:,3)
        v = prims(:,2)
        mu = prims(:,1)

        ! Lagrange Interpolation to obtain boundary values
        l0 = -1.0_dp/(6.0_dp*cgrid%dx**3)*(this%x_right-cgrid%coordinates_x(this%i_right-2))*(this%x_right-cgrid%coordinates_x(this%i_right-1))*(this%x_right-cgrid%coordinates_x(this%i_right))
        l1 = 1.0_dp/(2.0_dp*cgrid%dx**3)*(this%x_right-cgrid%coordinates_x(this%i_right-3))*(this%x_right-cgrid%coordinates_x(this%i_right-1))*(this%x_right-cgrid%coordinates_x(this%i_right))
        l2 = -1.0_dp/(2.0_dp*cgrid%dx**3)*(this%x_right-cgrid%coordinates_x(this%i_right-3))*(this%x_right-cgrid%coordinates_x(this%i_right-2))*(this%x_right-cgrid%coordinates_x(this%i_right))
        l3 = 1.0_dp/(6.0_dp*cgrid%dx**3)*(this%x_right-cgrid%coordinates_x(this%i_right-3))*(this%x_right-cgrid%coordinates_x(this%i_right-2))*(this%x_right-cgrid%coordinates_x(this%i_right-1))

        v_right = v(this%i_right-3)*l0 + v(this%i_right-2)*l1  + v(this%i_right-1)*l2  + v(this%i_right)*l3 
        alpha_right = alpha(this%i_right-3)*l0 + alpha(this%i_right-2)*l1  + alpha(this%i_right-1)*l2  + alpha(this%i_right)*l3 
        a_right = a(this%i_right-3)*l0 + a(this%i_right-2)*l1  + a(this%i_right-1)*l2  + a(this%i_right)*l3 


        ! In spherical symmetry only the right boundary is evolved
        dudt(1) = 0.0_dp
        dudt(2) = (alpha_right/a_right)*v_right


    end subroutine einstein_euler_compute_boundary_rhs


    !####################################################
    !# Compute RHS of PDE
    !####################################################

    subroutine einstein_euler_compute_rhs(this,t,u,cgrid,dudt)
        implicit none
        class(einstein_euler), intent(inout) :: this
        type(grid), intent(in) :: cgrid
        real(dp), intent(in) :: u(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,this%Nvar)
        real(dp), intent(in) :: t
        real(dp), intent(out) :: dudt(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,this%Nvar)


        ! Auxiliary Variables
        integer :: Nvar, Nx, Ngz
        real(dp) :: dx 
        real(dp) :: char_speed(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz)
        real(dp) :: prims(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp) :: cons(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp) :: grav_vars(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp) :: alpha(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz)

        ! Linear Reconstructions
        real(dp) :: prims_minus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp) :: prims_plus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp) :: cons_minus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp) :: cons_plus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp) :: grav_minus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp) :: grav_plus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)


        ! Fluxes 
        real(dp) :: flux_F1(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp) :: flux_F1_plus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2), &
                    flux_F1_minus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2), &
                    flux_pressure(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz)

        ! Sources
        real(dp) :: sources(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,3)

        ! Grid Variables for Spherical Symmetry
        real(dp) :: r_bdy(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz)

        ! Lagrange Interpolation Variables
        real(dp) :: l0, l1, l2, l3, v_right, alpha_right, a_right


        Ngz = cgrid%Ngz
        Nx = cgrid%Nx
        dx = cgrid%dx
        Nvar = this%Nvar
        
        ! Initialise time derivative
        dudt(:,:) = 0.0_dp


       
        !#---------------------------------------
        !# Recover Primitives
        !#---------------------------------------

        cons(:,1) = u(:,1) ! Tau
        cons(:,2) = u(:,2) ! S

        call conserved_to_primitive(this,cgrid,t,.False.,cons,prims)

        !#---------------------------------------
        !# Reconstruct Primitives at Cell Edges
        !#---------------------------------------

        ! call linear_reconstruct(2,cgrid,prims,prims_plus,prims_minus)
        call WENO5_reconstruct(2,cgrid,prims,prims_plus,prims_minus)
        
        ! Obtain reconstructed conservative variables
        call primitive_to_conserved(this,cgrid,t,prims_plus,cons_plus)
        call primitive_to_conserved(this,cgrid,t,prims_minus,cons_minus)


        !----------------------------------------------
        !# Reconstruct gravitational variables
        !----------------------------------------------

        grav_vars(:,1) = u(:,3) ! a
        grav_vars(:,2) = this%lapse_tmp ! alpha

        ! call linear_reconstruct(2,cgrid,grav_vars,grav_plus,grav_minus)
        call WENO5_reconstruct(2,cgrid,grav_vars,grav_plus,grav_minus)



        !#---------------------------------------
        !# Compute Characteristic Speeds
        !#---------------------------------------

        call compute_characteristic_speeds(this,cgrid,t,grav_plus,grav_minus,prims_plus,prims_minus,char_speed)

        ! Lagrange Interpolation to obtain boundary values
        l0 = -1.0_dp/(6.0_dp*cgrid%dx**3)*(this%x_right-cgrid%coordinates_x(this%i_right-2))*(this%x_right-cgrid%coordinates_x(this%i_right-1))*(this%x_right-cgrid%coordinates_x(this%i_right))
        l1 = 1.0_dp/(2.0_dp*cgrid%dx**3)*(this%x_right-cgrid%coordinates_x(this%i_right-3))*(this%x_right-cgrid%coordinates_x(this%i_right-1))*(this%x_right-cgrid%coordinates_x(this%i_right))
        l2 = -1.0_dp/(2.0_dp*cgrid%dx**3)*(this%x_right-cgrid%coordinates_x(this%i_right-3))*(this%x_right-cgrid%coordinates_x(this%i_right-2))*(this%x_right-cgrid%coordinates_x(this%i_right))
        l3 = 1.0_dp/(6.0_dp*cgrid%dx**3)*(this%x_right-cgrid%coordinates_x(this%i_right-3))*(this%x_right-cgrid%coordinates_x(this%i_right-2))*(this%x_right-cgrid%coordinates_x(this%i_right-1))

        v_right = prims(this%i_right-3,2)*l0 + prims(this%i_right-2,2)*l1  + prims(this%i_right-1,2)*l2  + prims(this%i_right,2)*l3 
        a_right = grav_vars(this%i_right-3,1)*l0 + grav_vars(this%i_right-2,1)*l1  + grav_vars(this%i_right-1,1)*l2  + grav_vars(this%i_right,1)*l3 
        alpha_right = grav_vars(this%i_right-3,2)*l0 + grav_vars(this%i_right-2,2)*l1  + grav_vars(this%i_right-1,2)*l2  + grav_vars(this%i_right,2)*l3 

        this%CS_x = max(maxval(char_speed(0:Nx)),abs(alpha_right*v_right/a_right))

        !#---------------------------------------
        !# Construct the Fluxes
        !#---------------------------------------

        call compute_flux_F1(this,cgrid,t,grav_plus,grav_minus,prims_plus,prims_minus,flux_F1_plus,flux_F1_minus)
        call compute_pressure_flux(this,cgrid,t,grav_plus,grav_minus,prims_plus,prims_minus,flux_pressure)


        call local_lax_friedrichs(cgrid,2,flux_F1_plus,flux_F1_minus,cons_plus,cons_minus,char_speed,flux_F1)



        !#---------------------------------------
        !# Construct the Sources
        !#---------------------------------------

        call compute_sources(this,cgrid,t,grav_vars,prims,sources)


        !#---------------------------------------
        !# Conmpute time derivatives
        !#---------------------------------------

        r_bdy = cgrid%coordinates_x + 0.5_dp*cgrid%dx

        ! Tau evolution
        dudt(1:Nx,1) = -3.0_dp/(r_bdy(1:Nx)**3 - r_bdy(0:Nx-1)**3)*(r_bdy(1:Nx)**2*flux_F1(1:Nx,1) - r_bdy(0:Nx-1)**2*flux_F1(0:Nx-1,1)) &
                        + sources(1:Nx,1)

        ! S evolution
        dudt(1:Nx,2) = -3.0_dp/(r_bdy(1:Nx)**3 - r_bdy(0:Nx-1)**3)*(r_bdy(1:Nx)**2*flux_F1(1:Nx,2) - r_bdy(0:Nx-1)**2*flux_F1(0:Nx-1,2)) &
                        - 1.0_dp/(dx)*(flux_pressure(1:Nx) - flux_pressure(0:Nx-1)) &
                        + sources(1:Nx,2)


        ! a evolution
        dudt(1:Nx,3) = sources(1:Nx,3)

        dudt(1:Nx,4) = 0.0_dp


    end subroutine einstein_euler_compute_rhs

    !------------------------------------------------
    ! Computational Domain Boundary Conditions
    !------------------------------------------------

    subroutine einstein_euler_apply_BCs(this,cgrid,u)
        implicit none
        class(einstein_euler), intent(in) :: this
        type(grid), intent(in) :: cgrid
        real(dp), intent(inout) :: u(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,this%Nvar)
        integer :: Nx, Ngz, i

        Nx = cgrid%Nx
        Ngz = cgrid%Ngz


        ! Right Boundary (First Order Extrapolation)
        do i = 1,cgrid%Ngz
            u(Nx+i,:) = u(Nx,:)
        end do

        ! Left Boundary
        do i = 1, cgrid%Ngz

            u(1-i,1) = u(i,1)
            u(1-i,2) = -u(i,2)
            u(1-i,3) = u(i,3)

        end do 


    end subroutine einstein_euler_apply_BCs


    !------------------------------------------------
    ! Fluid Domain Free Boundary Conditions
    !------------------------------------------------

    subroutine einstein_euler_apply_boundary_BCs(this,cgrid,t,u)
        implicit none
        class(einstein_euler), intent(inout) :: this
        type(grid), intent(in) :: cgrid
        real(dp), intent(in) :: t
        real(dp), intent(inout) :: u(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,this%Nvar)
        real(dp) :: prims(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp) :: cons(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,2)
        real(dp) :: mu(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz), v(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz), a(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz)
        real(dp) :: alpha, alpha_tilde, v_left, v_right, mu_left_D1, mu_right_D1, v_left_D1, v_right_D1,  mu_0, Gamma2_xR, Gamma2_xL
        real(dp) ::  v_right_D2, mu_right_D2
        real(dp) :: a_right
        real(dp) :: l0, l1, l2, l3, l0_D1, l1_D1, l2_D1, l3_D1, l0_D2, l1_D2, l2_D2, l3_D2
        real(dp) :: ell1_1,ell1_2, ell2_1, ell2_2
        real(dp) :: gamma_plus(4), gamma_minus(4)
        real(dp) :: gamma_plus_right, gamma_plus_left, gamma_minus_right, gamma_minus_left
        real(dp) :: gamma_plus_right_D1, gamma_plus_left_D1, gamma_minus_right_D1, gamma_minus_left_D1
        integer :: Nx, Ngz, i, j

        Nx = cgrid%Nx
        Ngz = cgrid%Ngz

        ! Value of density where pressure vanishes
        mu_0 = 1.0_dp/(this%K)

        ! Compute primitive fluid variables
        call conserved_to_primitive(this,cgrid,t,.True.,u,prims)
        mu = prims(:,1)
        v = prims(:,2)
        a = u(:,3)


        !######################## Boundary Conditions using Characteristic Decomposition ###############################

        !### Right Boundary

        ! Lagrange polynomials for zeroth and first derivatives
        l0 = -1.0_dp/(6.0_dp*cgrid%dx**3)*(this%x_right-cgrid%coordinates_x(this%i_right-2))*(this%x_right-cgrid%coordinates_x(this%i_right-1))*(this%x_right-cgrid%coordinates_x(this%i_right))
        l1 = 1.0_dp/(2.0_dp*cgrid%dx**3)*(this%x_right-cgrid%coordinates_x(this%i_right-3))*(this%x_right-cgrid%coordinates_x(this%i_right-1))*(this%x_right-cgrid%coordinates_x(this%i_right))
        l2 = -1.0_dp/(2.0_dp*cgrid%dx**3)*(this%x_right-cgrid%coordinates_x(this%i_right-3))*(this%x_right-cgrid%coordinates_x(this%i_right-2))*(this%x_right-cgrid%coordinates_x(this%i_right))
        l3 = 1.0_dp/(6.0_dp*cgrid%dx**3)*(this%x_right-cgrid%coordinates_x(this%i_right-3))*(this%x_right-cgrid%coordinates_x(this%i_right-2))*(this%x_right-cgrid%coordinates_x(this%i_right-1))

        ! 1st Derivative
        l0_D1 = -1.0_dp/(6.0_dp*cgrid%dx**3)*( (this%x_right-cgrid%coordinates_x(this%i_right-1))*(this%x_right-cgrid%coordinates_x(this%i_right)) &
                + (this%x_right-cgrid%coordinates_x(this%i_right-2))*(this%x_right-cgrid%coordinates_x(this%i_right))  &
                + (this%x_right-cgrid%coordinates_x(this%i_right-2))*(this%x_right-cgrid%coordinates_x(this%i_right-1)) )

        l1_D1 = 1.0_dp/(2.0_dp*cgrid%dx**3)*( (this%x_right-cgrid%coordinates_x(this%i_right-1))*(this%x_right-cgrid%coordinates_x(this%i_right)) &
                + (this%x_right-cgrid%coordinates_x(this%i_right-3))*(this%x_right-cgrid%coordinates_x(this%i_right))  &
                + (this%x_right-cgrid%coordinates_x(this%i_right-3))*(this%x_right-cgrid%coordinates_x(this%i_right-1)) )

        l2_D1 = -1.0_dp/(2.0_dp*cgrid%dx**3)*( (this%x_right-cgrid%coordinates_x(this%i_right-2))*(this%x_right-cgrid%coordinates_x(this%i_right)) &
                + (this%x_right-cgrid%coordinates_x(this%i_right-3))*(this%x_right-cgrid%coordinates_x(this%i_right))  &
                + (this%x_right-cgrid%coordinates_x(this%i_right-3))*(this%x_right-cgrid%coordinates_x(this%i_right-2)) )

        l3_D1 = 1.0_dp/(6.0_dp*cgrid%dx**3)*( (this%x_right-cgrid%coordinates_x(this%i_right-2))*(this%x_right-cgrid%coordinates_x(this%i_right-1)) &
                + (this%x_right-cgrid%coordinates_x(this%i_right-3))*(this%x_right-cgrid%coordinates_x(this%i_right-1))  &
                + (this%x_right-cgrid%coordinates_x(this%i_right-3))*(this%x_right-cgrid%coordinates_x(this%i_right-2)) )



        ! Extrapolate velocity to boundary and compute left eigenvectors at x_{R}
        v_right = v(this%i_right-3)*l0 + v(this%i_right-2)*l1  + v(this%i_right-1)*l2  + v(this%i_right)*l3 
        a_right = a(this%i_right-3)*l0 + a(this%i_right-2)*l1  + a(this%i_right-1)*l2  + a(this%i_right)*l3 

        ! ! l^{1} left eigenvector
        ! ell1_1 = (this%K-v_right**2)/((this%K - 1.0_dp)*v_right - (1.0_dp - v_right**2)*sqrt(this%K))
        ! ell1_2 = 1.0_dp

        ! ! l^{2} left eigenector
        ! ell2_1 = (v_right**2 - this%K)/((this%K - 1.0_dp)*v_right + (1.0_dp - v_right**2)*sqrt(this%K))
        ! ell2_2 = -1.0_dp


        ! l^{1} left eigenvector
        ell1_1 = (this%K*(1 - v_right**2)/(sqrt(this%K)*((this%K+1.0_dp)*mu_0 - 1.0_dp)))
        ell1_2 = 1.0_dp

        ! l^{2} left eigenector
        ell2_1 = -(this%K*(1 - v_right**2)/(sqrt(this%K)*((this%K+1.0_dp)*mu_0 - 1.0_dp)))
        ell2_2 = 1.0_dp
        

        ! Compute characteristic variables at interpolation nodes
        gamma_plus(1:4) = ell1_1*mu(this%i_right-3:this%i_right) + ell1_2*v(this%i_right-3:this%i_right) 
        ! gamma_minus(1:4) = ell2_1*mu(this%i_right-3:this%i_right) + ell2_2*v(this%i_right-3:this%i_right) 

        ! Extrapolate outgoing variables to the surface:
        gamma_plus_right = gamma_plus(1)*l0 + gamma_plus(2)*l1 + gamma_plus(3)*l2 + gamma_plus(4)*l3 
        gamma_plus_right_D1 =  gamma_plus(1)*l0_D1 + gamma_plus(2)*l1_D1 + gamma_plus(3)*l2_D1 + gamma_plus(4)*l3_D1

        ! Use boundary conditions to obtain ingoing variables at the surface
        Gamma2_xR = 1.0_dp/(1.0_dp - v_right**2)
        gamma_minus_right = ell2_1*mu_0 + ell2_2*v_right

        ! Solve for the first derivatives using the extrapolated outgoing values
        mu_right_D1 = (gamma_plus_right_D1 -ell1_2*(-2.0_dp*v_right/this%x_right + 0.5_dp*a_right**2*this%x_right*mu_0*v_right))/(ell1_1 + ell1_2*this%K*v_right/(Gamma2_xR*mu_0))

        v_right_D1 = this%K*v_right/(Gamma2_xR*mu_0)*mu_right_D1 -2.0_dp*v_right/this%x_right + 0.5_dp*a_right**2*this%x_right*mu_0*v_right

        ! Now that we have the values at the surface, we compute the Taylor expansion to fill ghost cells.
        do i = 1, Ngz+1
            v(this%i_right+i)  = v_right + (cgrid%coordinates_x(this%i_right+i)-this%x_right)*v_right_D1 
            mu(this%i_right+i) = mu_0  + (cgrid%coordinates_x(this%i_right+i)-this%x_right)*mu_right_D1  
        end do


        !### Left Boundary (mu is even, v is odd)
        do i = 1, cgrid%Ngz

            mu(1-i) = mu(i)
            v(1-i) = -v(i)

        end do 



        ! Convert back to conserved variables
        prims(:,1) = mu
        prims(:,2) = v


        call primitive_to_conserved(this,cgrid,t,prims,u)



    end subroutine einstein_euler_apply_boundary_BCs



end module einstein_euler_mod