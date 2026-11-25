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
        real(dp) :: a(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz), mu(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz), v(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz), Gamma2(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz), &
                    pressure(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz)
        real(dp) :: dx, f_i, f_i1
        integer :: i, Nx, Ngz

        dx = cgrid%dx 
        Nx = cgrid%Nx
        Ngz = cgrid%Ngz

        a = u(:,3)
        mu = prims(:,1)
        v = prims(:,2)
        Gamma2 = 1.0_dp/(1.0_dp - v**2)

        pressure = max(this%K*mu - 1.0_dp, 0.0_dp)

        lapse(Nx) = 1.0_dp/u(Nx,3)
        ! print *, mu
        do i = Nx-1,1,-1 ! Need to include step parameter for decreasing loop index

            ! ODE RHS at index i
            f_i = (a(i)**2 - 1.0_dp)/(2.0_dp*cgrid%coordinates_x(i)) &
                + (cgrid%coordinates_x(i)/2.0_dp)*a(i)**2*(pressure(i) + v(i)**2*Gamma2(i)*(mu(i) + pressure(i)))

            
            ! ODE RHS at index i+1
            f_i1 = (a(i+1)**2 - 1.0_dp)/(2.0_dp*cgrid%coordinates_x(i+1)) &
                    + (cgrid%coordinates_x(i+1)/2.0_dp)*a(i+1)**2*(pressure(i+1) + v(i+1)**2*Gamma2(i+1)*(mu(i+1) + pressure(i+1)))


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
        real(dp) :: v(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz), Gamma2, &
                    mu(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz)
        real(dp) :: K
        integer :: i, Nx, Ngz


        K = this%K

        Nx = cgrid%Nx
        Ngz = cgrid%Ngz


        !--------- Use Explicit Formula ------------

        ! Compute primitives
        v = 0.0_dp
        Gamma2 = 0.0_dp
        mu = 0.0_dp


        do i = 1-Ngz, Nx+Ngz

            if (cons(i,1) >= 0.5_dp*(1.0_dp/K) + 0.5_dp*sqrt((1.0_dp + 4.0_dp*K**2*cons(i,2)**2)/(K**2))) then
        
                mu(i) = (1.0_dp + (this%K - 1.0_dp)*cons(i,1) + sqrt((-1.0_dp + cons(i,1) + this%K*cons(i,1))**2 - 4.0_dp*this%K*cons(i,2)**2))/(2.0_dp*this%K)

                v(i) = (2.0_dp*cons(i,2))/(-1.0_dp + (this%K + 1.0_dp)*cons(i,1)+ sqrt((-1.0_dp + (this%K+1.0_dp)*cons(i,1))**2 - 4.0_dp*this%K*cons(i,2)**2))

            else 
                v(i) = cons(i,2)/cons(i,1)
                Gamma2 = 1.0_dp/(1.0_dp - v(i)**2)
                mu(i) = cons(i,1)/Gamma2

            end if


        end do



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
        real(dp) :: Gamma2, v, mu, pressure
        integer :: i

        do i = 1-cgrid%Ngz, cgrid%Nx+cgrid%Ngz

            mu = prims(i,1)
            v = prims(i,2)
            Gamma2 = 1.0_dp/(1.0_dp - v**(2))
            pressure = max(this%K*mu -1.0_dp,0.0_dp)

            cons(i,1) = Gamma2*(mu+pressure) - (pressure)

            cons(i,2) = Gamma2*(mu + pressure)*v


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
        real(dp) :: v_plus, v_minus, a_plus, a_minus, alpha_plus, alpha_minus, mu_plus, mu_minus
        real(dp) :: K
        real(dp) :: lam1_plus, lam1_minus, lam2_plus, lam2_minus, pressure_minus, pressure_plus

        Nx = cgrid%Nx
        Ngz = cgrid%Ngz

        K = this%K

        do i = 0, Nx

            v_plus = prims_plus(i,2)
            v_minus = prims_minus(i,2)

            mu_plus = prims_plus(i,1)
            mu_minus = prims_minus(i,1)

            a_plus = grav_plus(i,1)
            a_minus = grav_minus(i,1)

            alpha_plus = grav_plus(i,2)
            alpha_minus = grav_minus(i,2)

            if (mu_plus >= 1.0_dp/this%K) then
                pressure_plus = this%K
            else 
                pressure_plus = 0.0_dp
            end if

            if (mu_minus >= 1.0_dp/this%K) then
                pressure_minus = this%K
            else 
                pressure_minus = 0.0_dp
            end if

            lam1_plus = (alpha_plus/a_plus)*((-1.0_dp+pressure_plus)*v_plus + sqrt(pressure_plus)*(1.0_dp - v_plus**2))/(-1.0_dp + pressure_plus*v_plus**2)
            lam1_minus = (alpha_minus/a_minus)*((-1.0_dp+pressure_minus)*v_minus + sqrt(pressure_minus)*(1.0_dp - v_minus**2))/(-1.0_dp + pressure_minus*v_minus**2)


            lam2_plus = (alpha_plus/a_plus)*((-1.0_dp+pressure_plus)*v_plus - sqrt(pressure_plus)*(1.0_dp - v_plus**2))/(-1.0_dp + pressure_plus*v_plus**2)
            lam2_minus = (alpha_minus/a_minus)*((-1.0_dp+pressure_minus)*v_minus - sqrt(pressure_minus)*(1.0_dp - v_minus**2))/(-1.0_dp + pressure_minus*v_minus**2)

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
        real(dp) :: Gamma2_plus, Gamma2_minus, pressure_minus, pressure_plus

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

            if (mu_plus >= 1.0_dp/this%K) then
                pressure_plus = this%K*mu_plus -1.0_dp
            else 
                pressure_plus = 0.0_dp
            end if

            if (mu_minus >= 1.0_dp/this%K) then
                pressure_minus = this%K*mu_minus -1.0_dp
            else 
                pressure_minus = 0.0_dp
            end if
            
            ! Construct flux vector
            flux_F1_plus(i,1) = (alpha_plus/a_plus)*(Gamma2_plus*(mu_plus + pressure_plus)*v_plus)
            flux_F1_minus(i,1) = (alpha_minus/a_minus)*(Gamma2_minus*(mu_minus + pressure_minus)*v_minus)

            flux_F1_plus(i,2) = (alpha_plus/a_plus)*(Gamma2_plus*(mu_plus + pressure_plus)*v_plus**(2))
            flux_F1_minus(i,2) = (alpha_minus/a_minus)*(Gamma2_minus*(mu_minus + pressure_minus)*v_minus**(2))


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
        real(dp) :: Gamma2_plus, Gamma2_minus, pressure_minus, pressure_plus


        Nx = cgrid%Nx
        Ngz = cgrid%Ngz

        do i = 0, Nx


            mu_plus = prims_plus(i,1)
            mu_minus = prims_minus(i,1)

            a_plus = grav_plus(i,1)
            a_minus = grav_minus(i,1)

            alpha_plus = grav_plus(i,2)
            alpha_minus = grav_minus(i,2)

            if (mu_plus >= 1.0_dp/this%K) then
                pressure_plus = this%K*mu_plus -1.0_dp
            else 
                pressure_plus = 0.0_dp
            end if

            if (mu_minus >= 1.0_dp/this%K) then
                pressure_minus = this%K*mu_minus -1.0_dp
            else 
                pressure_minus = 0.0_dp
            end if
            
            ! Construct flux vector
            flux_pressure_plus(i) = (alpha_plus/a_plus)*(pressure_plus)
            flux_pressure_minus(i) = (alpha_minus/a_minus)*(pressure_minus)

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
        real(dp) :: mu, v, alpha, a, Gamma2, pressure
        integer :: i, Nx, Ngz

        Nx = cgrid%Nx
        Ngz = cgrid%Ngz

        do i = 1, Nx

            a = grav_vars(i,1)
            alpha = grav_vars(i,2)
            mu = prims(i,1)
            v = prims(i,2)
            Gamma2 = 1.0_dp/(1.0_dp - v**2)

            pressure = max(this%K*mu -1.0_dp, 0.0_dp)


            sources(i,1) = 0.0_dp

            sources(i,2) = ((a**2 - 1.0_dp)*alpha*(pressure - mu ))/(2.0_dp*a*cgrid%coordinates_x(i)) &
                            - a*alpha*cgrid%coordinates_x(i)*mu*pressure

            sources(i,3) = -(cgrid%coordinates_x(i)/2.0_dp)*a**2*alpha*Gamma2*v*(pressure + mu)


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
        real(dp) :: mu(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz)
        real(dp) :: alpha(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz), a(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz)
        real(dp) :: alpha_right, a_right
        integer :: i, Nx

        Nx = cgrid%Nx

        ! Boundary not tracked for atmosphere algorithm
        dudt(1) = 0.0_dp
        dudt(2) = 0.0_dp


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
        call compute_lapse(this,cgrid,t,u,prims,alpha)
        this%lapse_tmp = alpha

        !#---------------------------------------
        !# Reconstruct Primitives at Cell Edges
        !#---------------------------------------

        call linear_reconstruct(2,cgrid,prims,prims_plus,prims_minus)
        ! call WENO5_reconstruct(2,cgrid,prims,prims_plus,prims_minus)
        
        ! Obtain reconstructed conservative variables
        call primitive_to_conserved(this,cgrid,t,prims_plus,cons_plus)
        call primitive_to_conserved(this,cgrid,t,prims_minus,cons_minus)


        !----------------------------------------------
        !# Reconstruct gravitational variables
        !----------------------------------------------

        grav_vars(:,1) = u(:,3) ! a
        grav_vars(:,2) = this%lapse_tmp ! alpha

        call linear_reconstruct(2,cgrid,grav_vars,grav_plus,grav_minus)
        ! call WENO5_reconstruct(2,cgrid,grav_vars,grav_plus,grav_minus)



        !#---------------------------------------
        !# Compute Characteristic Speeds
        !#---------------------------------------

        call compute_characteristic_speeds(this,cgrid,t,grav_plus,grav_minus,prims_plus,prims_minus,char_speed)


        this%CS_x = maxval(char_speed(0:Nx))

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


    end subroutine einstein_euler_apply_boundary_BCs



end module einstein_euler_mod