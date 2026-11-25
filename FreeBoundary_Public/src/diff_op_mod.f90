module diff_op_mod
    use real_type_mod, only: dp
    use grid_mod
    implicit none

    contains

    !############################################
    !# Finite Difference Stencils
    !############################################

    !----------------------------------------------------
    ! Fourth-Order Central Stencil for 1st Derivative 
    !----------------------------------------------------

    subroutine FD_C4(dx,Nx,Ngz,u,du)
        ! Inputs:            
        ! dx - The spatial step size, real(dp)
        ! Nx - The number of grid points, integer
        ! Ngz - The number of ghost points, integer
        ! u - The function to be differentiated, real(dp) dimension(1,1-Ngz:Nx+Ngz)
        ! Output:
        ! du - The derivative, real(dp) dimension(1,1-Ngz:Nx+Ngz)
        implicit none
        real(dp), intent(in) :: dx 
        integer, intent(in) :: Nx, Ngz 
        real(dp), intent(in) :: u(1-Ngz:Nx+Ngz,1)
        real(dp), intent(out) :: du(1-Ngz:Nx+Ngz,1)
        ! integer :: i

        ! do i = 1, Nx
        !     du(1,i) = (u(1,i-2) - 8.0_dp*u(1,i-1) + 8.0_dp*u(1,i+1) - u(1, i+2))/(12.0_dp*dx)
        ! end do

        ! Initialise derivative array
        du(:,:) = 0.0_dp

        ! Compute derivative
        du(1:Nx,:) = (u(-1:Nx-2,:) - 8.0_dp*u(0:Nx-1,:) + 8.0_dp*u(2:Nx+1,:) - u(3:Nx+2,:))/(12.0_dp*dx)

    
    end subroutine FD_C4

    ! This version uses one-sided differences at the boundaries
    subroutine FD_OC4(dx,Nx,Ngz,Nvar,u,du)
        ! Inputs:            
        ! dx - The spatial step size, real(dp)
        ! Nx - The number of grid points, integer
        ! Ngz - The number of ghost points, integer
        ! u - The function to be differentiated, real(dp) dimension(1,1-Ngz:Nx+Ngz)
        ! Output:
        ! du - The derivative, real(dp) dimension(1,1-Ngz:Nx+Ngz)
        implicit none
        real(dp), intent(in) :: dx 
        integer, intent(in) :: Nx, Ngz, Nvar 
        real(dp), intent(in) :: u(1-Ngz:Nx+Ngz,Nvar)
        real(dp), intent(out) :: du(1-Ngz:Nx+Ngz,Nvar)
        ! integer :: i

        ! do i = 1, Nx
        !     du(1,i) = (u(1,i-2) - 8.0_dp*u(1,i-1) + 8.0_dp*u(1,i+1) - u(1, i+2))/(12.0_dp*dx)
        ! end do

        ! Initialise derivative array
        du(:,:) = 0.0_dp

        ! Compute derivative
        du(1:Nx,:) = (u(-1:Nx-2,:) - 8.0_dp*u(0:Nx-1,:) + 8.0_dp*u(2:Nx+1,:) - u(3:Nx+2,:))/(12.0_dp*dx)

        ! One-sided stencil at boundaries
        du(1,:) = (-25.0_dp*u(1,:) + 48.0_dp*u(2,:) - 36.0_dp*u(3,:) + 16.0_dp*u(4,:) - 3.0_dp*u(5,:))/(12.0_dp*dx)

        du(2,:) = (-25.0_dp*u(2,:) + 48.0_dp*u(3,:) - 36.0_dp*u(4,:) + 16.0_dp*u(5,:) - 3.0_dp*u(6,:))/(12.0_dp*dx)

        du(3,:) = (-25.0_dp*u(3,:) + 48.0_dp*u(4,:) - 36.0_dp*u(5,:) + 16.0_dp*u(6,:) - 3.0_dp*u(7,:))/(12.0_dp*dx)

        du(Nx,:) = (25.0_dp*u(Nx,:) - 48.0_dp*u(Nx-1,:) + 36.0_dp*u(Nx-2,:) - 16.0_dp*u(Nx-3,:) + 3.0_dp*u(Nx-4,:))/(12.0_dp*dx)

        du(Nx-1,:) = (25.0_dp*u(Nx-1,:) - 48.0_dp*u(Nx-2,:) + 36.0_dp*u(Nx-3,:) - 16.0_dp*u(Nx-4,:) + 3.0_dp*u(Nx-5,:))/(12.0_dp*dx)

        du(Nx-2,:) = (25.0_dp*u(Nx-2,:) - 48.0_dp*u(Nx-3,:) + 36.0_dp*u(Nx-4,:) - 16.0_dp*u(Nx-5,:) + 3.0_dp*u(Nx-6,:))/(12.0_dp*dx)

    
    end subroutine FD_OC4


    !----------------------------------------------------
    ! Fourth-Order Central Stencil for 3rd Derivative 
    !----------------------------------------------------

    ! Warning: This function is only used for saving norm of third derivative right now!
    ! Needs modification to use in PDEs etc.!
    subroutine FD3_C4(dx,Nx,Ngz,Nvar,u,du)
        ! Inputs:            
        ! dx - The spatial step size, real(dp)
        ! Nx - The number of grid points, integer
        ! Ngz - The number of ghost points, integer
        ! u - The function to be differentiated, real(dp) dimension(1,1-Ngz:Nx+Ngz)
        ! Output:
        ! du - The derivative, real(dp) dimension(1,1-Ngz:Nx+Ngz)
        implicit none
        real(dp), intent(in) :: dx 
        integer, intent(in) :: Nx, Ngz, Nvar
        real(dp), intent(inout) :: u(1-Ngz:Nx+Ngz,Nvar)
        real(dp), intent(out) :: du(1-Ngz:Nx+Ngz,Nvar)

        ! Initialise derivative array
        du(:,:) = 0.0_dp

        ! Centred stencil at interior
        du(1:Nx,:) = ((1.0_dp)/(8.0_dp)*u(-2:Nx-3,:) - u(-1:Nx-2,:) + (13.0_dp)/(8.0_dp)*u(0:Nx-1,:) - (13.0_dp)/(8.0_dp)*u(2:Nx+1,:) &
                        + 1.0_dp*u(3:Nx+2,:) - (1.0_dp)/(8.0_dp)*u(4:Nx+3,:))/(dx**(3))


        ! One-sided stencil at boundaries
        du(1,:) = (-(49.0_dp/8.0_dp)*u(1,:) + 29.0_dp*u(2,:) - (461.0_dp/8.0_dp)*u(3,:) &
                    + 62.0_dp*u(4,:) - (307.0_dp/8.0_dp)*u(5,:) + 13.0_dp*u(6,:) - (15.0_dp/8.0_dp)*u(7,:))/(dx**(3))

        du(2,:) = (-(49.0_dp/8.0_dp)*u(2,:) + 29.0_dp*u(3,:) - (461.0_dp/8.0_dp)*u(4,:) &
                    + 62.0_dp*u(5,:) - (307.0_dp/8.0_dp)*u(6,:) + 13.0_dp*u(7,:) - (15.0_dp/8.0_dp)*u(8,:))/(dx**(3))

        du(3,:) = (-(49.0_dp/8.0_dp)*u(3,:) + 29.0_dp*u(4,:) - (461.0_dp/8.0_dp)*u(5,:) &
                    + 62.0_dp*u(6,:) - (307.0_dp/8.0_dp)*u(7,:) + 13.0_dp*u(8,:) - (15.0_dp/8.0_dp)*u(9,:))/(dx**(3))

        du(Nx,:) = ((49.0_dp/8.0_dp)*u(Nx,:) - 29.0_dp*u(Nx-1,:) + (461.0_dp/8.0_dp)*u(Nx-2,:) &
                    - 62.0_dp*u(Nx-3,:) + (307.0_dp/8.0_dp)*u(Nx-4,:) - 13.0_dp*u(Nx-5,:) + (15.0_dp/8.0_dp)*u(Nx-6,:))/(dx**(3))

        du(Nx-1,:) = ((49.0_dp/8.0_dp)*u(Nx-1,:) - 29.0_dp*u(Nx-2,:) + (461.0_dp/8.0_dp)*u(Nx-3,:) &
                    - 62.0_dp*u(Nx-4,:) + (307.0_dp/8.0_dp)*u(Nx-5,:) - 13.0_dp*u(Nx-6,:) + (15.0_dp/8.0_dp)*u(Nx-7,:))/(dx**(3))

        du(Nx-2,:) = ((49.0_dp/8.0_dp)*u(Nx-2,:) - 29.0_dp*u(Nx-3,:) + (461.0_dp/8.0_dp)*u(Nx-4,:) &
                    - 62.0_dp*u(Nx-5,:) + (307.0_dp/8.0_dp)*u(Nx-6,:) - 13.0_dp*u(Nx-7,:) + (15.0_dp/8.0_dp)*u(Nx-8,:))/(dx**(3))

    
    end subroutine FD3_C4


    !############################################
    !# Polynomial Reconstruction Operators
    !############################################

    !------------------------------------------------
    ! Piecewise Linear Reconstruction (2nd Order)
    !------------------------------------------------

    ! This code should be changed to use loops like the WENO5 reconstruction
    ! so it doesn't have to initialise the arrays ratio_x/y and phi
    ! This might make it faster...

    subroutine linear_reconstruct(Nvar,cgrid,u,u_plus_x,u_minus_x)
        implicit none
        type(grid), intent(in) :: cgrid
        integer, intent(in) :: Nvar
        real(dp), intent(in) :: u(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,Nvar)
        real(dp), intent(out) :: u_plus_x(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,Nvar), &
                                u_minus_x(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,Nvar)
        real(dp) :: ratio_x(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,Nvar), &
                     phi(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,Nvar)
        integer :: Nx, Ngz 


        ! Initialise arrays to zero
        ratio_x = 0.0_dp
        phi = 0.0_dp
        u_plus_x = 0.0_dp
        u_minus_x = 0.0_dp

        Nx = cgrid%Nx
        Ngz = cgrid%Ngz

        ! Compute ratio of slopes for limiter along x-direction
        ! NB: We add small number to the denominator to 
        ! stop NaN issues when neighbouring points are close to equal
        ratio_x(0:Nx+1,:) = (u(0:Nx+1,:) - u(-1:Nx,:))/(u(1:Nx+2,:) - u(0:Nx+1,:) + 1e-16_dp) 

        ! Minmod limiter
        phi = max(0.0_dp,min(1.0_dp,ratio_x))

        ! Superbee limiter
        ! phi = max(0.0_dp,min(2*ratio_x,1.0_dp),min(ratio_x,2.0_dp))

        ! van Leer limiter
        ! phi = (ratio_x + abs(ratio_x))/(1+abs(ratio_x))

        ! Compute the slope-limited reconstruction
        u_minus_x(-1:Nx+1,:) = u(-1:Nx+1,:) + 0.5_dp*phi(-1:Nx+1,:)*(u(0:Nx+2,:) - u(-1:Nx+1,:))

        u_plus_x(-1:Nx,:) = u(0:Nx+1,:) - 0.5_dp*phi(0:Nx+1,:)*(u(1:Nx+2,:) - u(0:Nx+1,:))

        ! u_minus_x(1-Ngz:Nx+Ngz,:) = u(1-Ngz:Nx+Ngz,:)

        ! u_plus_x(1-Ngz:Nx+Ngz-1,:) = u(1-Ngz+1:Nx+Ngz,:)


    end subroutine linear_reconstruct


    !---------------------------------------------
    ! WENO5 Reconstruction (5th Order) 
    !---------------------------------------------
    
    ! Note that fifth-order reconstruction 
    ! requires 3 ghost points on each side!

    ! See Zhang, Y-T., and C-W. Shu. "ENO and WENO schemes",
    ! Handbook of numerical analysis, vol. 17, pp. 103-122. Elsevier, 2016.
    ! for details on stencil formulas


    subroutine WENO5_reconstruct(Nvar,cgrid,u,u_plus_x,u_minus_x)
        implicit none
        type(grid), intent(in) :: cgrid
        integer, intent(in) :: Nvar
        real(dp), intent(in) :: u(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,Nvar)
        real(dp), intent(out) :: u_plus_x(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,Nvar), &
                                u_minus_x(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,Nvar)

        ! Check number of ghost points
        if (cgrid%Ngz < 3) then
            print *, "ERROR: Three or more ghost points required for fifth-order reconstruction!"
            print *, "Terminating program..."
            stop
        end if

        call WENO5_left(Nvar,cgrid,u,u_minus_x)
        call WENO5_right(Nvar,cgrid,u,u_plus_x)

    end subroutine WENO5_reconstruct

    subroutine WENO5_left(Nvar,cgrid,u,u_minus_x)
        implicit none
        type(grid), intent(in) :: cgrid
        integer, intent(in) :: Nvar
        real(dp), intent(in) :: u(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,Nvar)
        real(dp), intent(out) :: u_minus_x(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,Nvar)

        ! Grid vars
        integer :: Nx, Ngz, i,  k

        ! Reconstruction Stencils                        
        real(dp) :: u0, u1, u2

        ! Weights                        
        real(dp) :: w0, w1, w2

        ! Smoothness Indicators and related parameters   
        real(dp) :: alpha0, alpha1, alpha2, alphasum

        real(dp) :: beta0, beta1, beta2
        real(dp), parameter :: d0 = 3.0_dp/10.0_dp, d1 = 3.0_dp/5.0_dp, d2 = 1.0_dp/10.0_dp

        ! Small parameter to avoid division by zero errors
        real(dp), parameter :: eps = 1e-6_dp

        ! Initialise arrays to zero
        u0 = 0.0_dp
        u1 = 0.0_dp
        u2 = 0.0_dp
        w0 = 0.0_dp
        w1 = 0.0_dp
        w2 = 0.0_dp
        alpha0 = 0.0_dp
        alpha1 = 0.0_dp
        alpha2 = 0.0_dp
        alphasum = 0.0_dp
        beta0 = 0.0_dp
        beta1 = 0.0_dp
        beta2 = 0.0_dp
        u_minus_x = 0.0_dp

        Nx = cgrid%Nx
        Ngz = cgrid%Ngz

        !######## Reconstruction along x-direction ###########

        ! Calculate Reconstruction Stencils

        
        do k = 1, Nvar
                do i = 0,Nx+1

                    u0 = 1.0_dp/3.0_dp*u(i,k) + 5.0_dp/6.0_dp*u(i+1,k) - 1.0_dp/6.0_dp*u(i+2,k)

                    u1 = -1.0_dp/6.0_dp*u(i-1,k) + 5.0_dp/6.0_dp*u(i,k) + 1.0_dp/3.0_dp*u(i+1,k)

                    u2 = 1.0_dp/3.0_dp*u(i-2,k) - 7.0_dp/6.0_dp*u(i-1,k) + 11.0_dp/6.0_dp*u(i,k)

                    ! Calculate smoothness indicators and weights

                    beta0 = 1.0_dp/4.0_dp*(3.0_dp*u(i,k) - 4.0_dp*u(i+1,k) + u(i+2,k))**(2) &
                                        + 13.0_dp/12.0_dp*(u(i,k) - 2.0_dp*u(i+1,k) + u(i+2,k))**(2) 

                    beta1 = 1.0_dp/4.0_dp*(u(i-1,k) - u(i+1,k))**(2) &
                                        + 13.0_dp/12.0_dp*(u(i-1,k) - 2.0_dp*u(i,k) + u(i+1,k))**(2) 


                    beta2 = 1.0_dp/4.0_dp*(u(i-2,k) -4.0_dp*u(i-1,k) + 3.0_dp*u(i,k))**(2) &
                                        + 13.0_dp/12.0_dp*(u(i-2,k) - 2.0_dp*u(i-1,k) + u(i,k))**(2) 

                                        
                    alpha0 = d0/(eps + beta0)**(2)
                    alpha1 = d1/(eps + beta1)**(2)
                    alpha2 = d2/(eps + beta2)**(2)

                    alphasum = alpha0 + alpha1 + alpha2

                    w0 = alpha0/alphasum
                    w1 = alpha1/alphasum
                    w2 = alpha2/alphasum

                    ! Calculate Reconstruction at interfaces

                    u_minus_x(i,k) = w0*u0 + w1*u1 + w2*u2

                end do
            end do

    
    end subroutine WENO5_left

    subroutine WENO5_right(Nvar,cgrid,u,u_plus_x)
        implicit none
        type(grid), intent(in) :: cgrid
        integer, intent(in) :: Nvar
        real(dp), intent(in) :: u(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,Nvar)
        real(dp), intent(out) :: u_plus_x(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,Nvar)

        ! Grid vars
        integer :: Nx, Ngz, i,  k

        ! Reconstruction Stencils                        
        real(dp) :: u0, u1, u2

        ! Weights                        
        real(dp) :: w0, w1, w2

        ! Smoothness Indicators and related parameters   
        real(dp) :: alpha0, alpha1, alpha2, alphasum

        real(dp) :: beta0, beta1, beta2
        real(dp), parameter :: d0 = 3.0_dp/10.0_dp, d1 = 3.0_dp/5.0_dp, d2 = 1.0_dp/10.0_dp

        ! Small parameter to avoid division by zero errors
        real(dp), parameter :: eps = 1e-6_dp

        ! Initialise arrays to zero
        u0 = 0.0_dp
        u1 = 0.0_dp
        u2 = 0.0_dp
        w0 = 0.0_dp
        w1 = 0.0_dp
        w2 = 0.0_dp
        alpha0 = 0.0_dp
        alpha1 = 0.0_dp
        alpha2 = 0.0_dp
        alphasum = 0.0_dp
        beta0 = 0.0_dp
        beta1 = 0.0_dp
        beta2 = 0.0_dp
        u_plus_x = 0.0_dp

        Nx = cgrid%Nx
        Ngz = cgrid%Ngz

        !######## Reconstruction along x-direction ###########

        ! Calculate Reconstruction Stencils

        
        do k = 1, Nvar
                do i = 0,Nx

                    u0 = 1.0_dp/3.0_dp*u(i+1,k) + 5.0_dp/6.0_dp*u(i,k) - 1.0_dp/6.0_dp*u(i-1,k)

                    u1 = -1.0_dp/6.0_dp*u(i+2,k) + 5.0_dp/6.0_dp*u(i+1,k) + 1.0_dp/3.0_dp*u(i,k)

                    u2 = 1.0_dp/3.0_dp*u(i+3,k) - 7.0_dp/6.0_dp*u(i+2,k) + 11.0_dp/6.0_dp*u(i+1,k)

                    ! Calculate smoothness indicators and weights

                    beta0 = 1.0_dp/4.0_dp*(3.0_dp*u(i+1,k) - 4.0_dp*u(i,k) + u(i-1,k))**(2) &
                                        + 13.0_dp/12.0_dp*(u(i+1,k) - 2.0_dp*u(i,k) + u(i-1,k))**(2) 

                    beta1 = 1.0_dp/4.0_dp*(u(i,k) - u(i+2,k))**(2) &
                                        + 13.0_dp/12.0_dp*(u(i+2,k) - 2.0_dp*u(i+1,k) + u(i,k))**(2) 


                    beta2 = 1.0_dp/4.0_dp*(u(i+3,k) -4.0_dp*u(i+2,k) + 3.0_dp*u(i+1,k))**(2) &
                                        + 13.0_dp/12.0_dp*(u(i+3,k) - 2.0_dp*u(i+2,k) + u(i+1,k))**(2) 

                                        
                    alpha0 = d0/(eps + beta0)**(2)
                    alpha1 = d1/(eps + beta1)**(2)
                    alpha2 = d2/(eps + beta2)**(2)

                    alphasum = alpha0 + alpha1 + alpha2

                    w0 = alpha0/alphasum
                    w1 = alpha1/alphasum
                    w2 = alpha2/alphasum

                    ! Calculate Reconstruction at interfaces

                    u_plus_x(i,k) = w0*u0 + w1*u1 + w2*u2

                end do
            end do


    end subroutine WENO5_right


    !###########################################################################
    !# Flux Computation
    !###########################################################################


    !---------------------------------------------
    ! Central Upwind/Local Lax-Friedrichs
    !---------------------------------------------

    subroutine local_lax_friedrichs(cgrid,Nvar_local,flux_plus,flux_minus,cons_plus,cons_minus,char_speed,flux)
        implicit none
        type(grid) ::  cgrid
        integer, intent(in) :: Nvar_local
        real(dp), intent(in) :: flux_plus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,Nvar_local), flux_minus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,Nvar_local)
        real(dp), intent(in) :: cons_plus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,Nvar_local), cons_minus(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,Nvar_local)
        real(dp), intent(in) :: char_speed(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz)
        real(dp), intent(out) :: flux(1-cgrid%Ngz:cgrid%Nx+cgrid%Ngz,Nvar_local)
        integer :: i, Nx

        flux = 0.0_dp
        Nx = cgrid%Nx

        do i = 1, Nvar_local
            flux(0:Nx,i) = 0.5_dp*( (flux_plus(0:Nx,i) + flux_minus(0:Nx,i))  - char_speed(0:Nx)*(cons_plus(0:Nx,i)  - cons_minus(0:Nx,i))) ! Future flux

            ! flux(0:Nx,i) = 0.5_dp*( (flux_plus(0:Nx,i) + flux_minus(0:Nx,i))  + char_speed(0:Nx)*(cons_plus(0:Nx,i)   - cons_minus(0:Nx,i))) ! Past flux
        end do

    end subroutine local_lax_friedrichs



end module diff_op_mod