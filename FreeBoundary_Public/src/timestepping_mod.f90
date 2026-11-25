module timestepping_mod
    use real_type_mod, only: dp
    use patch_mod
    use read_parameters_mod, only: Ngz, Nx, aux1
    implicit none

    real(dp), allocatable :: y0(:,:), y1(:,:), y2(:,:), y3(:,:)
    real(dp), allocatable :: k1(:,:), k2(:,:), k3(:,:), k4(:,:)
    real(dp), allocatable :: y0_bdy(:), y1_bdy(:), y2_bdy(:), y3_bdy(:)
    real(dp), allocatable :: k1_bdy(:), k2_bdy(:), k3_bdy(:), k4_bdy(:)


    contains


    !-------------------------------
    ! Initialise RK Arrays
    !--------------------------------
    
    ! This routine means we allocate the temporary arrays once
    ! at the start of the evolution loop. 
    subroutine initialise_rk_arrays(patches)
        implicit none
        type(patch) :: patches

        ! Substep arrays
        allocate(y0(1-Ngz:Nx+Ngz,patches%system%Nvar))
        allocate(y1(1-Ngz:Nx+Ngz,patches%system%Nvar))
        allocate(y2(1-Ngz:Nx+Ngz,patches%system%Nvar))
        allocate(y3(1-Ngz:Nx+Ngz,patches%system%Nvar))

        ! RHS arrays
        allocate(k1(1-Ngz:Nx+Ngz,patches%system%Nvar))
        allocate(k2(1-Ngz:Nx+Ngz,patches%system%Nvar))
        allocate(k3(1-Ngz:Nx+Ngz,patches%system%Nvar))
        allocate(k4(1-Ngz:Nx+Ngz,patches%system%Nvar))

        ! Substep boundary arrays
        allocate(y0_bdy(2))
        allocate(y1_bdy(2))
        allocate(y2_bdy(2))
        allocate(y3_bdy(2))

        ! RHS bounday arrays
        allocate(k1_bdy(2))
        allocate(k2_bdy(2))
        allocate(k3_bdy(2))
        allocate(k4_bdy(2))


    end subroutine initialise_rk_arrays



    !#######################################
    !# "Classic" RK4 Solver
    !#######################################
    subroutine rk4(patch_tmp,dt)
        implicit none
        type(patch), intent(inout) :: patch_tmp
        real(dp), intent(in) :: dt
        real(dp) :: t0

        t0 = patch_tmp%tl_current%t
        y0(:,:) = patch_tmp%tl_current%soln
        y0_bdy = patch_tmp%x_boundary

        call patch_tmp%system%apply_BCs(patch_tmp%cgrid,y0) ! Update External Boundary Conditions

        ! First, we evolve the boundaries of the star
        call patch_tmp%system%compute_boundary_rhs(t0,y0,y0_bdy,patch_tmp%cgrid,k1_bdy)
        y1_bdy = y0_bdy + 0.5_dp*dt*k1_bdy
        ! Update the star ghost points
        call patch_tmp%system%apply_boundary_BCs(patch_tmp%cgrid,t0,y0)
        call patch_tmp%system%compute_rhs(t0,y0(:,:),patch_tmp%cgrid,k1(:,:))
        y1(:,:) = y0(:,:) + 0.5_dp*dt*k1(:,:)

        call patch_tmp%system%apply_BCs(patch_tmp%cgrid,y1) ! Update External Boundary Conditions
        
        !Enforce Floor on Conserved Variables
        where ((y1(:,1) < 0.5_dp*(1.0_dp/aux1) + 0.5_dp*sqrt((1.0_dp + 4.0_dp*aux1**2*y1(:,2)**2)/(aux1**2))) .and. (patch_tmp%cgrid%coordinates_x < y1_bdy(2)) .and. (patch_tmp%cgrid%coordinates_x > 0.0_dp))
                y1(:,1) = 0.5_dp*(1.0_dp/aux1) + 0.5_dp*sqrt((1.0_dp + 4.0_dp*aux1**2*y1(:,2)**2)/(aux1**2)) + 1e-13_dp
        end where

        ! Zero all values outside of star
        where (patch_tmp%cgrid%coordinates_x>y1_bdy(2) )
                y1(:,1) = 0.0_dp
                y1(:,2) = 0.0_dp
        end where
        
        ! Step 2
        call patch_tmp%system%compute_boundary_rhs(t0+0.5_dp*dt,y1,y1_bdy,patch_tmp%cgrid,k2_bdy)
        y2_bdy = y0_bdy + 0.5_dp*dt*k2_bdy
        call patch_tmp%system%apply_boundary_BCs(patch_tmp%cgrid,t0+0.5_dp*dt,y1)
        call patch_tmp%system%compute_rhs(t0+0.5_dp*dt,y1(:,:),patch_tmp%cgrid,k2(:,:))
        y2(:,:) = y0(:,:) + 0.5_dp*dt*k2(:,:)

        call patch_tmp%system%apply_BCs(patch_tmp%cgrid,y2) ! Update External Boundary Conditions

        !Enforce Floor on Conserved Variables
        where (y2(:,1) < 0.5_dp*(1.0_dp/aux1) + 0.5_dp*sqrt((1.0_dp + 4.0_dp*aux1**2*y2(:,2)**2)/(aux1**2)) .and. (patch_tmp%cgrid%coordinates_x < y2_bdy(2)) .and. (patch_tmp%cgrid%coordinates_x > 0.0_dp))
                y2(:,1) = 0.5_dp*(1.0_dp/aux1) + 0.5_dp*sqrt((1.0_dp + 4.0_dp*aux1**2*y2(:,2)**2)/(aux1**2)) + 1e-13_dp
        end where

        ! Zero all values outside of star
        where (patch_tmp%cgrid%coordinates_x>y2_bdy(2) )
                y2(:,1) = 0.0_dp
                y2(:,2) = 0.0_dp
        end where

        ! Step 3
        call patch_tmp%system%compute_boundary_rhs(t0+0.5_dp*dt,y2,y2_bdy,patch_tmp%cgrid,k3_bdy)
        y3_bdy = y0_bdy + dt*k3_bdy
        call patch_tmp%system%apply_boundary_BCs(patch_tmp%cgrid,t0+0.5_dp*dt,y2)
        call patch_tmp%system%compute_rhs(t0+0.5_dp*dt,y2(:,:),patch_tmp%cgrid,k3(:,:))
        y3(:,:) = y0(:,:) + dt*k3(:,:)
 
        call patch_tmp%system%apply_BCs(patch_tmp%cgrid,y3) ! Update External Boundary Conditions

        !Enforce Floor on Conserved Variables
        where (y3(:,1) < 0.5_dp*(1.0_dp/aux1) + 0.5_dp*sqrt((1.0_dp + 4.0_dp*aux1**2*y3(:,2)**2)/(aux1**2)) .and. (patch_tmp%cgrid%coordinates_x < y3_bdy(2)) .and. (patch_tmp%cgrid%coordinates_x > 0.0_dp))
                y3(:,1) = 0.5_dp*(1.0_dp/aux1) + 0.5_dp*sqrt((1.0_dp + 4.0_dp*aux1**2*y3(:,2)**2)/(aux1**2)) + 1e-13_dp
        end where

        ! Zero all values outside of star
        where (patch_tmp%cgrid%coordinates_x>y3_bdy(2) )
                y3(:,1) = 0.0_dp
                y3(:,2) = 0.0_dp
        end where


        ! Step 4
        call patch_tmp%system%compute_boundary_rhs(t0+dt,y3,y3_bdy,patch_tmp%cgrid,k4_bdy)
        y0_bdy = y0_bdy +dt*((1.0_dp/6.0_dp)*k1_bdy(:) + (1.0_dp/3.0_dp)*k2_bdy(:) &
                        + (1.0_dp/3.0_dp)*k3_bdy(:)  + (1.0_dp/6.0_dp)*k4_bdy(:))
        call patch_tmp%system%apply_boundary_BCs(patch_tmp%cgrid,t0+dt,y3)
        call patch_tmp%system%compute_rhs(t0+dt,y3(:,:),patch_tmp%cgrid,k4(:,:))

        y0(:,:) = y0(:,:) + dt*((1.0_dp/6.0_dp)*k1(:,:) + (1.0_dp/3.0_dp)*k2(:,:) &
                        + (1.0_dp/3.0_dp)*k3(:,:)  + (1.0_dp/6.0_dp)*k4(:,:))

        call patch_tmp%system%apply_BCs(patch_tmp%cgrid,y0) ! Update External Boundary Conditions

        !Enforce Floor on Conserved Variables
        where (y0(:,1) < 0.5_dp*(1.0_dp/aux1) + 0.5_dp*sqrt((1.0_dp + 4.0_dp*aux1**2*y0(:,2)**2)/(aux1**2)) .and. (patch_tmp%cgrid%coordinates_x < y0_bdy(2)) .and. (patch_tmp%cgrid%coordinates_x > 0.0_dp))
                y0(:,1) = 0.5_dp*(1.0_dp/aux1) + 0.5_dp*sqrt((1.0_dp + 4.0_dp*aux1**2*y0(:,2)**2)/(aux1**2)) + 1e-13_dp
        end where

        ! Zero all values outside of star
        where (patch_tmp%cgrid%coordinates_x>y0_bdy(2) )
                y0(:,1) = 0.0_dp
                y0(:,2) = 0.0_dp
        end where
        
        ! Update current timelevel array
        patch_tmp%tl_current%soln = y0(:,:)
        patch_tmp%x_boundary = y0_bdy

        

    end subroutine rk4



    !#######################################
    !# SSP RK3 
    !#######################################
    subroutine rk3(patch_tmp,dt)
        implicit none
        type(patch), intent(inout) :: patch_tmp
        real(dp), intent(in) :: dt
        real(dp) :: t0

        t0 = patch_tmp%tl_current%t
        y0(:,:) = patch_tmp%tl_current%soln
        y0_bdy = patch_tmp%x_boundary

        call patch_tmp%system%apply_BCs(patch_tmp%cgrid,y0) ! Update External Boundary Conditions

        ! First, we evolve the boundaries of the star
        call patch_tmp%system%compute_boundary_rhs(t0,y0,y0_bdy,patch_tmp%cgrid,k1_bdy)
        y1_bdy = y0_bdy + dt*k1_bdy
        ! Update the star ghost points
        call patch_tmp%system%apply_boundary_BCs(patch_tmp%cgrid,t0,y0)
        call patch_tmp%system%compute_rhs(t0,y0(:,:),patch_tmp%cgrid,k1(:,:))
        y1(:,:) = y0(:,:) + dt*k1(:,:)

        call patch_tmp%system%apply_BCs(patch_tmp%cgrid,y1) ! Update External Boundary Conditions

        !Enforce Floor on Conserved Variables
        ! where ((y1(:,1) < 0.5_dp*(1.0_dp/aux1) + 0.5_dp*sqrt((1.0_dp + 4.0_dp*aux1**2*y1(:,2)**2)/(aux1**2))) .and. (patch_tmp%cgrid%coordinates_x < y1_bdy(2)) .and. (patch_tmp%cgrid%coordinates_x > 0.0_dp))
        !         y1(:,1) = 0.5_dp*(1.0_dp/aux1) + 0.5_dp*sqrt((1.0_dp + 4.0_dp*aux1**2*y1(:,2)**2)/(aux1**2)) + 1e-13_dp
        ! end where

        ! Zero all values outside of star
        where (patch_tmp%cgrid%coordinates_x>y1_bdy(2))
                y1(:,1) = 0.0_dp
                y1(:,2) = 0.0_dp
        end where
        
        ! Step 2
        call patch_tmp%system%compute_boundary_rhs(t0+dt,y1,y1_bdy,patch_tmp%cgrid,k2_bdy)
        y2_bdy = y0_bdy + 0.25_dp*dt*(k2_bdy + k1_bdy)
        call patch_tmp%system%apply_boundary_BCs(patch_tmp%cgrid,t0+dt,y1)
        call patch_tmp%system%compute_rhs(t0+dt,y1(:,:),patch_tmp%cgrid,k2(:,:))
        y2(:,:) = y0(:,:) + 0.25_dp*dt*(k2(:,:) + k1(:,:))

        call patch_tmp%system%apply_BCs(patch_tmp%cgrid,y2) ! Update External Boundary Conditions

        !Enforce Floor on Conserved Variables
        ! where (y2(:,1) < 0.5_dp*(1.0_dp/aux1) + 0.5_dp*sqrt((1.0_dp + 4.0_dp*aux1**2*y2(:,2)**2)/(aux1**2)) .and. (patch_tmp%cgrid%coordinates_x < y2_bdy(2)) .and. (patch_tmp%cgrid%coordinates_x > 0.0_dp))
        !         y2(:,1) = 0.5_dp*(1.0_dp/aux1) + 0.5_dp*sqrt((1.0_dp + 4.0_dp*aux1**2*y2(:,2)**2)/(aux1**2)) + 1e-13_dp
        ! end where

        ! Zero all values outside of star
        where (patch_tmp%cgrid%coordinates_x>y2_bdy(2) )
                y2(:,1) = 0.0_dp
                y2(:,2) = 0.0_dp
        end where

        ! Step 3
        call patch_tmp%system%compute_boundary_rhs(t0+0.5_dp*dt,y2,y2_bdy,patch_tmp%cgrid,k3_bdy)
        y0_bdy = y0_bdy + dt*(1.0_dp/6.0_dp*k1_bdy + 1.0_dp/6.0_dp*k2_bdy + 2.0_dp/3.0_dp*k3_bdy)
        call patch_tmp%system%apply_boundary_BCs(patch_tmp%cgrid,t0+0.5_dp*dt,y2)
        call patch_tmp%system%compute_rhs(t0+0.5_dp*dt,y2(:,:),patch_tmp%cgrid,k3(:,:))
        y0(:,:) = y0(:,:) + dt*(1.0_dp/6.0_dp*k1(:,:) + 1.0_dp/6.0_dp*k2(:,:) + 2.0_dp/3.0_dp*k3(:,:))
 
        call patch_tmp%system%apply_BCs(patch_tmp%cgrid,y0) ! Update External Boundary Conditions

        !Enforce Floor on Conserved Variables
        ! where (y0(:,1) < 0.5_dp*(1.0_dp/aux1) + 0.5_dp*sqrt((1.0_dp + 4.0_dp*aux1**2*y0(:,2)**2)/(aux1**2)) .and. (patch_tmp%cgrid%coordinates_x < y0_bdy(2)) .and. (patch_tmp%cgrid%coordinates_x > 0.0_dp))
        !         y0(:,1) = 0.5_dp*(1.0_dp/aux1) + 0.5_dp*sqrt((1.0_dp + 4.0_dp*aux1**2*y0(:,2)**2)/(aux1**2)) + 1e-13_dp
        ! end where

        ! Zero all values outside of star
        where (patch_tmp%cgrid%coordinates_x>y0_bdy(2) )
                y0(:,1) = 0.0_dp
                y0(:,2) = 0.0_dp
        end where
        
        ! Update current timelevel array
        patch_tmp%tl_current%soln = y0(:,:)
        patch_tmp%x_boundary = y0_bdy

        

    end subroutine rk3



    !#######################################
    !# Forward Euler (RK1)
    !#######################################

    subroutine rk1(patch_tmp,dt)
        implicit none
        type(patch), intent(inout) :: patch_tmp
        real(dp), intent(in) :: dt
        real(dp) :: t0

        t0 = patch_tmp%tl_current%t
        y0(:,:) = patch_tmp%tl_current%soln
        y0_bdy = patch_tmp%x_boundary

        call patch_tmp%system%apply_BCs(patch_tmp%cgrid,y0) ! Update External Boundary Conditions

        ! First, we evolve the boundaries of the star
        call patch_tmp%system%compute_boundary_rhs(t0,y0,y0_bdy,patch_tmp%cgrid,k1_bdy)
        y1_bdy = y0_bdy + dt*k1_bdy
        ! Update the star ghost points
        call patch_tmp%system%apply_boundary_BCs(patch_tmp%cgrid,t0,y0)
        call patch_tmp%system%compute_rhs(t0,y0(:,:),patch_tmp%cgrid,k1(:,:))
        y1(:,:) = y0(:,:) + dt*k1(:,:)

        call patch_tmp%system%apply_BCs(patch_tmp%cgrid,y1) ! Update External Boundary Conditions

        !Enforce Floor on Conserved Variables
        ! where (y1(:,1) < 0.5_dp*aux1 + 0.5_dp*sqrt((1.0_dp + 4.0_dp*aux1**2*y1(:,2)**2)/(aux1**2)))
        !         y1(:,1) = 0.5_dp*aux1 + 0.5_dp*sqrt((1.0_dp + 4.0_dp*aux1**2*y1(:,2)**2)/(aux1**2)) + 1e-13_dp
        ! end where

        ! Zero all values outside of star
        where (patch_tmp%cgrid%coordinates_x>y1_bdy(2) )
                y1(:,1) = 0.0_dp
                y1(:,2) = 0.0_dp
        end where
        
        ! Update current timelevel array
        patch_tmp%tl_current%soln = y1(:,:)
        patch_tmp%x_boundary = y1_bdy

        

    end subroutine rk1


end module timestepping_mod