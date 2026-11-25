module einstein_euler_initial_data_mod
    use real_type_mod, only: dp
    use patch_mod
    use read_parameters_mod, only: aux1, aux2, aux3
    implicit none

    contains

    subroutine initialise_euler_data(patches,t0)
        implicit none
        type(patch), intent(inout) :: patches
        real(dp), intent(in) :: t0
        integer :: i, Nx, Ngz
        real(dp) :: K, dx, a
        real(dp), allocatable :: v(:), mu(:), Gamma2(:)
        real(dp), parameter :: pi = acos(-1.0_dp)

        Nx = patches%cgrid%Nx
        Ngz = patches%cgrid%Ngz
        dx = patches%cgrid%dx

        K = aux1

        patches%x_boundary(1) = 0.0_dp
        patches%x_boundary(2) = 0.0_dp

        patches%tl_current%soln = 0.0_dp



        !-----------------------------------------------
        ! Load Perturbed TOV Initial Data
        !-----------------------------------------------

        allocate(v(1-Ngz:Nx+Ngz),mu(1-Ngz:Nx+Ngz),Gamma2(1-Ngz:Nx+Ngz))

        open(unit=2,file="rho_perturbed_ID.bin",access="stream",status="old",action="read")
        read(2) mu
        close(2)

        open(unit=2,file="v_perturbed_ID.bin",access="stream",status="old",action="read")
        read(2) v
        close(2)

        Gamma2 = 1.0_dp/(1.0_dp - v**2)

        ! Define tau and S
        patches%tl_current%soln(:,1) = ((K + 1.0_dp)*mu - 1.0_dp)*Gamma2 - (K*mu- 1.0_dp)
        patches%tl_current%soln(:,2) = ((K + 1.0_dp)*mu - 1.0_dp)*Gamma2*v
        
        ! Metric Function a
        open(unit=2,file="a_perturbed_ID.bin",access="stream",status="old",action="read")
        read(2) patches%tl_current%soln(:,3)
        close(2)

        ! Star boundary
        open(unit=2,file="star_bdy_perturbed_ID.bin",access="stream",status="old",action="read")
        read(2) patches%x_boundary(2)
        close(2)

        ! Apply floor outside of star
        where (patches%cgrid%coordinates_x > patches%x_boundary(2))

            patches%tl_current%soln(:,1) = 1e-10_dp
            patches%tl_current%soln(:,2) = 0.0_dp
            
            

        end where


        deallocate(v,mu,Gamma2)



        ! Apply boundary conditions
        call patches%system%apply_BCs(patches%cgrid,patches%tl_current%soln)


    end subroutine



end module einstein_euler_initial_data_mod