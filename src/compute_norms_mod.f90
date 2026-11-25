module compute_norms_mod
    use real_type_mod, only: dp
    use patch_mod
    use read_parameters_mod, only: Nx, Ngz

    contains

    
    !####################################################
    ! Compute the discrete L2 norm at a single timestep
    !####################################################
    subroutine compute_l2_norm(soln,dx,norm)
        implicit none
        real(dp), intent(in) :: soln(1-Ngz:Nx+Ngz,2) 
        real(dp), intent(in) :: dx
        real(dp), intent(out) :: norm(2)
        integer :: i

        ! allocate(norm(size(soln,dim=2))) ! norm has size Nvar

        norm = 0.0_dp
        do i = 1, Nx

            norm = norm + (dx*soln(i,:)**(2))

        end do

        norm = sqrt(norm)

    end subroutine compute_l2_norm



    !#############################################################
    ! Compute the discrete L Infinity norm at a single timestep
    !###########################################################
    subroutine compute_linf_norm(soln,norm)
        implicit none
        real(dp), intent(in) :: soln(1-Ngz:Nx+Ngz,2) 
        real(dp), intent(out) :: norm(2)

        ! allocate(norm(size(soln,dim=2))) ! norm has size Nvar

        norm = 0.0_dp
        
        norm = maxval(abs(soln),dim=1)

    end subroutine compute_linf_norm

    

end module compute_norms_mod