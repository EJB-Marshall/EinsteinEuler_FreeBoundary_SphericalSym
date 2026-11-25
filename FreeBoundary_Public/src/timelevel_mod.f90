module timelevel_mod
    use real_type_mod, only: dp
    implicit none

    !!################################ TIMELEVEL MOD ###############################################
    !!                       
    !! A timelevel stores the solution data at a single time
    !! 
    !! For the cubed sphere, the solution data is stored in the local patch coordinates!
    !!#########################################################################################


    type :: timelevel

    real(dp) :: t
    real(dp), dimension(:,:), allocatable :: soln

    end type timelevel

end module timelevel_mod