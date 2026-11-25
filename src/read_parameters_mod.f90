module read_parameters_mod
    use real_type_mod, only: dp
    implicit none
    integer :: Nx, Ngz
    real(dp) :: t0, tend, aux1, aux2, aux3 
    character(len=256) :: filename

    public t0, tend, Nx, Ngz, read_params, filename, aux1, aux2, aux3
    private 

    contains

    subroutine read_params 
        implicit none

        ! Load in parameters from namelist
        namelist /Grid_Params/ Nx, Ngz
        namelist /Sim_Params/ t0, tend, filename
        namelist /Aux_Params/ aux1, aux2, aux3

        open(file='params.nml',status='old',unit=10,action='read')
        read(10,nml=Grid_Params)
        read(10,nml=Sim_Params)
        read(10,nml=Aux_Params)
        close(10)

    
    end subroutine read_params


end module read_parameters_mod