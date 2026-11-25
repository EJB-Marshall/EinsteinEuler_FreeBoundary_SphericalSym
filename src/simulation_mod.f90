module simulation_mod
    use real_type_mod, only: dp
    use patch_mod
    use timestepping_mod, only: initialise_rk_arrays, rk4, rk1, rk3, rk4_atmosphere
    use read_parameters_mod, only: Nx, Ngz
    use hdf5_IO_mod
    use compute_norms_mod
    use diff_op_mod
    implicit none

    ! public
    ! private


    type :: simulation 

    real(dp) :: CFL, t0, tend
    type(patch) :: patches
    character(len=256) :: filename
    

    contains 
    
    procedure :: evolve
    

    end type simulation


    contains

    subroutine evolve(this)
        implicit none
        class(simulation), intent(inout) :: this
        integer :: Nvar 
        real(dp) :: dt, tsave, tinterval

        !# File and Storage Buffer Variables
        integer, parameter :: buffer_size = 1000
        real(dp), allocatable :: solution_buffer(:,:,:)
        real(dp) :: time_buffer(buffer_size)
        real(dp) :: boundary_buffer(2,buffer_size)
        integer :: buffer_index = 0
        ! Timing variables
        integer :: beginning
        real(dp) :: rate

        Nx = this%patches%cgrid%Nx
        Nvar = this%patches%system%Nvar

        !-----------------------------
        ! HDF Setup
        !-----------------------------

        ! Allocate solution buffer and save initial data
        allocate(solution_buffer(Nx, Nvar, buffer_size))
        buffer_index = 1
        solution_buffer(:,:, buffer_index) = this%patches%tl_current%soln(1:Nx,:)
        boundary_buffer(:,buffer_index) = this%patches%x_boundary

        ! Save initial time
        time_buffer(buffer_index) = this%t0


        ! Adjust how frequently we save data
        tsave = 0.1_dp ! First time to save
        tinterval = 0.1_dp ! Save every tinterval increments of coordinate time

        ! step_count = 10 ! Save every step_count amount of time steps

        ! Initialise HDF5 File
        call initialise_hdf(trim(this%filename),this%patches,buffer_size)

        
        !-----------------------------
        ! Set initial timestep
        !-----------------------------

        ! Set timestep - This needs to be updated to account for characteristic speeds etc.
        dt = this%CFL * this%patches%cgrid%dx ! Evolve towards future
        ! dt = -this%CFL * this%patches%cgrid%dx ! Evolve towards past
    
        !-----------------------------
        ! Main Evolution Loop
        !-----------------------------

        !Timing code
        call system_clock(beginning, rate)

        ! Initialise Runge-Kutta arrays
        call initialise_rk_arrays(this%patches)

        do while (this%t0 < this%tend) ! Future Evolution
        ! do while (this%t0 > this%tend) ! Past Evolution


            dt = min(this%CFL*this%patches%cgrid%dx/this%patches%system%CS_x, &
                    this%CFL*this%patches%cgrid%dx) ! To the future 

            ! dt = -min(this%CFL*this%patches%cgrid%dx/this%patches%system%CS_x, &
            !         this%CFL*this%patches%cgrid%dx) ! To the past 


            ! dt = this%CFL*this%patches%cgrid%dx



            !------ Force code to finish at correct time --------
        
            ! Future        
            if (this%t0 + dt > this%tend) then
                dt = this%tend - this%t0
            end if

            ! Force code to save every tsave interval
            if (this%t0 + dt > tsave) then
                dt = tsave - this%t0
            end if

            ! Past
            ! if (this%t0 + dt < this%tend) then
            !     dt = this%tend - this%t0
            ! end if

            ! ! Force code to save every tsave interval
            ! if (this%t0 + dt < tsave) then
            !     dt = tsave - this%t0
            ! end if


            ! Update solution using Runge Kutta routine
            call rk4(this%patches,dt)

            ! Update time after a single timestep
            ! This needs to be cleaned up! Don't need to store time in simulation and patches (why is each patch storing time twice?)
            this%t0 = this%t0 + dt
            this%patches%tl_current%t = this%patches%tl_current%t + dt
            this%patches%t = this%patches%t + dt

            
            !# Future Saving
            if (this%t0 >= tsave) then 
                buffer_index = buffer_index + 1
                solution_buffer(:,:, buffer_index) = this%patches%tl_current%soln(1:Nx,:)
                boundary_buffer(:, buffer_index) = this%patches%x_boundary

                time_buffer(buffer_index) = this%t0
                tsave = tsave + tinterval
                ! print *, "t is", this%t0
                ! print *, "dt is", dt
                ! print *, " "

                
                ! Save buffer array to file when full
                if (buffer_index == buffer_size) then
                    call write_to_hdf(solution_buffer,time_buffer,boundary_buffer,buffer_index,Nx,Nvar)
                    buffer_index = 0
                end if
                
            end if


            ! !# Past Saving
            ! if (this%t0 <= tsave) then 
            !     buffer_index = buffer_index + 1
            !     solution_buffer(:,:, buffer_index) = this%patches%tl_current%soln(1:Nx,:)

            !     time_buffer(buffer_index) = this%t0
            !     tsave = tsave - tinterval
            !     ! print *, "t is", this%t0
            !     ! print *, "dt is", dt
            !     ! print *, " "

            !     ! Save buffer array to file when full
            !     if (buffer_index == buffer_size) then
            !         call write_to_hdf(solution_buffer,time_buffer,boundary_buffer,buffer_index,Nx,Nvar)
            !         buffer_index = 0
            !     end if
                
            ! end if



        end do

        ! Update CSV if shock/stability outcome not reached 
        ! if (outcome == 3) then
        !             print *, "Simulation finished without stability outcome!"

        !             open(unit=10, file="parameter_sweep.csv", status="unknown", &
        !                 action="write", position="append")
        !             write(10, '(i0, ",", f10.6, ",", f10.6, ",", f10.6, ",", f10.6, ",", i0)') &
        !                     Nx, aux1, aux2, aux3, this%t0, outcome

        !             close(10)

        ! end if

        ! Save any remaining data in buffer array
        if (buffer_index > 0) then
            call write_to_hdf(solution_buffer,time_buffer,boundary_buffer,buffer_index,Nx,Nvar)
        end if

        ! Close HDF file
        call close_hdf()
        deallocate(solution_buffer)


    end subroutine evolve



    subroutine initialise_simulation(patches,CFL,t0, tend, filename, sim_obj)
        implicit none
        type(patch), intent(in) :: patches
        type(simulation), intent(out) :: sim_obj
        real(dp) :: CFL, t0, tend
        character(len=256) :: filename
        
        sim_obj%patches = patches
        sim_obj%CFL = CFL
        sim_obj%t0 = t0
        sim_obj%tend = tend
        sim_obj%filename = filename
        

        
    end subroutine initialise_simulation
    
    
end module simulation_mod