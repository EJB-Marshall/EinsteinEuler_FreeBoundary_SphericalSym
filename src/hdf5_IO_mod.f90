module hdf5_IO_mod
    use real_type_mod, only: dp
    use patch_mod
    use read_parameters_mod, only: aux1, aux2, aux3 
    use hdf5
    implicit none

    integer(hid_t) :: file_id, time_dset_id, boundary_dset_id
    integer(hid_t) :: patch_group_ids, soln_dset_ids
    integer :: error
    integer :: nt_written = 0 ! Number of timesteps saved to file


    ! Note: error *should* return a *negative* value when a call is unsuccessful
    ! This code should be updated to clarify what all the HDF routines do, add error checks etc.
    ! For further details on the HDF library/functions see: https://davis.lbl.gov/Manuals/HDF5-1.6.5/RM_H5S.html

    contains 

    !####################################################
    !# Create HDF5 File, Groups etc.
    !####################################################
    subroutine initialise_hdf(filename, patches, buffer_size)
        implicit none
        character(len=*), intent(in) :: filename
        type(patch), intent(in) :: patches
        integer, intent(in) :: buffer_size

        integer :: Nx, Nvar
        integer(hsize_t), dimension(3) :: dims_soln, maxdims_soln, chunk_dims_soln
        integer(hsize_t), dimension(1) :: dims_coords, grid_dims
        integer(hsize_t), dimension(1) :: dims_time, maxdims_time, chunk_dims_time
        integer(hsize_t), dimension(2) :: dims_boundary, maxdims_boundary, chunk_dims_boundary
        integer(hid_t) :: plist_soln, plist_time,  plist_boundary  ! Property list for chunking
        integer(hid_t) :: soln_dspace_id, time_dspace, boundary_dspace
        integer(hid_t) :: comp_grid_dspace_id, comp_grid_dset_id
        integer(hid_t) :: aux_dspace_id, aux_dset_id
        character(len=256) :: folder
        ! character(len=64) :: patch_name
        real(dp) :: aux(3)


        nt_written = 0 ! Initialise variable

        ! Get dimensions of each patch grid
        Nx = patches%cgrid%Nx
        Nvar = patches%system%Nvar

        ! First we initialise the Fortran boundary
        call h5open_f(error)

        ! Next, we create the HDF5 file
        folder = 'HDF_Files/'

        ! Create the HDF file
        call h5fcreate_f(trim(folder)//trim(filename), H5F_ACC_TRUNC_F, file_id, error)


        !##################################################
        !# Create groups and datasets for each patch
        !##################################################

        ! Dimension properties for dataspace of dataset
        dims_soln = (/ Nx, Nvar, 0 /)
        maxdims_soln = (/ int(Nx,8), int(Nvar,8),H5S_UNLIMITED_F /)
        chunk_dims_soln = (/ Nx, Nvar, buffer_size /)
        dims_coords = (/ Nx /)


        ! Create property list for chunking (required for unlimited dimensions)
        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_soln, error)
        call h5pset_chunk_f(plist_soln, 3, chunk_dims_soln, error)
        call h5pset_deflate_f(plist_soln, 5, error) ! Gzip compression (0 = no compression, 9 = most compression)


        ! Create groups
        call h5gcreate_f(file_id, trim('Solution'), patch_group_ids, error)

        ! Create data space for the solution arrays on each patch
        call h5screate_simple_f(3, dims_soln, soln_dspace_id, error, maxdims_soln)
        call h5dcreate_f(patch_group_ids, "Solution_Array", H5T_NATIVE_DOUBLE, &
                        soln_dspace_id, soln_dset_ids, error, plist_soln)
        call h5sclose_f(soln_dspace_id, error)
        

        ! Close the solution property list
        call h5pclose_f(plist_soln, error)


        !############################################
        !# Create fluid boundary dataset
        !############################################

        dims_boundary = (/2, 0 /)
        maxdims_boundary = (/2_8, H5S_UNLIMITED_F /)
        chunk_dims_boundary= (/ 2,buffer_size /)

        ! Create data space for the boundary data
        call h5screate_simple_f(2, dims_boundary, boundary_dspace, error, maxdims_boundary)

        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_boundary, error)
        call h5pset_chunk_f(plist_boundary, 2, chunk_dims_boundary, error)
        call h5pset_deflate_f(plist_boundary, 6, error) 

        ! Create the boundary dataset 
        call h5dcreate_f(file_id, "Fluid_boundaries", H5T_NATIVE_DOUBLE, boundary_dspace, boundary_dset_id, error, plist_boundary)

        ! Close the dataspace
        call h5sclose_f(boundary_dspace, error)

        ! Close property list
        call h5pclose_f(plist_boundary, error)


        !##################################################
        !# Create time dataset
        !##################################################
        
        dims_time = (/ 0 /)
        maxdims_time = (/ H5S_UNLIMITED_F /)
        chunk_dims_time = (/ buffer_size /)

        ! Create data space for the time data
        call h5screate_simple_f(1, dims_time, time_dspace, error, maxdims_time)

        call h5pcreate_f(H5P_DATASET_CREATE_F, plist_time, error)
        call h5pset_chunk_f(plist_time, 1, chunk_dims_time, error)
        call h5pset_deflate_f(plist_time, 6, error) 

        ! Create the time dataset (Note this is *not* in a group)
        call h5dcreate_f(file_id, "Time", H5T_NATIVE_DOUBLE, time_dspace, time_dset_id, error, plist_time)

        ! Close the dataspace
        call h5sclose_f(time_dspace, error)

        ! Close property list
        call h5pclose_f(plist_time, error)

        !##################################################
        !# Store Computational Grid
        !##################################################
        ! Create dataspace for computational grid - since each patch has the same comp coordinates, only save this once
        grid_dims = (/ int(Nx,8)  /)
        call h5screate_simple_f(1, grid_dims, comp_grid_dspace_id, error)
        
        ! Create grid datasets, write data, then close dataset
        call h5dcreate_f(file_id, "x_coordinates", H5T_NATIVE_DOUBLE, comp_grid_dspace_id, &
                        comp_grid_dset_id, error)
        call h5dwrite_f(comp_grid_dset_id, H5T_NATIVE_DOUBLE, patches%cgrid%coordinates_x(1:Nx), grid_dims, error)
        call h5dclose_f(comp_grid_dset_id, error)

        

        !##################################################
        !# Store Auxiliary Parameters
        !##################################################
        ! Create vector with auxiliary quantities
        aux(1) = aux1
        aux(2) = aux2
        aux(3) = aux3

        ! Create dataspace for aux vector
        call h5screate_simple_f(1, (/ 3_8 /), aux_dspace_id, error)
        
        ! Create aux dataset 
        call h5dcreate_f(file_id, "Aux", H5T_NATIVE_DOUBLE, aux_dspace_id, &
                        aux_dset_id, error)
        
        ! Write aux vector 
        call h5dwrite_f(aux_dset_id, H5T_NATIVE_DOUBLE, aux, (/ 3_8 /), error)
        
        ! Close datasets
        call h5dclose_f(aux_dset_id, error)
        call h5sclose_f(aux_dspace_id, error)

    end subroutine initialise_hdf

    !####################################################
    !# Write solution buffer to HDF file
    !####################################################
    subroutine write_to_hdf(soln_buffer,time_buffer,boundary_buffer,buffer_count,Nx,Nvar)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: soln_buffer !(Nx, Nvar, buffer_count)
        real(dp), dimension(:), intent(in) :: time_buffer
        real(dp), dimension(:,:), intent(in) :: boundary_buffer
        integer, intent(in) :: buffer_count, Nx,  Nvar

        ! Local variables
        integer :: hstat
        integer(hsize_t), dimension(3) :: new_dims_soln, start_soln, count_soln
        integer(hsize_t), dimension(2) :: new_dims_boundary, start_boundary, count_boundary
        integer(hsize_t), dimension(1) :: new_dims_time, start_time, count_time
        integer(hid_t) :: memspace_soln, filespace_soln, memspace_time, filespace_time, &
                        memspace_boundary, filespace_boundary

        !############################
        !# Write Solution Data
        !###########################

        new_dims_soln  = (/ int(Nx,8), int(Nvar,8), int(nt_written + buffer_count,8) /) ! Append new data after the last set of buffered timesteps
        start_soln  = (/ 0_8, 0_8,  int(nt_written,8) /)         
        count_soln  = (/ int(Nx,8), int(Nvar,8), int(buffer_count,8) /)
        
        ! Memory space must match buffer shape
        call h5screate_simple_f(3, count_soln, memspace_soln, error)

        ! Extend dataset
        call h5dset_extent_f(soln_dset_ids, new_dims_soln, error)

        ! Select region in file
        call h5dget_space_f(soln_dset_ids, filespace_soln, error)
        call h5sselect_hyperslab_f(filespace_soln, H5S_SELECT_SET_F, start_soln, count_soln, error)

        ! Write the whole (Nx, buffer_count) block
        call h5dwrite_f(soln_dset_ids, H5T_NATIVE_DOUBLE, soln_buffer(:,:,1:buffer_count), &
                        count_soln, error, memspace_soln, filespace_soln)

        call h5sclose_f(filespace_soln, error)

        call h5sclose_f(memspace_soln, error)



        !############################
        !# Write boundary Data
        !###########################
        new_dims_boundary = (/ 2_8,int(nt_written + buffer_count,8) /)
        start_boundary = (/ 0_8,int(nt_written,8) /)
        count_boundary = (/ 2_8,int(buffer_count,8) /)

        ! Save location of fluid boundarys
        call h5dset_extent_f(boundary_dset_id, new_dims_boundary, hstat)
        call h5dget_space_f(boundary_dset_id, filespace_boundary, hstat)
        call h5sselect_hyperslab_f(filespace_boundary, H5S_SELECT_SET_F, start_boundary, count_boundary, hstat)

        call h5screate_simple_f(2, count_boundary, memspace_boundary, hstat)
        call h5dwrite_f(boundary_dset_id, H5T_NATIVE_DOUBLE, boundary_buffer(:,1:buffer_count), &
                        count_boundary, hstat, mem_space_id=memspace_boundary, file_space_id=filespace_boundary)

        call h5sclose_f(memspace_boundary, hstat)
        call h5sclose_f(filespace_boundary, hstat)

        
        !############################
        !# Write Time Data
        !###########################
        new_dims_time = (/ int(nt_written + buffer_count,8) /)
        start_time  = (/ int(nt_written,8) /)
        count_time  = (/ int(buffer_count,8) /)

        call h5dset_extent_f(time_dset_id, new_dims_time, hstat)
        call h5dget_space_f(time_dset_id, filespace_time, hstat)
        call h5sselect_hyperslab_f(filespace_time, H5S_SELECT_SET_F, start_time, count_time, hstat)

        call h5screate_simple_f(1, count_time, memspace_time, hstat)
        call h5dwrite_f(time_dset_id, H5T_NATIVE_DOUBLE, time_buffer(1:buffer_count), &
                        count_time, hstat, mem_space_id=memspace_time, file_space_id=filespace_time)

        call h5sclose_f(memspace_time, hstat)
        call h5sclose_f(filespace_time, hstat)


        ! Advance the write pointer
        nt_written = nt_written + buffer_count

    end subroutine write_to_hdf


    !####################################################
    !# Close HDF File
    !####################################################
    subroutine close_hdf()
        implicit none
    
        ! Close all datasets
        call h5dclose_f(soln_dset_ids, error)
        call h5gclose_f(patch_group_ids, error)

        ! Close boundary dataset
        call h5dclose_f(boundary_dset_id, error)

        ! Close time dataset
        call h5dclose_f(time_dset_id, error)
        
        ! Close the file
        call h5fclose_f(file_id, error)
        
        ! Close the Fortran boundary
        call h5close_f(error)

    end subroutine close_hdf


end module hdf5_IO_mod