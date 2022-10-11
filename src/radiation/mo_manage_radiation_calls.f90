! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_manage_radiation_calls

  ! This module manages the calls to the radiation routines.

  use mo_definitions,      only: wp,t_grid,t_state,t_diag,t_radiation
  use mo_grid_nml,         only: n_cells,n_layers,n_levels
  use mo_rad_nml,          only: n_cells_rad,n_rad_blocks,rad_config
  use mo_constituents_nml, only: n_constituents
  use mo_surface_nml,      only: nsoillays
  use mo_held_suarez,      only: held_suar
  use mo_rrtmgp_coupler,   only: calc_radiative_flux_convergence

  implicit none
  
  contains
  
  subroutine update_rad_fluxes(state,diag,grid,time_coordinate)
    
    ! This subroutine manages the calls to RTE+RRTMGP and the Held-Suarez interface.
  
    type(t_state), intent(in)    :: state           ! state variables to use for updating the radiative fluxes
    type(t_diag),  intent(inout) :: diag            ! diagnostic quantities to use for updating the radiative fluxes
    type(t_grid),  intent(in)    :: grid            ! grid quantities
    real(wp),      intent(in)    :: time_coordinate ! epoch timestamp (needed for computing the zenith angle)
  
    ! local variables
    integer           :: rad_block_index ! radiation block index (for OMP parallelization)
    type(t_radiation) :: radiation       ! radiation state
    
    if (rad_config==1) then
      write(*,*) "Starting update of radiative fluxes ..."
    endif
    
    ! loop over all radiation blocks
    !$omp parallel do private(rad_block_index,radiation)
    do rad_block_index=1,n_rad_blocks
    
      allocate(radiation%lat_scal(n_cells_rad))
      allocate(radiation%lon_scal(n_cells_rad))
      allocate(radiation%sfc_sw_in(n_cells_rad))
      allocate(radiation%sfc_lw_out(n_cells_rad))
      allocate(radiation%sfc_albedo(n_cells_rad))
      allocate(radiation%temp_sfc(n_cells_rad))
      allocate(radiation%z_scal(n_cells_rad,n_layers))
      allocate(radiation%z_vect(n_cells_rad,n_levels))
      allocate(radiation%rho(n_cells_rad,n_layers,n_constituents))
      allocate(radiation%temp(n_cells_rad,n_layers))
      allocate(radiation%rad_tend(n_cells_rad,n_layers))
      
      radiation%lat_scal = 0._wp
      radiation%lon_scal = 0._wp
      radiation%sfc_sw_in = 0._wp
      radiation%sfc_lw_out = 0._wp
      radiation%sfc_albedo = 0._wp
      radiation%temp_sfc = 0._wp
      radiation%z_scal = 0._wp
      radiation%z_vect = 0._wp
      radiation%rho = 0._wp
      radiation%temp = 0._wp
      radiation%rad_tend = 0._wp
    
      ! remapping all the arrays
      call create_rad_array_scalar_h(grid%lat_c,radiation%lat_scal,rad_block_index)
      call create_rad_array_scalar_h(grid%lon_c,radiation%lon_scal,rad_block_index)
      call create_rad_array_scalar_h(state%temperature_soil(:,1),radiation%temp_sfc,rad_block_index)
      call create_rad_array_scalar_h(grid%sfc_albedo,radiation%sfc_albedo,rad_block_index)
      call create_rad_array_scalar(grid%z_scalar,radiation%z_scal,rad_block_index)
      call create_rad_array_vector(grid%z_vector_v,radiation%z_vect,rad_block_index)
      call create_rad_array_mass_den(state%rho,radiation%rho,rad_block_index)
      call create_rad_array_scalar(diag%temperature,radiation%temp,rad_block_index)
      ! calling the radiation routine
      ! RTE+RRTMGP
      if (rad_config==1) then
        
        call calc_radiative_flux_convergence(radiation%lat_scal,radiation%lon_scal,radiation%z_scal,radiation%z_vect, &
                                             radiation%rho,radiation%temp,radiation%rad_tend,radiation%temp_sfc, &
                                             radiation%sfc_sw_in,radiation%sfc_lw_out,radiation%sfc_albedo,time_coordinate)
      endif
      ! Held-Suarez
      if (rad_config==2) then
        call held_suar(radiation%lat_scal,radiation%rho,radiation%temp,radiation%rad_tend)
      endif
      
      ! filling the actual radiation tendency
      call remap_to_original(radiation%rad_tend,diag%radiation_tendency,rad_block_index)
      call remap_to_original_scalar_h(radiation%sfc_sw_in,diag%sfc_sw_in,rad_block_index)
      call remap_to_original_scalar_h(radiation%sfc_lw_out,diag%sfc_lw_out,rad_block_index)
      
      deallocate(radiation%lat_scal)
      deallocate(radiation%lon_scal)
      deallocate(radiation%sfc_sw_in)
      deallocate(radiation%sfc_lw_out)
      deallocate(radiation%sfc_albedo)
      deallocate(radiation%temp_sfc)
      deallocate(radiation%z_scal)
      deallocate(radiation%z_vect)
      deallocate(radiation%rho)
      deallocate(radiation%temp)
      deallocate(radiation%rad_tend)
      
    enddo
    !$omp end parallel do
    
    if (rad_config==1) then
      write(*,*) "Update of radiative fluxes completed."
    endif
  
  end subroutine update_rad_fluxes
  
  subroutine create_rad_array_scalar(in_array,out_array,rad_block_index)

    ! This subroutine cuts out a slice of a scalar field for hand-over to the radiation routine (done for RAM efficiency reasons).
    
    real(wp), intent(in)  :: in_array(n_cells,n_layers)
    real(wp), intent(out) :: out_array(n_cells_rad,n_layers)
    integer               :: rad_block_index
    
    ! local variables
    integer :: ji ! horizontal loop index
    
    ! loop over all elements of the resulting array
    do ji=1,n_cells_rad
      out_array(ji,:) = in_array((rad_block_index-1)*n_cells_rad+ji,:)
    enddo
  
  end subroutine create_rad_array_scalar
  
  subroutine create_rad_array_scalar_h(in_array,out_array,rad_block_index)

    ! This subroutine cuts out a slice of a horizontal scalar field for hand-over to the radiation routine (done for RAM efficiency reasons).
    
    real(wp), intent(in)  :: in_array(n_cells)
    real(wp), intent(out) :: out_array(n_cells_rad)
    integer               :: rad_block_index
    
    ! local variables
    integer :: ji ! horizontal loop index
    
    ! loop over all elements of the resulting array
    do ji=1,n_cells_rad
      out_array(ji) = in_array((rad_block_index-1)*n_cells_rad+ji)
    enddo
  
  end subroutine create_rad_array_scalar_h

  subroutine create_rad_array_vector(in_array,out_array,rad_block_index)

    ! This subroutine cuts out a slice of a vector field for hand-over to the radiation routine (done for RAM efficiency reasons).
  ! Only the vertical vector points are taken into account since only they are needed by the radiation.

    real(wp), intent(in)  :: in_array(n_cells,n_levels)
    real(wp), intent(out) :: out_array(n_cells_rad,n_levels)
    integer               :: rad_block_index
    
    ! local variables
    integer :: ji ! horizontal loop index
    
    ! loop over all elements of the resulting array
    do ji=1,n_cells_rad
      out_array(ji,:) = in_array((rad_block_index-1)*n_cells_rad+ji,:)
    enddo
  
  end subroutine create_rad_array_vector
  
  subroutine create_rad_array_mass_den(in_array,out_array,rad_block_index)

    ! This subroutine does same thing as create_rad_array_scalar,only for a mass density field.
    
    real(wp), intent(in)  :: in_array(n_cells,n_layers,n_constituents)
    real(wp), intent(out) :: out_array(n_cells_rad,n_layers,n_constituents)
    integer               :: rad_block_index
    
    ! local variables
    integer :: ji
    
    ! loop over all cells of the resulting array
    do ji=1,n_cells_rad
      out_array(ji,:,:) = in_array((rad_block_index-1)*n_cells_rad+ji,:,:)
    enddo
  
  end subroutine create_rad_array_mass_den

  subroutine remap_to_original(in_array,out_array,rad_block_index)

    ! This subroutine reverses what create_rad_array_scalar has done.
    
    real(wp), intent(in)  :: in_array(n_cells_rad,n_layers)
    real(wp), intent(out) :: out_array(n_cells,n_layers)
    integer               :: rad_block_index
    
    ! local variables
    integer :: ji
    
    ! loop over all cells of the resulting array
    do ji=1,n_cells_rad
      out_array((rad_block_index-1)*n_cells_rad+ji,:) = in_array(ji,:)
    enddo
  
  end subroutine remap_to_original

  subroutine remap_to_original_scalar_h(in_array,out_array,rad_block_index)

    ! This subroutine reverses what create_rad_array_scalar_h has done.
    
    real(wp), intent(in)  :: in_array(n_cells_rad)
    real(wp), intent(out) :: out_array(n_cells)
    integer               :: rad_block_index
    
    ! local variables
    integer :: ji
    
    ! loop over all elements of the resulting array
    do ji=1,n_cells_rad
      out_array((rad_block_index-1)*n_cells_rad+ji) = in_array(ji)
    enddo
  
  
  end subroutine remap_to_original_scalar_h

end module mo_manage_radiation_calls
















