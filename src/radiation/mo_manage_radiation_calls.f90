! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_manage_radiation_calls

  ! This module manages the calls to the radiation routines.

  use mo_definitions,      only: wp,t_grid,t_state,t_diag
  use mo_grid_nml,         only: n_cells,n_layers,n_levels
  use mo_rad_nml,          only: n_cells_rad,n_rad_blocks,rad_config,n_no_cond_rad_layers
  use mo_constituents_nml, only: lmoist,n_condensed_constituents,n_constituents
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
    integer               :: rad_block_index ! radiation block index (for OMP parallelization)
    real(wp), allocatable :: lat_scal(:)
    real(wp), allocatable :: lon_scal(:)
    real(wp), allocatable :: sfc_sw_in(:)
    real(wp), allocatable :: sfc_lw_out(:)
    real(wp), allocatable :: sfc_albedo(:)
    real(wp), allocatable :: temp_sfc(:)
    real(wp), allocatable :: z_scal(:,:)
    real(wp), allocatable :: z_vect(:,:)
    real(wp), allocatable :: rho(:,:,:)
    real(wp), allocatable :: temp(:,:)
    real(wp), allocatable :: rad_tend(:,:)
    
    if (rad_config==1) then
      write(*,*) "Starting update of radiative fluxes ..."
    endif
    
    ! loop over all radiation blocks
    !$omp parallel do private(rad_block_index,lat_scal,lon_scal,sfc_sw_in,sfc_lw_out,sfc_albedo, &
    !$omp temp_sfc,z_scal,z_vect,rho,temp,rad_tend)
    do rad_block_index=1,n_rad_blocks
    
      ! allocating memory for the arrays that will be handed over to the radiation subroutine
      allocate(lat_scal(n_cells_rad))
      allocate(lon_scal(n_cells_rad))
      allocate(sfc_sw_in(n_cells_rad))
      allocate(sfc_lw_out(n_cells_rad))
      allocate(sfc_albedo(n_cells_rad))
      allocate(temp_sfc(n_cells_rad))
      allocate(z_scal(n_cells_rad,n_layers))
      allocate(z_vect(n_cells_rad,n_levels))
      allocate(rho(n_cells_rad,n_layers,n_constituents))
      allocate(temp(n_cells_rad,n_layers))
      allocate(rad_tend(n_cells_rad,n_layers))
      
      ! initializing all radiation-specific arrays with zeros
      lat_scal = 0._wp
      lon_scal = 0._wp
      sfc_sw_in = 0._wp
      sfc_lw_out = 0._wp
      sfc_albedo = 0._wp
      temp_sfc = 0._wp
      z_scal = 0._wp
      z_vect = 0._wp
      rho = 0._wp
      temp = 0._wp
      rad_tend = 0._wp
    
      ! remapping all the arrays
      call create_rad_array_scalar_h(grid%lat_c,lat_scal,rad_block_index)
      call create_rad_array_scalar_h(grid%lon_c,lon_scal,rad_block_index)
      call create_rad_array_scalar_h(state%temperature_soil(:,1),temp_sfc,rad_block_index)
      call create_rad_array_scalar_h(grid%sfc_albedo,sfc_albedo,rad_block_index)
      call create_rad_array_scalar(grid%z_scalar,z_scal,rad_block_index)
      call create_rad_array_vector(grid%z_vector_v,z_vect,rad_block_index)
      call create_rad_array_mass_den(state%rho,rho,rad_block_index)
      call create_rad_array_scalar(diag%temperature,temp,rad_block_index)
      
      ! calling the radiation routine
      ! RTE+RRTMGP
      if (rad_config==1) then
        ! this is necessary for stability for now
        if (lmoist) then
          rho(:,1:n_no_cond_rad_layers,1:n_condensed_constituents) = 0._wp
        endif
        call calc_radiative_flux_convergence(lat_scal,lon_scal,z_scal,z_vect,rho,temp,rad_tend,temp_sfc, &
                                             sfc_sw_in,sfc_lw_out,sfc_albedo,time_coordinate)
      endif
      ! Held-Suarez
      if (rad_config==2) then
        call held_suar(lat_scal,rho,temp,rad_tend)
      endif
      
      ! filling the actual radiation tendency
      call remap_to_original(rad_tend,diag%radiation_tendency,rad_block_index)
      call remap_to_original_scalar_h(sfc_sw_in,diag%sfc_sw_in,rad_block_index)
      call remap_to_original_scalar_h(sfc_lw_out,diag%sfc_lw_out,rad_block_index)
      
      deallocate(lat_scal)
      deallocate(lon_scal)
      deallocate(sfc_sw_in)
      deallocate(sfc_lw_out)
      deallocate(sfc_albedo)
      deallocate(temp_sfc)
      deallocate(z_scal)
      deallocate(z_vect)
      deallocate(rho)
      deallocate(temp)
      deallocate(rad_tend)
      
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
















