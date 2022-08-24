! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_manage_radiation_calls

  ! This module manages the calls to the radiation routines.

  use iso_c_binding
  use mo_definitions,      only: wp,t_radiation
  use mo_grid_nml,         only: n_scalars,n_scalars_h,n_vectors_per_layer,n_vectors,n_layers
  use mo_rad_nml,          only: n_scals_rad_h,n_scals_rad,n_rad_blocks,rad_config
  use mo_constituents_nml, only: n_constituents,n_condensed_constituents
  use mo_surface_nml,      only: nsoillays
  use mo_held_suarez,      only: held_suar
  use mo_rrtmgp_coupler,   only: calc_radiative_flux_convergence

  implicit none
  
  contains
  
  subroutine call_radiation(latitude_scalar,longitude_scalar,temperature_soil,sfc_albedo,z_scalar, &
                            z_vector,rho,temperature,radiation_tendency,sfc_sw_in,sfc_lw_out,time_coordinate) &
  bind(c,name = "call_radiation")
  
    real(wp), intent(in)  :: latitude_scalar(n_scalars_h),longitude_scalar(n_scalars_h), &
                             temperature_soil(nsoillays*n_scalars_h),sfc_albedo(n_scalars_h), &
                             z_scalar(n_scalars),z_vector(n_vectors),rho(n_constituents*n_scalars), &
                             temperature(n_scalars),time_coordinate
    real(wp), intent(out) :: radiation_tendency(n_scalars),sfc_sw_in(n_scalars_h),sfc_lw_out(n_scalars_h)
  
    ! local variables
    integer           :: rad_block_index
    type(t_radiation) :: radiation
    
    if (rad_config==1) then
      write(*,*) "Starting update of radiative fluxes ..."
    endif
    ! loop over all radiation blocks
    !$omp parallel do private(rad_block_index,radiation)
    do rad_block_index=0,n_rad_blocks-1
    
      allocate(radiation%lat_scal(n_scals_rad_h))
      allocate(radiation%lon_scal(n_scals_rad_h))
      allocate(radiation%sfc_sw_in(n_scals_rad_h))
      allocate(radiation%sfc_lw_out(n_scals_rad_h))
      allocate(radiation%sfc_albedo(n_scals_rad_h))
      allocate(radiation%temp_sfc(n_scals_rad_h))
      allocate(radiation%z_scal(n_scals_rad))
      allocate(radiation%z_vect(n_scals_rad+n_scals_rad_h))
      allocate(radiation%rho(n_constituents*n_scals_rad))
      allocate(radiation%temp(n_scals_rad))
      allocate(radiation%rad_tend(n_scals_rad))
      
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
      call create_rad_array_scalar_h(latitude_scalar,radiation%lat_scal,rad_block_index)
      call create_rad_array_scalar_h(longitude_scalar,radiation%lon_scal,rad_block_index)
      call create_rad_array_scalar_h(temperature_soil,radiation%temp_sfc,rad_block_index)
      call create_rad_array_scalar_h(sfc_albedo,radiation%sfc_albedo,rad_block_index)
      call create_rad_array_scalar(z_scalar,radiation%z_scal,rad_block_index)
      call create_rad_array_vector(z_vector,radiation%z_vect,rad_block_index)
      call create_rad_array_mass_den(rho,radiation%rho,rad_block_index)
      call create_rad_array_scalar(temperature,radiation%temp,rad_block_index)
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
      call remap_to_original(radiation%rad_tend,radiation_tendency,rad_block_index)
      call remap_to_original_scalar_h(radiation%sfc_sw_in,sfc_sw_in,rad_block_index)
      call remap_to_original_scalar_h(radiation%sfc_lw_out,sfc_lw_out,rad_block_index)
      
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
  
  end subroutine call_radiation
  
  subroutine create_rad_array_scalar(in_array,out_array,rad_block_index) &
  bind(c,name = "create_rad_array_scalar")

    ! This subroutine cuts out a slice of a scalar field for hand-over to the radiation routine (done for RAM efficiency reasons).
    
    real(wp),intent(in)  :: in_array(n_scalars)
    real(wp),intent(out) :: out_array(n_scals_rad)
    integer              :: rad_block_index
    
    ! local variables
    integer :: ji,layer_index,h_index
    
    ! loop over all elements of the resulting array
    do ji=1,n_scals_rad
      layer_index = (ji-1)/n_scals_rad_h
      h_index = ji - layer_index*n_scals_rad_h
      out_array(ji) = in_array(rad_block_index*n_scals_rad_h + h_index + layer_index*n_scalars_h)
    enddo
  
  end subroutine create_rad_array_scalar
  
  subroutine create_rad_array_scalar_h(in_array,out_array,rad_block_index) &
  bind(c,name = "create_rad_array_scalar_h")

    ! This subroutine cuts out a slice of a horizontal scalar field for hand-over to the radiation routine (done for RAM efficiency reasons).
    
    real(wp),intent(in)  :: in_array(n_scalars_h)
    real(wp),intent(out) :: out_array(n_scals_rad_h)
    integer              :: rad_block_index
    
    ! local variables
    integer :: ji
    
    ! loop over all elements of the resulting array
    do ji=1,n_scals_rad_h
      out_array(ji) = in_array(rad_block_index*n_scals_rad_h + ji)
    enddo
  
  end subroutine create_rad_array_scalar_h

  subroutine create_rad_array_vector(in_array,out_array,rad_block_index) &
  bind(c,name = "create_rad_array_vector")

    ! This subroutine cuts out a slice of a vector field for hand-over to the radiation routine (done for RAM efficiency reasons).
  ! Only the vertical vector points are taken into account since only they are needed by the radiation.
    real(wp),intent(in)  :: in_array(n_vectors)
    real(wp),intent(out) :: out_array(n_scals_rad+n_scals_rad_h)
    integer              :: rad_block_index
    
    ! local variables
    integer :: ji,layer_index,h_index
    
    ! loop over all elements of the resulting array
    do ji=1,n_scals_rad+n_scals_rad_h
      layer_index = (ji-1)/n_scals_rad_h
      h_index = ji - layer_index*n_scals_rad_h
      out_array(ji) = in_array(rad_block_index*n_scals_rad_h + h_index + layer_index*n_vectors_per_layer)
    enddo
  
  end subroutine create_rad_array_vector
  
  subroutine create_rad_array_mass_den(in_array,out_array,rad_block_index) &
  bind(c,name = "create_rad_array_mass_den")

    ! This subroutine does same thing as create_rad_array_scalar,only for a mass density field.
    
    real(wp),intent(in)  :: in_array(n_constituents*n_scalars)
    real(wp),intent(out) :: out_array(n_constituents*n_scals_rad)
    integer              :: rad_block_index
    
    ! local variables
    integer :: const_id,ji,layer_index,h_index
    
    ! loop over all constituents
    do const_id=1,n_constituents
       ! loop over all elements of the resulting array
      do ji=1,n_scals_rad
        layer_index = (ji-1)/n_scals_rad_h
        h_index = ji - layer_index*n_scals_rad_h
        out_array((const_id-1)*n_scals_rad+ji) = in_array((const_id-1)*n_scalars+rad_block_index*n_scals_rad_h &
                                                          +h_index+layer_index*n_scalars_h)
      enddo
    enddo
  
  end subroutine create_rad_array_mass_den

  subroutine remap_to_original(in_array,out_array,rad_block_index) &
  bind(c,name = "remap_to_original")

    ! This subroutine reverses what create_rad_array_scalar has done.
    
    real(wp),intent(in)  :: in_array(n_scals_rad)
    real(wp),intent(out) :: out_array(n_scalars)
    integer              :: rad_block_index
    
    ! local variables
    integer :: ji,layer_index,h_index
    
    ! loop over all elements of the resulting array
    do ji=1,n_scals_rad
      layer_index = (ji-1)/n_scals_rad_h
      h_index = ji - layer_index*n_scals_rad_h
      out_array(rad_block_index*n_scals_rad_h+layer_index*n_scalars_h+h_index) = in_array(ji)
    enddo
  
  end subroutine remap_to_original

  subroutine remap_to_original_scalar_h(in_array,out_array,rad_block_index) &
  bind(c,name = "remap_to_original_scalar_h")

    ! This subroutine reverses what create_rad_array_scalar_h has done.
    
    real(wp),intent(in)  :: in_array(n_scals_rad_h)
    real(wp),intent(out) :: out_array(n_scalars_h)
    integer              :: rad_block_index
    
    ! local variables
    integer :: ji
    
    ! loop over all elements of the resulting array
    do ji=1,n_scals_rad_h
      out_array(rad_block_index*n_scals_rad_h+ji) = in_array(ji)
    enddo
  
  end subroutine remap_to_original_scalar_h

end module mo_manage_radiation_calls















