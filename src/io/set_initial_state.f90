! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module set_initial_state

  ! In this module, the initial state of the simulation is set.

  use iso_c_binding
  use netcdf
  use definitions,        only: wp
  use constants,          only: p_0,r_d,c_d_p,m_d,m_v
  use grid_nml,           only: n_scalars,n_scalars_h,n_vectors
  use constituents_nml,   only: n_condensed_constituents,n_constituents
  use surface_nml,        only: nsoillays
  use derived_quantities, only: rel_humidity
  use various_helpers,    only: nc_check
  
  implicit none
  
  contains

  subroutine set_soil_temp(is_land,init_state_file,t_const_soil,z_t_const,temperature,temperature_soil) &
  bind(c,name = "set_soil_temp")
    
    ! This subroutine sets the soil and SST temperature.
    
    integer,  intent(in)         :: is_land(n_scalars_h)
    character(len=1), intent(in) :: init_state_file
    real(wp), intent(in)         :: t_const_soil(n_scalars_h),temperature(n_scalars),z_t_const
    real(wp), intent(out)        :: temperature_soil(nsoillays*n_scalars_h)
    
    ! local variables
    integer               :: soil_layer_index,ji,ncid,sst_id,soil_index,sst_avail,t_soil_avail,soil_id
    real(wp)              :: z_soil,t_sfc
    real(wp), allocatable :: sst(:)
    
    ! figuring out if the SST is included in the initialization file and reading it if it exists (important for NWP)
    allocate(sst(n_scalars_h))
    
    sst_avail = 0
    if (len(init_state_file)/=0) then
      call nc_check(nf90_open(init_state_file,NF90_CLOBBER,ncid))
      
      ! figuring out if the netcdf file contains SST
      if (nf90_inq_varid(ncid,"sst",sst_id)==0) then
        sst_avail = 1
        write(*,*) "SST found in initialization file."
      else
        write(*,*) "SST not found in initialization file."
      endif
      
      ! reading the SST data if it is present in the netcdf file
      if (sst_avail==1) then
        call nc_check(nf90_get_var(ncid,sst_id,sst))
      endif
      
      ! we do not need the netcdf file any further
      call nc_check(nf90_close(ncid))
    endif
    
    ! figuring out if the soil temperature is included in the initialization file and reading it if it exists (important for NWP)
    t_soil_avail = 0
    if (len(init_state_file)/=0) then
      call nc_check(nf90_open(init_state_file,NF90_CLOBBER,ncid))
      
      ! figuring out if the netcdf file contains the soil temperature
      if (nf90_inq_varid(ncid,"t_soil",soil_id)==0) then
        t_soil_avail = 1
        write(*,*) "Soil temperature found in initialization file."
      else
        write(*,*) "Soil temperature not found in initialization file."
      endif
      
      ! reading the soil temperature if it is present in the netcdf file
      if (t_soil_avail==1) then
        call nc_check(nf90_get_var(ncid,soil_id,temperature_soil))
      endif
      
      ! we do not need the netcdf file any further
      call nc_check(nf90_close(ncid))
    endif
    
    ! setting what has not yet been set
    !$omp parallel do private(ji,soil_layer_index,soil_index,z_soil,t_sfc)
    do ji=1,n_scalars_h
      ! sea surface temperature if SST is available
      if (is_land(ji)==0 .and. sst_avail==1) then
        ! loop over all soil layers
        do soil_layer_index=0,nsoillays-1
          temperature_soil(ji+soil_layer_index*n_scalars_h) = sst(ji)
        enddo
      endif
      
      ! if the soil temperature over land or the SST over water is not available in the initialization
      ! state file, we obtain it by linearly interpolating between the surface
      ! and the depth of constant temperature    
      if ((is_land(ji)==1 .and. t_soil_avail==0) .or. (is_land(ji)==0 .and. sst_avail==0)) then
        ! setting the surface temperature identical to the air temperature in the lowest layer
        t_sfc = temperature(n_scalars-n_scalars_h+ji)
        
        ! loop over all soil layers
        do soil_layer_index=0,nsoillays-1
          ! index of this soil grid point
          soil_index = ji+soil_layer_index*n_scalars_h
          z_soil = z_t_const/nsoillays*(0.5_wp+soil_layer_index)
          temperature_soil(soil_index) = t_sfc + (t_const_soil(ji) - t_sfc)*z_soil/z_t_const
        enddo
      endif
    enddo
    !$omp end parallel do
    
    deallocate(sst)
  
  end subroutine set_soil_temp

  subroutine read_init_data(init_state_file,rho,wind,rhotheta_v,theta_v_pert,exner_pert,t_const_soil,tke, &
                            temperature_soil,is_land,z_t_const,theta_v_bg,exner_bg) &
  bind(c,name = "read_init_data")
    
    ! This subroutine reads the initial state of the model atmosphere from a netCDF file.
    
    character(len=1), intent(in)  :: init_state_file
    real(wp),         intent(out) :: rho(n_constituents*n_scalars),wind(n_vectors),temperature_soil(nsoillays*n_scalars_h), &
                                     rhotheta_v(n_scalars),theta_v_pert(n_scalars),exner_pert(n_scalars),tke(n_scalars)
    real(wp),         intent(in)  :: t_const_soil(n_scalars_h),z_t_const,theta_v_bg(n_scalars),exner_bg(n_scalars)
    integer,          intent(in)  :: is_land(n_scalars_h)
    
    ! local variables
    integer               :: ji,ncid,tke_id,tke_avail,densities_id,temperature_id,wind_id
    real(wp)              :: pressure,pot_temp_v
    real(wp), allocatable :: temperature(:),temperature_v(:)
    
    allocate(temperature(n_scalars))
    call nc_check(nf90_open(init_state_file,NF90_CLOBBER,ncid))
    call nc_check(nf90_inq_varid(ncid,"densities",densities_id))
    call nc_check(nf90_inq_varid(ncid,"temperature",temperature_id))
    call nc_check(nf90_inq_varid(ncid,"wind",wind_id))
    tke_avail = 0
    if (nf90_inq_varid(ncid,"tke",tke_id)==0) then
      tke_avail = 1
      write(*,*) "TKE found in initialization file."
    else
      write(*,*) "TKE not found in initialization file. TKE set to zero."
    endif
    call nc_check(nf90_get_var(ncid,densities_id,rho))
    call nc_check(nf90_get_var(ncid,temperature_id,temperature))
    call nc_check(nf90_get_var(ncid,wind_id,wind))
    if (tke_avail==1) then
      call nc_check(nf90_get_var(ncid,tke_id,tke))
    endif
    call nc_check(nf90_close(ncid))
    
    ! resricting the maximum relative humidity to 100 %
    if (n_condensed_constituents==4) then
      !$omp parallel do private(ji)
      do ji=1,n_scalars
        if (rel_humidity(rho((n_condensed_constituents+1)*n_scalars+ji),temperature(ji))>1._wp) then
          rho((n_condensed_constituents+1)*n_scalars + ji) = rho((n_condensed_constituents + 1)*n_scalars+ji) &
          /rel_humidity(rho((n_condensed_constituents+1)*n_scalars+ji),temperature(ji))
        endif
      enddo
      !$omp end parallel do
    endif
    
    ! diagnostic thermodynamical quantities
    allocate(temperature_v(n_scalars))
    !$omp parallel do private(ji,pressure,pot_temp_v)
    do ji=1,n_scalars
      temperature_v(ji) = temperature(ji) &
      *(1._wp+rho((n_condensed_constituents+1)*n_scalars+ji)/rho(n_condensed_constituents*n_scalars+ji)*(m_d/m_v-1._wp))
      pressure = rho(n_condensed_constituents*n_scalars+ji)*r_d*temperature_v(ji)
      pot_temp_v = temperature_v(ji)*(p_0/pressure)**(r_d/c_d_p)
      rhotheta_v(ji) = rho(n_condensed_constituents*n_scalars+ji)*pot_temp_v
      ! calculating the virtual potential temperature perturbation
      theta_v_pert(ji) = pot_temp_v - theta_v_bg(ji)
      ! calculating the Exner pressure perturbation
      exner_pert(ji) = temperature_v(ji)/(theta_v_bg(ji)+theta_v_pert(ji)) - exner_bg(ji)
    enddo
    !$omp end parallel do
    
    deallocate(temperature_v)
    
    ! checking for negative densities
    !$omp parallel do private(ji)
    do ji=1,n_constituents*n_scalars
      if (rho(ji)<0._wp) then
        write(*,*) "Negative density found."
        write(*,*) "Aborting."
        call exit(1)
      endif
    enddo
    !$omp end parallel do
    
    ! setting the soil temperature
    call set_soil_temp(is_land,init_state_file,t_const_soil,z_t_const,temperature,temperature_soil)
    
    deallocate(temperature)
    
  end subroutine read_init_data
  
end module 




