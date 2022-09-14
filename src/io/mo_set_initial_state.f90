! This source file is part of the Geophysical Fluids Modeling Framework (GAME)which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_set_initial_state

  ! In this modulethe initial state of the simulation is set.

  use netcdf
  use mo_definitions,      only: wp,t_grid,t_state,t_diag
  use mo_constants,        only: p_0,r_d,c_d_p,m_d,m_v
  use mo_grid_nml,         only: n_scalars,n_scalars_h,n_vectors,n_vectors_per_layer,n_vectors_h,n_levels,n_dual_vectors, &
                                 n_layers,oro_id,res_id,n_dual_v_vectors,n_dual_scalars_h
  use mo_constituents_nml, only: n_condensed_constituents,n_constituents,lmoist
  use mo_surface_nml,      only: nsoillays
  use mo_derived,          only: rel_humidity
  use mo_various_helpers,  only: nc_check,int2string
  use mo_grid_setup,       only: radius_rescale,z_t_const
  use mo_io_nml,           only: ideal_input_id,init_state_file
  use mo_rad_nml,          only: rad_config
  use baroclinic_wave,     only: baroclinic_wave_test
  use mo_multiplications,  only: scalar_times_vector
  use mo_vorticities,      only: calc_pot_vort
  use mo_vorticity_flux,   only: vorticity_flux
  use mo_inner_product,    only: inner_product
  
  implicit none
  
  contains
  
  subroutine set_ideal_init(state,diag,grid)
  
    ! This subroutine sets the initial state of the model atmosphere for idealized test cases.
    
    type(t_state), intent(inout) :: state
    type(t_diag),  intent(inout) :: diag
    type(t_grid),  intent(in)    :: grid
    
    ! local variables
    integer               :: ji,jl,jc,layer_index,h_index,scalar_index,ncid_grid,latitude_vector_id, &
                             longitude_vector_id
    real(wp)              :: dummy_1,dummy_2,dummy_3,dummy_4,dummy_5,dummy_6,dummy_7,lat,lon,z_height,u,v, &
                             pressure_value,specific_humidity,dry_density,b,c,small_atmos_rescale
    real(wp), allocatable :: pressure(:),temperature(:),temperature_v(:),water_vapour_density(:),latitude_vector(:), &
                             longitude_vector(:)
    character(len=128)    :: grid_file_name
    
    ! determining the grid file
    grid_file_name = "../../grid_generator/grids/RES" // trim(int2string(res_id)) // "_L" // &
                     trim(int2string(n_layers)) // "_ORO" // trim(int2string(oro_id)) // ".nc"
    
    small_atmos_rescale = 1._wp/radius_rescale
    
    ! dummy arguments
    dummy_1 = 0._wp
    dummy_2 = 0._wp
    dummy_3 = 0._wp
    dummy_4 = 0._wp
    dummy_5 = 0._wp
    dummy_6 = 0._wp
    dummy_7 = 0._wp
    
    allocate(pressure(n_scalars))
    allocate(temperature(n_scalars))
    allocate(temperature_v(n_scalars))
    allocate(water_vapour_density(n_scalars))
    
    !$omp parallel workshare
    pressure = 0._wp
    temperature = 0._wp
    temperature_v = 0._wp
    water_vapour_density = 0._wp
    !$omp end parallel workshare
    
    ! 3D scalar fields determined hereapart from density
    !$omp parallel do private(ji,layer_index,h_index,lat,lon,z_height,dry_density,specific_humidity)
    do ji=1,n_scalars
      layer_index = (ji-1)/n_scalars_h
      h_index = ji - layer_index*n_scalars_h
        lat = grid%lat_c(h_index)
        lon = grid%lon_c(h_index)
        z_height = grid%z_scalar(ji)
        ! standard atmosphere
        if (ideal_input_id==0) then
          temperature(ji) = grid%theta_v_bg(ji)*grid%exner_bg(ji)
          temperature_v(ji) = temperature(ji)
          pressure(ji) = p_0*grid%exner_bg(ji)**(c_d_p/r_d)
        endif
        ! dry Ullrich test
        if (ideal_input_id==1) then
          call baroclinic_wave_test(1,0,1,small_atmos_rescale,lon,lat,pressure(ji),z_height,1,dummy_1,dummy_2,temperature(ji), &
                                    dummy_3,dummy_4,dummy_5,dummy_6,dummy_7)
          temperature_v(ji) = temperature(ji)
        endif
        ! moist Ullrich test
        if (ideal_input_id==2) then
          call baroclinic_wave_test(1,1,1,small_atmos_rescale,lon,lat,pressure(ji),z_height,1,dummy_1,dummy_2,temperature(ji), &
                                    dummy_3,dummy_5,dummy_6,dry_density,specific_humidity)
          temperature_v(ji) = temperature(ji)*(1._wp+specific_humidity*(m_d/m_v-1._wp))
          water_vapour_density(ji) = dry_density*specific_humidity/(1._wp-specific_humidity)
        endif
    enddo
    !$omp end parallel do
    
    ! resricting the maximum relative humidity to 100 %
    if (n_condensed_constituents==4) then
      !$omp parallel do private(ji)
      do ji=1,n_scalars
        if (rel_humidity(water_vapour_density(ji),temperature(ji))>1._wp) then
          water_vapour_density(ji) = water_vapour_density(ji)/rel_humidity(water_vapour_density(ji),temperature(ji))
        endif
      enddo
      !$omp end parallel do
    endif

    ! horizontal wind fields are determind here
    ! reading the grid properties which are not part of the struct grid
    allocate(latitude_vector(n_vectors_h))
    allocate(longitude_vector(n_vectors_h))
    call nc_check(nf90_open(grid_file_name,NF90_CLOBBER,ncid_grid))
    call nc_check(nf90_inq_varid(ncid_grid,"lat_e",latitude_vector_id))
    call nc_check(nf90_inq_varid(ncid_grid,"lon_e",longitude_vector_id))
    call nc_check(nf90_get_var(ncid_grid,latitude_vector_id,latitude_vector))
    call nc_check(nf90_get_var(ncid_grid,longitude_vector_id,longitude_vector))
    call nc_check(nf90_close(ncid_grid))
    !$omp parallel do private(ji,jl,lat,lon,z_height,u,v,dummy_1,dummy_2,dummy_3,dummy_4,dummy_5,dummy_6,dummy_7)
    do ji=1,n_vectors_h
      do jl=0,n_layers-1
        lat = latitude_vector(ji)
        lon = longitude_vector(ji)
        z_height = grid%z_vector(n_scalars_h + ji + jl*n_vectors_per_layer)
        ! standard atmosphere: no wind
        if (ideal_input_id==0) then
          state%wind(n_scalars_h + jl*n_vectors_per_layer + ji) = 0._wp          
                
          ! adding a "random" perturbation to the horizontal wind in the case of the Held-Suarez test case
          if (rad_config==2) then
            state%wind(n_scalars_h + jl*n_vectors_per_layer + ji) &
            = state%wind(n_scalars_h + jl*n_vectors_per_layer + ji) + 0.1_wp*mod(ji,17)/16._wp
          endif
        endif
        ! dry Ullrich test
        if (ideal_input_id==1) then
          call baroclinic_wave_test(1,0,1,small_atmos_rescale,lon,lat,dummy_1,z_height,1, &
                                    u,v,dummy_2,dummy_3,dummy_4,dummy_5,dummy_6,dummy_7)
          state%wind(n_scalars_h + jl*n_vectors_per_layer + ji) = u*cos(grid%direction(ji)) + v*sin(grid%direction(ji))
        endif
        ! moist Ullrich test
        if (ideal_input_id==2) then
          call baroclinic_wave_test(1,1,1,small_atmos_rescale,lon,lat,dummy_1,z_height,1, &
                                    u,v,dummy_2,dummy_3,dummy_4,dummy_5,dummy_6,dummy_7)
          state%wind(n_scalars_h + jl*n_vectors_per_layer + ji) = u*cos(grid%direction(ji)) + v*sin(grid%direction(ji))
        endif
      enddo
    enddo
    !$omp end parallel do
    
    deallocate(latitude_vector)
    deallocate(longitude_vector)
    
    ! setting the vertical wind field equal to zero
    !$omp parallel do private(ji,jl)
    do ji=1,n_scalars_h
      do jl=0,n_levels-1
        state%wind(jl*n_vectors_per_layer + ji) = 0._wp
      enddo
    enddo
    !$omp end parallel do
    
    ! this is the moist air density which has not yet been hydrostatically balanced
    !$omp parallel workshare
    diag%scalar_field_placeholder = pressure/(r_d*temperature_v)
    !$omp end parallel workshare
    
    call scalar_times_vector(diag%scalar_field_placeholder,state%wind,diag%flux_density,grid)
    ! Nowthe potential vorticity is evaluated.
    call calc_pot_vort(state,diag%scalar_field_placeholder,diag,grid)
    ! Nowthe generalized Coriolis term is evaluated.
    call vorticity_flux(diag,grid)
    
    ! Kinetic energy is prepared for the gradient term of the Lamb transformation.
    call inner_product(state%wind,state%wind,diag%v_squared,grid)
    
    ! density is determined out of the hydrostatic equation
    ! theta_v_pert and exner_pert are a misuse of name herethey contain the full values here
    !$omp parallel do private(ji,scalar_index,b,c,pressure_value)
    do ji=1,n_scalars_h
      ! integrating from bottom to top
      do jl=n_layers-1,0,-1
        scalar_index = jl*n_scalars_h + ji
        ! lowest layer
        if (jl==n_layers-1) then
          pressure_value = pressure(scalar_index)
          state%exner_pert(scalar_index) = (pressure_value/p_0)**(r_d/c_d_p)
        ! other layers
        else
          ! solving a quadratic equation for the Exner pressure
          b = -0.5_wp*state%exner_pert(scalar_index + n_scalars_h)/temperature_v(scalar_index + n_scalars_h) &
          *(temperature_v(scalar_index) - temperature_v(scalar_index + n_scalars_h) &
          + 2._wp/c_d_p*(grid%gravity_potential(scalar_index) - grid%gravity_potential(scalar_index + n_scalars_h) &
          + 0.5_wp*diag%v_squared(scalar_index) - 0.5_wp*diag%v_squared(scalar_index + n_scalars_h) &
          - (grid%z_scalar(scalar_index) - grid%z_scalar(scalar_index + n_scalars_h)) &
          *diag%pot_vort_tend(ji + (jl+1)*n_vectors_per_layer)))
          c = state%exner_pert(scalar_index + n_scalars_h)**2*temperature_v(scalar_index)/temperature_v(scalar_index + n_scalars_h)
          state%exner_pert(scalar_index) = b + (b**2+c)**0.5_wp
        endif
        ! this is the full virtual potential temperature here
        state%theta_v_pert(scalar_index) = temperature_v(scalar_index)/state%exner_pert(scalar_index)
        
        ! scalar_field_placeholder is the moist air gas density here
        diag%scalar_field_placeholder(scalar_index) = p_0*state%exner_pert(scalar_index)**(c_d_p/r_d) &
                                                      /(r_d*temperature_v(scalar_index))
        
        ! setting rhotheta_v according to its definition
        state%rhotheta_v(scalar_index) = diag%scalar_field_placeholder(scalar_index)*state%theta_v_pert(scalar_index)
      enddo
    enddo
    !$omp end parallel do
    
    deallocate(pressure)
    deallocate(temperature_v)
    
    ! substracting the background state
    !$omp parallel workshare
    state%exner_pert = state%exner_pert - grid%exner_bg
    state%theta_v_pert = state%theta_v_pert - grid%theta_v_bg
    !$omp end parallel workshare
    
    !$omp parallel do private(ji,jc)
    do ji=1,n_scalars
      do jc=0,n_condensed_constituents-1
        ! condensed densities are zero in all test states
        state%rho(jc*n_scalars + ji) = 0._wp
      enddo
      ! the moist air density
      state%rho(n_condensed_constituents*n_scalars + ji) = diag%scalar_field_placeholder(ji)
      ! water vapour density
      if (n_condensed_constituents==4) then
        state%rho((n_condensed_constituents+1)*n_scalars + ji) = water_vapour_density(ji)
      endif
    enddo
    !$omp end parallel do
    
    deallocate(water_vapour_density)
    
    ! setting the soil temperature
    call set_soil_temp(state,temperature,grid)
    deallocate(temperature)
    
  end subroutine set_ideal_init

  subroutine read_init_data(state,diag,grid)
    
    ! This subroutine reads the initial state of the model atmosphere from a netCDF file.
    
    type(t_state),      intent(inout) :: state ! state variables
    type(t_diag),       intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),       intent(in)    :: grid  ! grid quantities
    
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
    call nc_check(nf90_get_var(ncid,densities_id,state%rho))
    call nc_check(nf90_get_var(ncid,temperature_id,temperature))
    call nc_check(nf90_get_var(ncid,wind_id,state%wind))
    if (tke_avail==1) then
      call nc_check(nf90_get_var(ncid,tke_id,diag%tke))
    else
      !$omp parallel workshare
      diag%tke = 0._wp
      !$omp end parallel workshare
    endif
    call nc_check(nf90_close(ncid))
    
    ! resricting the maximum relative humidity to 100 %
    if (lmoist) then
      !$omp parallel do private(ji)
      do ji=1,n_scalars
        if (rel_humidity(state%rho((n_condensed_constituents+1)*n_scalars+ji),temperature(ji))>1._wp) then
          state%rho((n_condensed_constituents+1)*n_scalars + ji) = state%rho((n_condensed_constituents+1)*n_scalars+ji) &
          /rel_humidity(state%rho((n_condensed_constituents+1)*n_scalars+ji),temperature(ji))
        endif
      enddo
      !$omp end parallel do
    endif
    
    ! diagnostic thermodynamical quantities
    allocate(temperature_v(n_scalars))
    !$omp parallel do private(ji,pressure,pot_temp_v)
    do ji=1,n_scalars
      temperature_v(ji) = temperature(ji) &
      *(1._wp+state%rho((n_condensed_constituents+1)*n_scalars+ji)/state%rho(n_condensed_constituents*n_scalars+ji)*(m_d/m_v-1._wp))
      pressure = state%rho(n_condensed_constituents*n_scalars+ji)*r_d*temperature_v(ji)
      pot_temp_v = temperature_v(ji)*(p_0/pressure)**(r_d/c_d_p)
      state%rhotheta_v(ji) = state%rho(n_condensed_constituents*n_scalars+ji)*pot_temp_v
      ! calculating the virtual potential temperature perturbation
      state%theta_v_pert(ji) = pot_temp_v - grid%theta_v_bg(ji)
      ! calculating the Exner pressure perturbation
      state%exner_pert(ji) = temperature_v(ji)/(grid%theta_v_bg(ji)+state%theta_v_pert(ji)) - grid%exner_bg(ji)
    enddo
    !$omp end parallel do
    
    deallocate(temperature_v)
    
    ! checking for negative densities
    !$omp parallel do private(ji)
    do ji=1,n_constituents*n_scalars
      if (state%rho(ji)<0._wp) then
        write(*,*) "Negative density found."
        write(*,*) "Aborting."
        call exit(1)
      endif
    enddo
    !$omp end parallel do
    
    ! setting the soil temperature
    call set_soil_temp(state,temperature,grid)
    
    deallocate(temperature)
    
  end subroutine read_init_data

  subroutine set_soil_temp(state,temperature,grid)
    
    ! This subroutine sets the soil and SST temperature.
    
    type(t_state),    intent(inout) :: state                    ! state to which to write
    real(wp),         intent(in)    :: temperature(n_scalars)   ! air temperature
    type(t_grid),     intent(in)    :: grid                     ! grid quantities
    
    ! local variables
    integer               :: soil_layer_index,ji,ncid,sst_id,soil_index,sst_avail,t_soil_avail,soil_id
    real(wp)              :: z_soil,t_sfc
    real(wp), allocatable :: sst(:)
    
    ! figuring out if the SST is included in the initialization file and reading it if it exists (important for NWP)
    allocate(sst(n_scalars_h))
    
    sst_avail = 0
    if (ideal_input_id==-1) then
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
    if (ideal_input_id==-1) then
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
        call nc_check(nf90_get_var(ncid,soil_id,state%temperature_soil))
      endif
      
      ! we do not need the netcdf file any further
      call nc_check(nf90_close(ncid))
    endif
    
    ! setting what has not yet been set
    !$omp parallel do private(ji,soil_layer_index,soil_index,z_soil,t_sfc)
    do ji=1,n_scalars_h
      ! sea surface temperature if SST is available
      if (grid%is_land(ji)==0 .and. sst_avail==1) then
        ! loop over all soil layers
        do soil_layer_index=0,nsoillays-1
          state%temperature_soil(ji+soil_layer_index*n_scalars_h) = sst(ji)
        enddo
      endif
      
      ! if the soil temperature over land or the SST over water is not available in the initialization
      ! state filewe obtain it by linearly interpolating between the surface
      ! and the depth of constant temperature    
      if ((grid%is_land(ji)==1 .and. t_soil_avail==0) .or. (grid%is_land(ji)==0 .and. sst_avail==0)) then
        ! setting the surface temperature identical to the air temperature in the lowest layer
        t_sfc = temperature(n_scalars-n_scalars_h+ji)
        
        ! loop over all soil layers
        do soil_layer_index=0,nsoillays-1
          ! index of this soil grid point
          soil_index = ji+soil_layer_index*n_scalars_h
          z_soil = z_t_const/nsoillays*(0.5_wp+soil_layer_index)
          state%temperature_soil(soil_index) = t_sfc + (grid%t_const_soil(ji) - t_sfc)*z_soil/z_t_const
        enddo
      endif
    enddo
    !$omp end parallel do
    
    deallocate(sst)
  
  end subroutine set_soil_temp
  
end module mo_set_initial_state




