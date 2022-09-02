! This source file is part of the Geophysical Fluids Modeling Framework (GAME)which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_set_initial_state

  ! In this modulethe initial state of the simulation is set.

  use netcdf
  use mo_definitions,      only: wp
  use mo_constants,        only: p_0,r_d,c_d_p,m_d,m_v
  use mo_grid_nml,         only: n_scalars,n_scalars_h,n_vectors,n_vectors_per_layer,n_vectors_h,n_levels,n_dual_vectors, &
                                 n_layers,oro_id,res_id,n_dual_v_vectors,n_dual_scalars_h
  use mo_constituents_nml, only: n_condensed_constituents,n_constituents
  use mo_surface_nml,      only: nsoillays
  use mo_derived,          only: rel_humidity
  use mo_various_helpers,  only: nc_check,int2string
  use mo_grid_setup,       only: radius_rescale,z_t_const
  use mo_run_nml,          only: ideal_input_id
  use mo_rad_nml,          only: rad_config
  use baroclinic_wave,     only: baroclinic_wave_test
  use mo_multiplications,  only: scalar_times_vector
  use mo_vorticities,      only: calc_pot_vort
  use mo_vorticity_flux,   only: vorticity_flux
  use mo_inner_product,    only: inner_product
  
  implicit none
  
  contains
  
  subroutine set_ideal_init(exner_pert,theta_v_pert,scalar_field_placeholder,exner_bg,theta_v_bg,adjacent_vector_indices_h, &
                            area_dual,density_to_rhombi_indices,density_to_rhombi_weights,f_vec,flux_density, &
                            from_index,to_index,from_index_dual,to_index_dual,rho,inner_product_weights,normal_distance, &
                            pot_vort_tend,z_scalar,rhotheta_v,wind,v_squared,direction,latitude_scalar,longitude_scalar, &
                            z_vector,slope,gravity_potential,pot_vort,rel_vort,rel_vort_on_triangles,trsk_indices,trsk_weights, &
                            trsk_modified_curl_indices,z_vector_dual,vorticity_indices_triangles,vorticity_signs_triangles, &
                            t_const_soil,is_land,temperature_soil)
  
    ! This subroutine sets the initial state of the model atmosphere for idealized test cases.
  
    real(wp), intent(in)  :: exner_bg(n_scalars),theta_v_bg(n_scalars),area_dual(n_dual_vectors), &
                             density_to_rhombi_weights(4*n_vectors_h),f_vec(2*n_vectors_h),z_scalar(n_scalars), &
                             inner_product_weights(8*n_scalars),normal_distance(n_vectors), &
                             direction(n_vectors_h),latitude_scalar(n_scalars_h),longitude_scalar(n_scalars_h), &
                             z_vector(n_vectors),slope(n_vectors),gravity_potential(n_scalars),trsk_weights(10*n_vectors_h), &
                             z_vector_dual(n_dual_vectors),t_const_soil(n_scalars_h)
    real(wp), intent(out) :: exner_pert(n_scalars),theta_v_pert(n_scalars),scalar_field_placeholder(n_scalars), &
                             flux_density(n_vectors),rho(n_constituents*n_scalars),rhotheta_v(n_scalars),wind(n_vectors), &
                             v_squared(n_scalars),pot_vort((2*n_layers+1)*n_vectors_h),rel_vort((2*n_layers+1)*n_vectors_h), &
                             rel_vort_on_triangles(n_dual_v_vectors),temperature_soil(nsoillays*n_scalars_h), &
                             pot_vort_tend(n_vectors)
    integer,  intent(in)  :: adjacent_vector_indices_h(6*n_scalars_h),density_to_rhombi_indices(4*n_vectors_h), &
                             from_index(n_vectors_h),to_index(n_vectors_h),from_index_dual(n_vectors_h), &
                             to_index_dual(n_vectors_h),trsk_indices(10*n_vectors_h),trsk_modified_curl_indices(10*n_vectors_h), &
                             vorticity_indices_triangles(3*n_dual_scalars_h),vorticity_signs_triangles(3*n_dual_scalars_h), &
                             is_land(n_scalars_h)
    
    ! local variables
    integer               :: ji,jl,jc,layer_index,h_index,scalar_index,ncid_grid,latitude_vector_id, &
                             longitude_vector_id
    real(wp)              :: dummy_0,dummy_1,dummy_2,dummy_3,dummy_4,dummy_5,dummy_6,lat,lon,z_height,u,v, &
                             pressure_value,specific_humidity,dry_density,b,c,small_atmos_rescale
    real(wp), allocatable :: pressure(:),temperature(:),temperature_v(:),water_vapour_density(:),latitude_vector(:), &
                             longitude_vector(:)
    character(len=128)    :: grid_file_name
    
    ! determining the grid file
    grid_file_name = "../../grid_generator/grids/RES" // trim(int2string(res_id)) // "_L" // &
                     trim(int2string(n_layers)) // "_ORO" // trim(int2string(oro_id)) // ".nc"
    
    small_atmos_rescale = 1._wp/radius_rescale
    
    ! dummy arguments
    dummy_0 = 0._wp
    dummy_1 = 0._wp
    dummy_2 = 0._wp
    dummy_3 = 0._wp
    dummy_4 = 0._wp
    dummy_5 = 0._wp
    dummy_6 = 0._wp
    
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
        lat = latitude_scalar(h_index)
        lon = longitude_scalar(h_index)
        z_height = z_scalar(ji)
        ! standard atmosphere
        if (ideal_input_id==0) then
          temperature(ji) = theta_v_bg(ji)* exner_bg(ji)
          temperature_v(ji) = temperature(ji)
          pressure(ji) = p_0*exner_bg(ji)**(c_d_p/r_d)
        endif
        ! dry Ullrich test
        if (ideal_input_id==1) then
          call baroclinic_wave_test(1,0,1,small_atmos_rescale,lon,lat,pressure(ji),z_height,1,dummy_0,dummy_1,temperature(ji), &
                                    dummy_2,dummy_3,dummy_4,dummy_5,dummy_6)
          temperature_v(ji) = temperature(ji)
        endif
        ! moist Ullrich test
        if (ideal_input_id==2) then
          call baroclinic_wave_test(1,1,1,small_atmos_rescale,lon,lat,pressure(ji),z_height,1,dummy_0,dummy_1,temperature(ji), &
                                    dummy_2,dummy_4,dummy_5,dry_density,specific_humidity)
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
    call nc_check(nf90_inq_varid(ncid_grid,"latitude_vector",latitude_vector_id))
    call nc_check(nf90_inq_varid(ncid_grid,"longitude_vector",longitude_vector_id))
    call nc_check(nf90_get_var(ncid_grid,latitude_vector_id,latitude_vector))
    call nc_check(nf90_get_var(ncid_grid,longitude_vector_id,longitude_vector))
    call nc_check(nf90_close(ncid_grid))
    !$omp parallel do private(lat,lon,z_height,u,v,dummy_0,dummy_1,dummy_2,dummy_3,dummy_4,dummy_5,dummy_6)
    do ji=1,n_vectors_h
      do jl=0,n_layers-1
        lat = latitude_vector(ji)
        lon = longitude_vector(ji)
        z_height = z_vector(n_scalars_h + ji + jl*n_vectors_per_layer)
        ! standard atmosphere: no wind
        if (ideal_input_id==0) then
          wind(n_scalars_h + jl*n_vectors_per_layer + ji) = 0._wp          
                
          ! adding a "random" perturbation to the horizontal wind in the case of the Held-Suarez test case
          if (rad_config==2) then
            wind(n_scalars_h + jl*n_vectors_per_layer + ji) &
            = wind(n_scalars_h + jl*n_vectors_per_layer + ji) + 0.1_wp*mod(ji,17)/16._wp
          endif
        endif
        ! dry Ullrich test
        if (ideal_input_id==1) then
          call baroclinic_wave_test(1,0,1,small_atmos_rescale,lon,lat,dummy_0,z_height,1, &
                                    u,v,dummy_1,dummy_2,dummy_3,dummy_4,dummy_5,dummy_6)
          wind(n_scalars_h + jl*n_vectors_per_layer + ji) = u*cos(direction(ji)) + v*sin(direction(ji))
        endif
        ! moist Ullrich test
        if (ideal_input_id==2) then
          call baroclinic_wave_test(1,1,1,small_atmos_rescale,lon,lat,dummy_0,z_height,1, &
                                    u,v,dummy_1,dummy_2,dummy_3,dummy_4,dummy_5,dummy_6)
          wind(n_scalars_h + jl*n_vectors_per_layer + ji) = u*cos(direction(ji)) + v*sin(direction(ji))
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
        wind(jl*n_vectors_per_layer + ji) = 0._wp
      enddo
    enddo
    !$omp end parallel do
    
    ! this is the moist air density which has not yet been hydrostatically balanced
    !$omp parallel workshare
    scalar_field_placeholder = pressure/(r_d*temperature_v)
    !$omp end parallel workshare
    
    call scalar_times_vector(scalar_field_placeholder,wind,flux_density,from_index,to_index)
    ! Nowthe potential vorticity is evaluated.
    call calc_pot_vort(wind,rel_vort_on_triangles,z_vector,z_vector_dual,rel_vort, &
                       vorticity_indices_triangles,vorticity_signs_triangles,normal_distance, &
                       area_dual,from_index,to_index,from_index_dual,to_index_dual,inner_product_weights, &
                       slope,f_vec,pot_vort,density_to_rhombi_indices,density_to_rhombi_weights, &
                       scalar_field_placeholder)
    ! Nowthe generalized Coriolis term is evaluated.
    call vorticity_flux(from_index,to_index,pot_vort_tend,trsk_indices,trsk_modified_curl_indices,trsk_weights, &
                        flux_density,pot_vort,inner_product_weights,adjacent_vector_indices_h)
    
    ! Kinetic energy is prepared for the gradient term of the Lamb transformation.
    call inner_product(wind,wind,v_squared,adjacent_vector_indices_h,inner_product_weights)
    
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
          exner_pert(scalar_index) = (pressure_value/p_0)**(r_d/c_d_p)
        ! other layers
        else
          ! solving a quadratic equation for the Exner pressure
          b = -0.5_wp*exner_pert(scalar_index + n_scalars_h)/temperature_v(scalar_index + n_scalars_h) &
          *(temperature_v(scalar_index) - temperature_v(scalar_index + n_scalars_h) &
          + 2._wp/c_d_p*(gravity_potential(scalar_index) - gravity_potential(scalar_index + n_scalars_h) &
          + 0.5_wp*v_squared(scalar_index) - 0.5_wp*v_squared(scalar_index + n_scalars_h) &
          - (z_scalar(scalar_index) - z_scalar(scalar_index + n_scalars_h))*pot_vort_tend(ji + (jl+1)*n_vectors_per_layer)))
          c = exner_pert(scalar_index + n_scalars_h)**2*temperature_v(scalar_index)/temperature_v(scalar_index + n_scalars_h)
          exner_pert(scalar_index) = b + (b**2+c)**0.5_wp
        endif
        ! this is the full virtual potential temperature here
        theta_v_pert(scalar_index) = temperature_v(scalar_index)/exner_pert(scalar_index)
        
        ! scalar_field_placeholder is the moist air gas density here
        scalar_field_placeholder(scalar_index) = p_0*exner_pert(scalar_index)**(c_d_p/r_d)/(r_d*temperature_v(scalar_index))
        
        ! setting rhotheta_v according to its definition
        rhotheta_v(scalar_index) = scalar_field_placeholder(scalar_index)*theta_v_pert(scalar_index)
      enddo
    enddo
    !$omp end parallel do
    
    deallocate(pressure)
    deallocate(temperature_v)
    
    ! substracting the background state
    !$omp parallel workshare
    exner_pert = exner_pert - exner_bg
    theta_v_pert = theta_v_pert - theta_v_bg
    !$omp end parallel workshare
    
    !$omp parallel do private(ji,jc)
    do ji=1,n_scalars
      do jc=0,n_condensed_constituents-1
        ! condensed densities are zero in all test states
        rho(jc*n_scalars + ji) = 0._wp
      enddo
      ! the moist air density
      rho(n_condensed_constituents*n_scalars + ji) = scalar_field_placeholder(ji)
      ! water vapour density
      if (n_condensed_constituents==4) then
        rho((n_condensed_constituents+1)*n_scalars + ji) = water_vapour_density(ji)
      endif
    enddo
    !$omp end parallel do
    
    deallocate(water_vapour_density)
    
    ! setting the soil temperature
    call set_soil_temp(is_land,"NONE",t_const_soil,temperature,temperature_soil)
    deallocate(temperature)
    
  end subroutine set_ideal_init

  subroutine read_init_data(init_state_file,rho,wind,rhotheta_v,theta_v_pert,exner_pert,t_const_soil,tke, &
                            temperature_soil,is_land,theta_v_bg,exner_bg)
    
    ! This subroutine reads the initial state of the model atmosphere from a netCDF file.
    
    character(len=128), intent(in)  :: init_state_file
    real(wp),           intent(out) :: rho(n_constituents*n_scalars),wind(n_vectors),temperature_soil(nsoillays*n_scalars_h), &
                                       rhotheta_v(n_scalars),theta_v_pert(n_scalars),exner_pert(n_scalars),tke(n_scalars)
    real(wp),           intent(in)  :: t_const_soil(n_scalars_h),theta_v_bg(n_scalars),exner_bg(n_scalars)
    integer,            intent(in)  :: is_land(n_scalars_h)
    
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
          rho((n_condensed_constituents+1)*n_scalars + ji) = rho((n_condensed_constituents+1)*n_scalars+ji) &
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
    call set_soil_temp(is_land,init_state_file,t_const_soil,temperature,temperature_soil)
    
    deallocate(temperature)
    
  end subroutine read_init_data

  subroutine set_soil_temp(is_land,init_state_file,t_const_soil,temperature,temperature_soil)
    
    ! This subroutine sets the soil and SST temperature.
    
    integer,          intent(in)  :: is_land(n_scalars_h)
    character(len=*), intent(in)  :: init_state_file
    real(wp),         intent(in)  :: t_const_soil(n_scalars_h),temperature(n_scalars)
    real(wp),         intent(out) :: temperature_soil(nsoillays*n_scalars_h)
    
    ! local variables
    integer               :: soil_layer_index,ji,ncid,sst_id,soil_index,sst_avail,t_soil_avail,soil_id
    real(wp)              :: z_soil,t_sfc
    real(wp), allocatable :: sst(:)
    
    ! figuring out if the SST is included in the initialization file and reading it if it exists (important for NWP)
    allocate(sst(n_scalars_h))
    
    sst_avail = 0
    if (init_state_file/="NONE") then
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
    if (init_state_file/="NONE") then
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
      ! state filewe obtain it by linearly interpolating between the surface
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
  
end module mo_set_initial_state




