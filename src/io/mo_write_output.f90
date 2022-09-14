! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_write_output

  ! In this module,the output is written to netCDF files and integrals are written to text files if configured that way.
  ! In addition to that,some postprocessing diagnostics are also calculated here.
  
  use netcdf
  use mo_definitions,            only: wp,t_grid,t_state,t_diag
  use mo_constants,              only: c_d_p,r_d,t_0,EPSILON_SECURITY,c_d_v,r_v,p_0,M_PI,gravity
  use mo_grid_nml,               only: n_scalars,n_scalars_h,n_vectors_per_layer,n_vectors_h,n_dual_vectors, &
                                       n_lat_io_points,n_lon_io_points,n_vectors,n_levels,n_layers,n_dual_scalars_h, &
                                       n_dual_v_vectors,n_latlon_io_points
  use mo_various_helpers,        only: nc_check,find_min_index,int2string
  use mo_geodesy,                only: passive_turn
  use mo_constituents_nml,       only: n_constituents,n_condensed_constituents,lmoist,rain_velocity,snow_velocity
  use mo_io_nml,                 only: n_pressure_levels,pressure_levels,time_to_next_analysis_min,lmodel_level_output, &
                                       lsurface_output,lpressure_level_output,n_output_steps_10m_wind,ideal_input_id
  use mo_dictionary,             only: saturation_pressure_over_ice,saturation_pressure_over_water
  use mo_surface_nml,            only: nsoillays,lprog_soil_temp,pbl_scheme,lsfc_sensible_heat_flux,lsfc_phase_trans
  use mo_derived,                only: rel_humidity,gas_constant_diagnostics,temperature_diagnostics
  use mo_run_nml,                only: run_id,start_date,start_hour
  use mo_vorticities,            only: calc_rel_vort,calc_pot_vort
  use mo_spatial_ops_for_output, only: inner_product_tangential,epv_diagnostics,edges_to_cells_lowest_layer, &
                                       calc_uv_at_edge,edges_to_cells,curl_field_to_cells
  use mo_divergences,            only: div_h
  use mo_inner_product,          only: inner_product
  
  implicit none
  
  contains

  subroutine interpolate_to_ll(in_field,out_field,grid)
  
    ! This subroutine interpolates a single-layer scalar field to a lat-lon grid.
    
    real(wp),     intent(in)  :: in_field(n_scalars_h)                      ! the horizontal scalar field to interpolate
    real(wp),     intent(out) :: out_field(n_lat_io_points,n_lon_io_points) ! the resulting 2D scalar field on a latlon grid
    type(t_grid), intent(in)  :: grid                                       ! grid quantities
    
    ! local variables
    integer :: ji,jk,jm
    
    ! loop over all output points
    !$omp parallel do private(ji,jk,jm)
    do ji=1,n_lat_io_points
      do jk=1,n_lon_io_points
        ! initializing the result with zero
        out_field(ji,jk) = 0._wp
        ! 1/r-average
        do jm=1,5
          if (in_field(grid%latlon_interpol_indices(5*(jk-1 + n_lon_io_points*(ji-1))+jm))/=9999) then
            out_field(ji,jk) = out_field(ji,jk) + grid%latlon_interpol_weights(5*(jk-1+n_lon_io_points*(ji-1))+jm) &
                               *in_field(grid%latlon_interpol_indices(5*(jk-1+n_lon_io_points*(ji-1))+jm))
          else
            out_field(ji,jk) = 9999
            exit
          endif
        enddo
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine interpolate_to_ll

  subroutine write_out_integrals(state,diag,grid,time_since_init)
    
    ! This subroutine writes out fundamental integral properties of the atmosphere to a text file.
   
    type(t_grid),  intent(in) :: grid            ! grid properties
    type(t_diag),  intent(in) :: diag            ! diagnostic quantities
    type(t_state), intent(in) :: state           ! the state to use for writing the integrals
    real(wp),      intent(in) :: time_since_init ! the time since model initialization
    
    ! local variables
    integer               :: ji,const_id
    real(wp)              :: global_integral,kinetic_integral,potential_integral,internal_integral
    real(wp), allocatable :: int_energy_density(:),pot_energy_density(:),e_kin_density(:)
    
    global_integral = 0._wp
    
    ! masses
    if (time_since_init==0._wp) then
      open(1,file="masses",action="write")
    else
      open(1,file="masses",status="old",position="append",action="write")
    endif
    write(1,fmt="(F20.3)",advance="no") time_since_init
    do const_id=1,n_constituents
      global_integral = 0._wp
      do ji=1,n_scalars
        global_integral = global_integral + state%rho((const_id-1)*n_scalars + ji)*grid%volume(ji)
      enddo
      if (const_id==n_constituents) then
        write(1,fmt="(F30.3)") global_integral
      else
        write(1,fmt="(F30.3)",advance="no") global_integral
      endif
    enddo
    close(1)
        
    ! density times virtual potential temperature
    if (time_since_init==0._wp) then
      open(1,file="potential_temperature_density",action="write")
    else
      open(1,file="potential_temperature_density",status="old",position="append",action="write")
    endif
    global_integral = 0._wp
    do ji=1,n_scalars
      global_integral = global_integral + state%rhotheta_v(ji)*grid%volume(ji)
    enddo
    write(1,fmt="(F20.3,F30.3)") time_since_init,global_integral
    close(1)
        
    ! energies
    if (time_since_init==0._wp) then
      open(1,file="energy",action="write")
    else
      open(1,file="energy",status="old",position="append",action="write")
    endif
    allocate(e_kin_density(n_scalars))
    call inner_product(state%wind,state%wind,e_kin_density,grid)
    !$omp parallel do private(ji)
    do ji=1,n_scalars
      e_kin_density(ji) = state%rho(n_condensed_constituents*n_scalars+ji)*e_kin_density(ji)
    enddo
    !$omp end parallel do
    kinetic_integral = 0._wp
    do ji=1,n_scalars
      kinetic_integral = kinetic_integral + e_kin_density(ji)*grid%volume(ji)
    enddo
    deallocate(e_kin_density)
    allocate(pot_energy_density(n_scalars))
    !$omp parallel do private(ji)
    do ji=1,n_scalars
      pot_energy_density(ji) = state%rho(n_condensed_constituents*n_scalars + ji)*grid%gravity_potential(ji)
    enddo
    !$omp end parallel do
    potential_integral = 0._wp
    do ji=1,n_scalars
      potential_integral = potential_integral + pot_energy_density(ji)*grid%volume(ji)
     enddo
    deallocate(pot_energy_density)
    allocate(int_energy_density(n_scalars))
    !$omp parallel do private(ji)
    do ji=1,n_scalars
      int_energy_density(ji) = state%rho(n_condensed_constituents*n_scalars + ji)*diag%temperature(ji)
    enddo
    !$omp end parallel do
    internal_integral = 0._wp
    do ji=1,n_scalars
      internal_integral = internal_integral+int_energy_density(ji)*grid%volume(ji)
    enddo
    write(1,fmt="(F20.3,F30.3,F30.3,F30.3)") time_since_init,0.5_wp*kinetic_integral,potential_integral, &
                                             c_d_v*internal_integral
    deallocate(int_energy_density)
    close(1)
    
  end subroutine write_out_integrals

  subroutine write_out(state,diag,grid,wind_h_lowest_layer_array,t_init,t_write,ltotally_first_step)
  
    ! This subroutine is the central subroutine for writing the output.
  
    type(t_state), intent(in)    :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    real(wp), intent(in)  :: t_init,t_write, &
                             wind_h_lowest_layer_array(n_output_steps_10m_wind*n_vectors_h)
    logical,       intent(in)    :: ltotally_first_step
  
    ! local variables
    logical               :: lcontains_nan
    integer               :: ji,jk,jl,jm,lat_lon_dimids(2),ncid,single_int_dimid,lat_dimid,lon_dimid,start_date_id,start_hour_id, &
                             lat_id,lon_id,layer_index,closest_index,second_closest_index,temperature_ids(n_layers), &
                             pressure_ids(n_layers),rel_hum_ids(n_layers),wind_u_ids(n_layers),wind_v_ids(n_layers), &
                             rel_vort_ids(n_layers),div_h_ids(n_layers),wind_w_ids(n_levels), &
                             time_since_init_min,mslp_id,sp_id,rprate_id,sprate_id, &
                             cape_id,tcc_id,t2_id,u10_id,v10_id,gusts_id,sfc_sw_down_id,gh_ids(n_pressure_levels), &
                             temp_p_ids(n_pressure_levels),rh_p_ids(n_pressure_levels),wind_u_p_ids(n_pressure_levels), &
                             wind_v_p_ids(n_pressure_levels),epv_p_ids(n_pressure_levels),rel_vort_p_ids(n_pressure_levels), &
                             scalar_dimid,soil_dimid,vector_dimid,densities_dimid,densities_id,temperature_id,wind_id, &
                             tke_id,soil_id,time_step_10_m_wind,pressure_level_hpa
    real(wp)              :: delta_latitude,delta_longitude,lat_vector(n_lat_io_points),lon_vector(n_lon_io_points), &
                             min_precip_rate_mmh,min_precip_rate,cloud_water2cloudiness,temp_lowest_layer, &
                             pressure_value,mslp_factor,sp_factor,temp_mslp,temp_surface,z_height,theta_v, &
                             cape_integrand,delta_z,temp_closest,temp_second_closest,delta_z_temp,temperature_gradient, &
                             theta_e,u_850_surrogate,u_950_surrogate,u_850_proxy_height,u_950_proxy_height, &
                             wind_tangential,wind_u_value,wind_v_value, &
                             roughness_length_extrapolation,actual_roughness_length,z_sfc,z_agl,rescale_factor, &
                             cloud_water_content,vector_to_minimize(n_layers),closest_weight,z_tropopause, &
                             standard_vert_lapse_rate
    real(wp), allocatable :: wind_10_m_mean_u(:),wind_10_m_mean_v(:),mslp(:),sp(:),t2(:),tcc(:),rprate(:),sprate(:),cape(:), &
                             sfc_sw_down(:),geopotential_height(:,:),t_on_p_levels(:,:),rh_on_p_levels(:,:), &
                             epv_on_p_levels(:,:),u_on_p_levels(:,:),v_on_p_levels(:,:),zeta_on_p_levels(:,:), &
                             wind_10_m_mean_u_at_cell(:),wind_10_m_mean_v_at_cell(:),wind_10_m_gusts_speed_at_cell(:), &
                             div_h_all_layers(:),rel_vort_scalar_field(:),rh(:),epv(:),pressure(:),lat_lon_output_field(:,:), &
                             u_at_cell(:),v_at_cell(:),u_at_edge(:),v_at_edge(:)
    character(len=64)     :: output_file,output_file_p_level,varname
  
    write(*,*) "Writing output ..."
    
    ! checking for nan values
    !$omp parallel workshare
    lcontains_nan = any(isnan(state%exner_pert))
    !$omp end parallel workshare
    if (lcontains_nan) then
      write(*,*) "Congratulations, the model crashed."
      call exit(1)
    endif
    
    ! latitude resolution of the grid
    delta_latitude = M_PI/n_lat_io_points
    ! longitude resolution of the grid
    delta_longitude = 2._wp*M_PI/n_lon_io_points
    
    do ji=1,n_lat_io_points
      lat_vector(ji) = M_PI/2._wp - 0.5_wp*delta_latitude - (ji-1)*delta_latitude
    enddo
    do ji=1,n_lon_io_points
      lon_vector(ji) = (ji-1)*delta_longitude
    enddo
    
    ! precipitation rates smaller than this value are set to zero to not confuse users
    min_precip_rate_mmh = 0.01_wp
    min_precip_rate = min_precip_rate_mmh/(1000._wp*3600._wp/1024._wp)
    ! this heuristic coefficient converts the cloud water content to cloud cover
    cloud_water2cloudiness = 10._wp
    
    allocate(lat_lon_output_field(n_lat_io_points,n_lon_io_points))
    
    ! diagnosing the temperature
    call temperature_diagnostics(state,diag,grid)
    time_since_init_min = int(t_write - t_init)
    time_since_init_min = time_since_init_min/60
    
    ! Surface output including diagnostics.
    ! -------------------------------------
    if (lsurface_output) then
    
      allocate(mslp(n_scalars_h))
      allocate(sp(n_scalars_h))
      allocate(t2(n_scalars_h))
      allocate(tcc(n_scalars_h))
      allocate(rprate(n_scalars_h))
      allocate(sprate(n_scalars_h))
      allocate(cape(n_scalars_h))
      allocate(sfc_sw_down(n_scalars_h))
      z_tropopause = 12e3_wp
      standard_vert_lapse_rate = 0.0065_wp
      !$omp parallel do private(ji,jl,temp_lowest_layer,pressure_value,mslp_factor,sp_factor,temp_mslp,temp_surface, &
      !$omp z_height,theta_v,cape_integrand,delta_z,temp_closest,temp_second_closest,delta_z_temp,temperature_gradient, &
      !$omp theta_e,layer_index,closest_index,second_closest_index,cloud_water_content,vector_to_minimize)
      do ji=1,n_scalars_h
        ! Now the aim is to determine the value of the mslp.
        temp_lowest_layer = diag%temperature((n_layers-1)*n_scalars_h+ji)
        pressure_value = state%rho(n_condensed_constituents*n_scalars + (n_layers-1)*n_scalars_h + ji) &
        *gas_constant_diagnostics(state%rho,(n_layers-1)*n_scalars_h+ji-1)*temp_lowest_layer
        temp_mslp = temp_lowest_layer + standard_vert_lapse_rate*grid%z_scalar((n_layers-1)*n_scalars_h + ji)
        mslp_factor = (1._wp - (temp_mslp - temp_lowest_layer)/temp_mslp)**(grid%gravity_m((n_layers-1)*n_vectors_per_layer + ji)/ &
        (gas_constant_diagnostics(state%rho,(n_layers-1)*n_scalars_h+ji-1)*standard_vert_lapse_rate))
        mslp(ji) = pressure_value/mslp_factor
        
        ! Now the aim is to determine the value of the surface pressure.
        temp_surface = temp_lowest_layer + standard_vert_lapse_rate &
        *(grid%z_scalar(ji + (n_layers-1)*n_scalars_h) - grid%z_vector(n_vectors - n_scalars_h+ji))
        sp_factor = (1._wp - (temp_surface - temp_lowest_layer)/temp_surface) &
                     **(grid%gravity_m((n_layers-1)*n_vectors_per_layer + ji)/ &
                    (gas_constant_diagnostics(state%rho,(n_layers-1)*n_scalars_h+ji-1)*standard_vert_lapse_rate))
        sp(ji) = pressure_value/sp_factor
        
        ! Now the aim is to calculate the 2 m temperature.
        do jl=1,n_layers
          vector_to_minimize(jl) = abs(grid%z_vector(n_layers*n_vectors_per_layer+ji)+2._wp - grid%z_scalar(ji+(jl-1)*n_scalars_h))
        enddo
        closest_index = find_min_index(vector_to_minimize,n_layers)
        temp_closest = diag%temperature((closest_index-1)*n_scalars_h+ji)
        delta_z_temp = grid%z_vector(n_layers*n_vectors_per_layer+ji)+2._wp - grid%z_scalar(ji + (closest_index-1)*n_scalars_h)
        ! real radiation
        if (lprog_soil_temp) then
          temperature_gradient = (temp_closest - state%temperature_soil(ji))/ &
                                 (grid%z_scalar(ji+(closest_index-1)*n_scalars_h) - grid%z_vector(n_layers*n_vectors_per_layer+ji))
        ! no real radiation
        else
          second_closest_index = closest_index-1
          if (grid%z_scalar(ji + (closest_index-1)*n_scalars_h)>grid%z_vector(n_layers*n_vectors_per_layer + ji)+2._wp &
              .and. closest_index<n_layers) then
            second_closest_index = closest_index+1
          endif
          temp_second_closest = diag%temperature((second_closest_index-1)*n_scalars_h+ji)
          ! calculating the vertical temperature gradient that will be used for the extrapolation
          temperature_gradient = (temp_closest - temp_second_closest) &
                                 /(grid%z_scalar(ji+(closest_index-1)*n_scalars_h) &
                                 - grid%z_scalar(ji+(second_closest_index-1)*n_scalars_h))
        endif
        ! performing the interpolation / extrapolation to two meters above the surface
        t2(ji) = temp_closest + delta_z_temp*temperature_gradient
        
        ! diagnozing CAPE
        ! initializing CAPE with zero
        cape(ji) = 0._wp
        layer_index = n_layers - 1
        z_height = grid%z_scalar(layer_index*n_scalars_h+ji)
        ! pseduovirtual potential temperature of the particle in the lowest layer
        theta_e = pseudopotential_temperature(state,diag,(layer_index-1)*n_scalars_h+ji,grid)
        do while (z_height<z_tropopause)
          ! full virtual potential temperature in the grid box
          theta_v = grid%theta_v_bg(layer_index*n_scalars_h + ji) + state%theta_v_pert(layer_index*n_scalars_h + ji)
          ! thickness of the gridbox
          delta_z = grid%layer_thickness(layer_index*n_scalars_h + ji)
          ! this is the candidate that we might want to add to the integral
          cape_integrand = grid%gravity_m(layer_index*n_vectors_per_layer+ji)*(theta_e - theta_v)/theta_v
          ! we do not add negative values to CAPE (see the definition of CAPE)
          if (cape_integrand>0._wp) then
            cape(ji) = cape(ji) + cape_integrand*delta_z
          endif
          layer_index = layer_index-1
          z_height = grid%z_scalar(layer_index*n_scalars_h + ji)
        enddo
        
        sfc_sw_down(ji) = diag%sfc_sw_in(ji)/(1._wp-grid%sfc_albedo(ji)+EPSILON_SECURITY)
        
        ! Now come the hydrometeors.
        ! Calculation of the total cloud cover
        if (n_condensed_constituents==4) then
          ! calculating the cloud water content in this column
          cloud_water_content = 0._wp
          do jl=1,n_layers
            if (grid%z_scalar((jl-1)*n_scalars_h+ji)<z_tropopause) then
              cloud_water_content = cloud_water_content &
              + (state%rho(2*n_scalars + (jl-1)*n_scalars_h + ji) + state%rho(3*n_scalars + (jl-1)*n_scalars_h+ji)) &
              *(grid%z_vector(ji + (jl-1)*n_vectors_per_layer) - grid%z_vector(ji + jl*n_vectors_per_layer))
            endif
          enddo
          ! some heuristic ansatz for the total cloud cover
          tcc(ji) = min(cloud_water2cloudiness*cloud_water_content,1._wp)
          ! conversion of the total cloud cover into a percentage
          tcc(ji) = 100._wp*tcc(ji)
          ! setting too small values to zero to not confuse users
          if (tcc(ji)<0._wp) then
            tcc(ji) = 0._wp
          else
            tcc(ji) = 0._wp
          endif
          ! solid precipitation rate
          sprate(ji) = 0._wp
          if (n_condensed_constituents==4) then
            sprate(ji) = snow_velocity*state%rho((n_layers-1)*n_scalars_h+ji)
          endif
          ! liquid precipitation rate
          rprate(ji) = 0._wp
          if (n_condensed_constituents==4) then
            rprate(ji) = rain_velocity*state%rho(n_scalars + (n_layers-1)*n_scalars_h+ji)
          endif
          ! setting very small values to zero
          if (rprate(ji)<min_precip_rate) then
            rprate(ji) = 0._wp
          endif
          ! setting very small values to zero
          if (sprate(ji)<min_precip_rate) then
            sprate(ji) = 0._wp
          endif
        endif
      enddo
      !$omp end parallel do
      
      ! 10 m wind diagnostics
      ! ---------------------
      
      allocate(wind_10_m_mean_u(n_vectors_h))
      allocate(wind_10_m_mean_v(n_vectors_h))
      ! temporal average over the ten minutes output interval
      !$omp parallel do private(ji,jk,time_step_10_m_wind,wind_tangential,wind_u_value,wind_v_value)
      do ji=1,n_vectors_h
        ! initializing the means with zero
        wind_10_m_mean_u(ji) = 0._wp
        wind_10_m_mean_v(ji) = 0._wp
        ! loop over the time steps
        do time_step_10_m_wind=1,n_output_steps_10m_wind
          jk = (time_step_10_m_wind-1)*n_vectors_h + ji
          wind_tangential = 0._wp
          do jk=1,10
            wind_tangential = wind_tangential + grid%trsk_weights(10*(ji-1)+jk) &
                              *wind_h_lowest_layer_array((time_step_10_m_wind-1)*n_vectors_h+1+grid%trsk_indices(10*(ji-1)+jk))
          enddo
          wind_10_m_mean_u(ji) = wind_10_m_mean_u(ji) + 1._wp/n_output_steps_10m_wind*wind_h_lowest_layer_array(jk)
          wind_10_m_mean_v(ji) = wind_10_m_mean_v(ji) + 1._wp/n_output_steps_10m_wind*wind_tangential
        enddo
        ! passive turn to obtain the u- and v-components of the wind
        call passive_turn(wind_10_m_mean_u(ji),wind_10_m_mean_v(ji),-grid%direction(ji),wind_u_value,wind_v_value)
        wind_10_m_mean_u(ji) = wind_u_value
        wind_10_m_mean_v(ji) = wind_v_value
      enddo
      !$omp end parallel do
      ! vertically extrapolating to ten meters above the surface
      !$omp parallel do private(ji,roughness_length_extrapolation,actual_roughness_length,z_sfc,z_agl,rescale_factor)
      do ji=1,n_vectors_h
        actual_roughness_length = 0.5_wp*(grid%roughness_length(1+grid%from_index(ji)) + grid%roughness_length(1+grid%to_index(ji)))
        ! roughness length of grass according to WMO
        roughness_length_extrapolation = 0.02_wp
        if (grid%is_land(grid%from_index(ji)+1)==0) then
          roughness_length_extrapolation = actual_roughness_length
        endif
        z_sfc = 0.5_wp*(grid%z_vector(n_vectors - n_scalars_h + grid%from_index(ji)) &
                        + grid%z_vector(n_vectors - n_scalars_h + grid%to_index(ji)))
        z_agl = grid%z_vector(n_vectors - n_vectors_per_layer+ji) - z_sfc
        
        ! rescale factor for computing the wind in a height of 10 m
        rescale_factor = log(10._wp/roughness_length_extrapolation)/log(z_agl/actual_roughness_length)
        
        wind_10_m_mean_u(ji) = rescale_factor*wind_10_m_mean_u(ji)
        wind_10_m_mean_v(ji) = rescale_factor*wind_10_m_mean_v(ji)
      enddo
      !$omp end parallel do
      
      ! averaging the wind quantities to cell centers for output
      allocate(wind_10_m_mean_u_at_cell(n_scalars_h))
      call edges_to_cells_lowest_layer(wind_10_m_mean_u,wind_10_m_mean_u_at_cell,grid)
      deallocate(wind_10_m_mean_u)
      allocate(wind_10_m_mean_v_at_cell(n_scalars_h))
      call edges_to_cells_lowest_layer(wind_10_m_mean_v,wind_10_m_mean_v_at_cell,grid)
      deallocate(wind_10_m_mean_v)
      
      ! gust diagnostics
      u_850_proxy_height = 8000._wp*log(1000._wp/850._wp)
      u_950_proxy_height = 8000._wp*log(1000._wp/950._wp)
      allocate(wind_10_m_gusts_speed_at_cell(n_scalars_h))
      !$omp parallel do private(ji,jl,vector_to_minimize,closest_index,second_closest_index,u_850_surrogate,u_950_surrogate)
      do ji=1,n_scalars_h
      
        ! This is the normal case.
        if ((lsfc_sensible_heat_flux .or. lsfc_phase_trans .or. pbl_scheme==1) &
            .and. abs(diag%monin_obukhov_length(ji))>EPSILON_SECURITY) then
          ! This follows IFS DOCUMENTATION – Cy43r1 - Operational implementation 22 Nov 2016 - PART IV: PHYSICAL PROCESSES.
          wind_10_m_gusts_speed_at_cell(ji) = sqrt(wind_10_m_mean_u_at_cell(ji)**2 + wind_10_m_mean_v_at_cell(ji)**2) &
          + 7.71_wp*diag%roughness_velocity(ji)* &
          (max(1._wp - 0.5_wp/12._wp*1000._wp/diag%monin_obukhov_length(ji),0._wp))**(1._wp/3._wp)
          ! calculating the wind speed in a height representing 850 hPa
          do jl=1,n_layers
            vector_to_minimize(jl) = abs(grid%z_scalar((jl-1)*n_scalars_h+ji) &
                                         - (grid%z_vector(n_vectors-n_scalars_h+ji) + u_850_proxy_height))
          enddo
          closest_index = find_min_index(vector_to_minimize,n_layers)
          second_closest_index = closest_index-1
          if (closest_index<n_layers-1 &
              .and. grid%z_scalar((closest_index-1)*n_scalars_h+ji)-grid%z_vector(n_vectors-n_scalars_h+ji)>u_850_proxy_height) then
            second_closest_index = closest_index+1
          endif
          u_850_surrogate = sqrt(diag%v_squared(ji + (closest_index-1)*n_scalars_h)) &
          + (sqrt(diag%v_squared(ji+(closest_index-1)*n_scalars_h))-sqrt(diag%v_squared(ji+(second_closest_index-1)*n_scalars_h))) &
          /(grid%z_scalar(ji + (closest_index-1)*n_scalars_h) - grid%z_scalar(ji + (second_closest_index-1)*n_scalars_h)) &
          *(grid%z_vector(n_vectors - n_scalars_h+ji) + u_850_proxy_height - grid%z_scalar(ji + (closest_index-1)*n_scalars_h))
          ! calculating the wind speed in a height representing 950 hPa
          do jl=1,n_layers
            vector_to_minimize(jl) = abs(grid%z_scalar((jl-1)*n_scalars_h + ji) &
                                         - (grid%z_vector(n_vectors - n_scalars_h + ji) + u_950_proxy_height))
          enddo
          closest_index = find_min_index(vector_to_minimize,n_layers)
          second_closest_index = closest_index-1
          if (closest_index<n_layers &
              .and. grid%z_scalar((closest_index-1)*n_scalars_h+ji)-grid%z_vector(n_vectors-n_scalars_h+ji)>u_950_proxy_height) then
            second_closest_index = closest_index+1
          endif
          u_950_surrogate = sqrt(diag%v_squared(ji + (closest_index-1)*n_scalars_h)) &
          + (sqrt(diag%v_squared(ji+(closest_index-1)*n_scalars_h))-sqrt(diag%v_squared(ji+(second_closest_index-1)*n_scalars_h))) &
          /(grid%z_scalar(ji + (closest_index-1)*n_scalars_h) - grid%z_scalar(ji + (second_closest_index-1)*n_scalars_h)) &
          *(grid%z_vector(n_vectors - n_scalars_h+ji) + u_950_proxy_height - grid%z_scalar(ji + (closest_index-1)*n_scalars_h))
          ! adding the baroclinic and convective component to the gusts
          wind_10_m_gusts_speed_at_cell(ji) = wind_10_m_gusts_speed_at_cell(ji) &
                                              + 0.6_wp*max(0._wp,u_850_surrogate - u_950_surrogate)
          wind_10_m_gusts_speed_at_cell(ji) = min(wind_10_m_gusts_speed_at_cell(ji), &
                                                  3._wp*sqrt(wind_10_m_mean_u_at_cell(ji)**2 + wind_10_m_mean_v_at_cell(ji)**2))
        ! This is used if the turbulence quantities are not populated.
        else
          wind_10_m_gusts_speed_at_cell(ji) = 1.67_wp*sqrt(wind_10_m_mean_u_at_cell(ji)**2 + wind_10_m_mean_v_at_cell(ji)**2)
        endif

      enddo
      !$omp end parallel do
      
      output_file = trim(run_id) // "+" // trim(int2string(time_since_init_min)) // "min_surface.nc"
      
      call nc_check(nf90_create(output_file,NF90_CLOBBER,ncid))
      call nc_check(nf90_def_dim(ncid,"single_int_index",1,single_int_dimid))
      call nc_check(nf90_def_dim(ncid,"lat_index",n_lat_io_points,lat_dimid))
      call nc_check(nf90_def_dim(ncid,"lon_index",n_lon_io_points,lon_dimid))
        
      lat_lon_dimids(1) = lat_dimid
      lat_lon_dimids(2) = lon_dimid
      
      ! defining the variables
      call nc_check(nf90_def_var(ncid,"start_date",NF90_INT,single_int_dimid,start_date_id))
      call nc_check(nf90_def_var(ncid,"start_hour",NF90_INT,single_int_dimid,start_hour_id))
      call nc_check(nf90_def_var(ncid,"lat",NF90_REAL,lat_dimid,lat_id))
      call nc_check(nf90_def_var(ncid,"lon",NF90_REAL,lon_dimid,lon_id))
      call nc_check(nf90_def_var(ncid,"mslp",NF90_REAL,lat_lon_dimids,mslp_id))
      call nc_check(nf90_put_att(ncid,mslp_id,"units","Pa"))
      call nc_check(nf90_def_var(ncid,"sp",NF90_REAL,lat_lon_dimids,sp_id))
      call nc_check(nf90_put_att(ncid,sp_id,"units","Pa"))
      call nc_check(nf90_def_var(ncid,"t2",NF90_REAL,lat_lon_dimids,t2_id))
      call nc_check(nf90_put_att(ncid,t2_id,"units","K"))
      call nc_check(nf90_def_var(ncid,"tcc",NF90_REAL,lat_lon_dimids,tcc_id))
      call nc_check(nf90_put_att(ncid,tcc_id,"units","%"))
      call nc_check(nf90_def_var(ncid,"rprate",NF90_REAL,lat_lon_dimids,rprate_id))
      call nc_check(nf90_put_att(ncid,rprate_id,"units","kg/(m^2s)"))
      call nc_check(nf90_def_var(ncid,"sprate",NF90_REAL,lat_lon_dimids,sprate_id))
      call nc_check(nf90_put_att(ncid,sprate_id,"units","kg/(m^2s)"))
      call nc_check(nf90_def_var(ncid,"cape",NF90_REAL,lat_lon_dimids,cape_id))
      call nc_check(nf90_put_att(ncid,cape_id,"units","J/kg"))
      call nc_check(nf90_def_var(ncid,"sfc_sw_down",NF90_REAL,lat_lon_dimids,sfc_sw_down_id))
      call nc_check(nf90_put_att(ncid,sfc_sw_down_id,"units","W/m^2"))
      call nc_check(nf90_def_var(ncid,"u10",NF90_REAL,lat_lon_dimids,u10_id))
      call nc_check(nf90_put_att(ncid,u10_id,"units","m/s"))
      call nc_check(nf90_def_var(ncid,"v10",NF90_REAL,lat_lon_dimids,v10_id))
      call nc_check(nf90_put_att(ncid,v10_id,"units","m/s"))
      call nc_check(nf90_def_var(ncid,"gusts10",NF90_REAL,lat_lon_dimids,gusts_id))
      call nc_check(nf90_put_att(ncid,gusts_id,"units","m/s"))
      call nc_check(nf90_enddef(ncid))
      
      ! writing the variables
      call nc_check(nf90_put_var(ncid,start_date_id,start_date))
      call nc_check(nf90_put_var(ncid,start_hour_id,start_hour))
      call nc_check(nf90_put_var(ncid,lat_id,lat_vector))
      call nc_check(nf90_put_var(ncid,lon_id,lon_vector))
      
      call interpolate_to_ll(mslp,lat_lon_output_field,grid)
      call nc_check(nf90_put_var(ncid,mslp_id,lat_lon_output_field))
      
      call interpolate_to_ll(sp,lat_lon_output_field,grid)
      call nc_check(nf90_put_var(ncid,sp_id,lat_lon_output_field))
      
      call interpolate_to_ll(t2,lat_lon_output_field,grid)
      call nc_check(nf90_put_var(ncid,t2_id,lat_lon_output_field))
      
      call interpolate_to_ll(tcc,lat_lon_output_field,grid)
      call nc_check(nf90_put_var(ncid,tcc_id,lat_lon_output_field))
      
      call interpolate_to_ll(rprate,lat_lon_output_field,grid)
      call nc_check(nf90_put_var(ncid,rprate_id,lat_lon_output_field))
      
      call interpolate_to_ll(sprate,lat_lon_output_field,grid)
      call nc_check(nf90_put_var(ncid,sprate_id,lat_lon_output_field))
      
      call interpolate_to_ll(cape,lat_lon_output_field,grid)
      call nc_check(nf90_put_var(ncid,cape_id,lat_lon_output_field))
      
      call interpolate_to_ll(sfc_sw_down,lat_lon_output_field,grid)
      call nc_check(nf90_put_var(ncid,sfc_sw_down_id,lat_lon_output_field))
      
      call interpolate_to_ll(wind_10_m_mean_u_at_cell,lat_lon_output_field,grid)
      call nc_check(nf90_put_var(ncid,u10_id,lat_lon_output_field))
      
      call interpolate_to_ll(wind_10_m_mean_v_at_cell,lat_lon_output_field,grid)
      call nc_check(nf90_put_var(ncid,v10_id,lat_lon_output_field))
      
      call interpolate_to_ll(wind_10_m_gusts_speed_at_cell,lat_lon_output_field, grid)
      call nc_check(nf90_put_var(ncid,gusts_id,lat_lon_output_field))
      
      ! closing the netcdf file
      call nc_check(nf90_close(ncid))
      
      deallocate(wind_10_m_mean_u_at_cell)
      deallocate(wind_10_m_mean_v_at_cell)
      deallocate(wind_10_m_gusts_speed_at_cell)
      deallocate(t2)
      deallocate(mslp)
      deallocate(sp)
      deallocate(rprate)
      deallocate(sprate)
      deallocate(tcc)
      deallocate(cape)
      deallocate(sfc_sw_down)
      
    endif
    
    ! Diagnostics of quantities that are not surface-specific.    
    allocate(div_h_all_layers(n_scalars))
    call div_h(state%wind,div_h_all_layers,grid)
    call calc_rel_vort(state,diag,grid)
    allocate(rel_vort_scalar_field(n_scalars))
    call curl_field_to_cells(diag%rel_vort,rel_vort_scalar_field,grid)
    
    ! Diagnozing the u and v wind components at the vector points.
    allocate(u_at_edge(n_vectors))
    allocate(v_at_edge(n_vectors))
    call calc_uv_at_edge(state%wind,u_at_edge,v_at_edge,grid)
    ! Averaging to cell centers for output.
    allocate(u_at_cell(n_scalars))
    allocate(v_at_cell(n_scalars))
    call edges_to_cells(u_at_edge,u_at_cell,grid)
    call edges_to_cells(v_at_edge,v_at_cell,grid)
    allocate(rh(n_scalars))
    allocate(epv(n_scalars))
    allocate(pressure(n_scalars))
    !$omp parallel do private(ji)
    do ji=1,n_scalars
      if (n_constituents>=4) then
        rh(ji) = 100._wp*rel_humidity(state%rho((n_condensed_constituents + 1)*n_scalars+ji),diag%temperature(ji))
      endif
      pressure(ji) = state%rho(n_condensed_constituents*n_scalars+ji)*gas_constant_diagnostics(state%rho,ji-1)*diag%temperature(ji)
    enddo
    !$omp end parallel do
    
    !$omp parallel do private(ji)
    do ji=1,n_scalars
      diag%scalar_field_placeholder(ji) = state%rho(n_condensed_constituents*n_scalars+ji)
    enddo
    !$omp end parallel do
    
    call calc_pot_vort(state,diag%scalar_field_placeholder,diag,grid)
    call epv_diagnostics(state,diag,epv,grid)
    
    ! pressure level output
    if (lpressure_level_output) then
      ! allocating memory for the variables on pressure levels
      allocate(geopotential_height(n_scalars_h,n_pressure_levels))
      allocate(t_on_p_levels(n_scalars_h,n_pressure_levels))
      allocate(rh_on_p_levels(n_scalars_h,n_pressure_levels))
      allocate(epv_on_p_levels(n_scalars_h,n_pressure_levels))
      allocate(u_on_p_levels(n_scalars_h,n_pressure_levels))
      allocate(v_on_p_levels(n_scalars_h,n_pressure_levels))
      allocate(zeta_on_p_levels(n_scalars_h,n_pressure_levels))
      
      ! vertical interpolation to the pressure levels
      !$omp parallel do private(ji,jl,jm,vector_to_minimize,closest_index,second_closest_index,closest_weight)
      do ji=1,n_scalars_h
        do jl=1,n_pressure_levels
          do jm=1,n_layers
            
            ! It is approx. p = p_0exp(-z/H) => log(p) = log(p_0) - z/H => z/H = log(p_0) - log(p) = log(p_0/p) => z = H*log(p_0/p).
            ! This leads to abs(z_2 - z_1) = abs(H*log(p_2/p) - H*log(p_1/p)) = H*abs(log(p_2/p) - log(p_1/p)) = H*abs(log(p_2/p_1))
            ! propto abs(log(p_2/p_1)).
            
            vector_to_minimize(jm) = abs(log(pressure_levels(jl)/(pressure((jm-1)*n_scalars_h+ji))))
          enddo
          ! finding the model layer that is the closest to the desired pressure level
          closest_index = find_min_index(vector_to_minimize,n_layers)
          ! first guess for the other layer that will be used for the interpolation
          second_closest_index = closest_index + 1
          ! in this case,the layer above the closest layer will be used for the interpolation
          if (pressure_levels(jl)<pressure((closest_index-1)*n_scalars_h+ji)) then
            second_closest_index = closest_index - 1
          endif
          ! in this case,a missing value will be written
          if ((closest_index==n_layers .and. second_closest_index==n_layers+1) &
              .or. (closest_index<0 .or. second_closest_index<0)) then
            geopotential_height(ji,jl) = 9999
            t_on_p_levels(ji,jl) = 9999
            rh_on_p_levels(ji,jl) = 9999
            epv_on_p_levels(ji,jl) = 9999
            zeta_on_p_levels(ji,jl) = 9999
            u_on_p_levels(ji,jl) = 9999
            v_on_p_levels(ji,jl) = 9999
          else
            ! this is the interpolation weight:
            ! closest_weight = 1 - abs((delta z)_{closest})/(abs(z_{closest} - z_{other}))
            
            closest_weight = 1._wp - vector_to_minimize(closest_index)/ &
            (abs(log(pressure((closest_index-1)*n_scalars_h+ji)/pressure((second_closest_index-1)*n_scalars_h+ji))) &
            +EPSILON_SECURITY)
            geopotential_height(ji,jl) = closest_weight*grid%gravity_potential((closest_index-1)*n_scalars_h+ji) &
            + (1._wp - closest_weight)*grid%gravity_potential((second_closest_index-1)*n_scalars_h+ji)
            geopotential_height(ji,jl) = geopotential_height(ji,jl)/gravity
            t_on_p_levels(ji,jl) = closest_weight*diag%temperature((closest_index-1)*n_scalars_h+ji) &
            + (1._wp - closest_weight)*diag%temperature((second_closest_index-1)*n_scalars_h+ji)
            rh_on_p_levels(ji,jl) = closest_weight*rh((closest_index-1)*n_scalars_h+ji) &
            + (1._wp - closest_weight)*rh((second_closest_index-1)*n_scalars_h+ji)
            epv_on_p_levels(ji,jl) = closest_weight*epv((closest_index-1)*n_scalars_h+ji) &
            + (1._wp - closest_weight)*epv((second_closest_index-1)*n_scalars_h+ji)
            zeta_on_p_levels(ji,jl) = closest_weight*rel_vort_scalar_field((closest_index-1)*n_scalars_h+ji) &
            + (1._wp - closest_weight)*rel_vort_scalar_field((second_closest_index-1)*n_scalars_h+ji)
            u_on_p_levels(ji,jl) = closest_weight*u_at_cell((closest_index-1)*n_scalars_h+ji) &
            + (1._wp - closest_weight)*u_at_cell((second_closest_index-1)*n_scalars_h+ji)
            v_on_p_levels(ji,jl) = closest_weight*v_at_cell((closest_index-1)*n_scalars_h+ji) &
            + (1._wp - closest_weight)*v_at_cell((second_closest_index-1)*n_scalars_h+ji)
          endif
        enddo
      enddo
      !$omp end parallel do
      
      output_file_p_level = trim(run_id) // "+" // trim(int2string(time_since_init_min)) // "min_pressure_levels.nc"
      
      
      call nc_check(nf90_create(output_file_p_level,NF90_CLOBBER,ncid))
      call nc_check(nf90_def_dim(ncid,"single_int_index",1,single_int_dimid))
      call nc_check(nf90_def_dim(ncid,"lat_index",n_lat_io_points,lat_dimid))
      call nc_check(nf90_def_dim(ncid,"lon_index",n_lon_io_points,lon_dimid))
      
      lat_lon_dimids(1) = lat_dimid
      lat_lon_dimids(2) = lon_dimid
      
      ! defining the variables
      call nc_check(nf90_def_var(ncid,"start_date",NF90_INT,single_int_dimid,start_date_id))
      call nc_check(nf90_def_var(ncid,"start_hour",NF90_INT,single_int_dimid,start_hour_id))
      call nc_check(nf90_def_var(ncid,"lat",NF90_REAL,lat_dimid,lat_id))
      call nc_check(nf90_def_var(ncid,"lon",NF90_REAL,lon_dimid,lon_id))
      
      do jl=1,n_pressure_levels
        pressure_level_hpa = pressure_levels(jl)/100
        
        varname = "geopot_layer_" // trim(int2string(pressure_level_hpa))
        call nc_check(nf90_def_var(ncid,varname,NF90_REAL,lat_lon_dimids,gh_ids(jl)))
        call nc_check(nf90_put_att(ncid,gh_ids(jl),"units","gpm"))
        
        varname = "temperature_layer_" // trim(int2string(pressure_level_hpa))
        call nc_check(nf90_def_var(ncid,varname,NF90_REAL,lat_lon_dimids,temp_p_ids(jl)))
        call nc_check(nf90_put_att(ncid,temp_p_ids(jl),"units","K"))
        
        varname = "rel_hum_layer_" // trim(int2string(pressure_level_hpa))
        call nc_check(nf90_def_var(ncid,varname,NF90_REAL,lat_lon_dimids,rh_p_ids(jl)))
        call nc_check(nf90_put_att(ncid,rh_p_ids(jl),"units","%"))
        
        varname = "wind_u_layer_" // trim(int2string(pressure_level_hpa))
        call nc_check(nf90_def_var(ncid,varname,NF90_REAL,lat_lon_dimids,wind_u_p_ids(jl)))
        call nc_check(nf90_put_att(ncid,wind_u_p_ids(jl),"units","m/s"))
        
        varname = "wind_v_layer_" // trim(int2string(pressure_level_hpa))
        call nc_check(nf90_def_var(ncid,varname,NF90_REAL,lat_lon_dimids,wind_v_p_ids(jl)))
        call nc_check(nf90_put_att(ncid,wind_v_p_ids(jl),"units","m/s"))
        
        varname = "rel_vort_layer_" // trim(int2string(pressure_level_hpa))
        call nc_check(nf90_def_var(ncid,varname,NF90_REAL,lat_lon_dimids,rel_vort_p_ids(jl)))
        call nc_check(nf90_put_att(ncid,rel_vort_p_ids(jl),"units","1/s"))
        
        varname = "epv_layer_" // trim(int2string(pressure_level_hpa))
        call nc_check(nf90_def_var(ncid,varname,NF90_REAL,lat_lon_dimids,epv_p_ids(jl)))
        call nc_check(nf90_put_att(ncid,epv_p_ids(jl),"units","PVU"))
        
      enddo
      
      call nc_check(nf90_enddef(ncid))
      
      ! writing the variables
      call nc_check(nf90_put_var(ncid,start_date_id,start_date))
      call nc_check(nf90_put_var(ncid,start_hour_id,start_hour))
      call nc_check(nf90_put_var(ncid,lat_id,lat_vector))
      call nc_check(nf90_put_var(ncid,lon_id,lon_vector))
      
      do jl=1,n_pressure_levels
        
        call interpolate_to_ll(geopotential_height(:,jl),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,gh_ids(jl),lat_lon_output_field))
        
        call interpolate_to_ll(t_on_p_levels(:,jl),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,temp_p_ids(jl),lat_lon_output_field))
        
        call interpolate_to_ll(rh_on_p_levels(:,jl),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,rh_p_ids(jl),lat_lon_output_field))
        
        call interpolate_to_ll(u_on_p_levels(:,jl),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,wind_u_p_ids(jl),lat_lon_output_field))
        
        call interpolate_to_ll(v_on_p_levels(:,jl),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,wind_v_p_ids(jl),lat_lon_output_field))
        
        call interpolate_to_ll(zeta_on_p_levels(:,jl),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,rel_vort_p_ids(jl),lat_lon_output_field))
        
        call interpolate_to_ll(epv_on_p_levels(:,jl),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,epv_p_ids(jl),lat_lon_output_field))
      
      enddo
      
      ! closing the netcdf file
      call nc_check(nf90_close(ncid))
      
      deallocate(geopotential_height)
      deallocate(t_on_p_levels)
      deallocate(rh_on_p_levels)
      deallocate(u_on_p_levels)
      deallocate(v_on_p_levels)
      deallocate(epv_on_p_levels)
    endif

    ! model level output
    if (lmodel_level_output) then
    
      output_file = trim(run_id) // "+" // trim(int2string(time_since_init_min)) // "min.nc"
      
      call nc_check(nf90_create(output_file,NF90_CLOBBER,ncid))
      call nc_check(nf90_def_dim(ncid,"single_int_index",1,single_int_dimid))
      call nc_check(nf90_def_dim(ncid,"lat_index",n_lat_io_points,lat_dimid))
      call nc_check(nf90_def_dim(ncid,"lon_index",n_lon_io_points,lon_dimid))
        
      lat_lon_dimids(1) = lat_dimid
      lat_lon_dimids(2) = lon_dimid
      
      ! defining the variables
      call nc_check(nf90_def_var(ncid,"start_date",NF90_INT,single_int_dimid,start_date_id))
      call nc_check(nf90_def_var(ncid,"start_hour",NF90_INT,single_int_dimid,start_hour_id))
      call nc_check(nf90_def_var(ncid,"lat",NF90_REAL,lat_dimid,lat_id))
      call nc_check(nf90_def_var(ncid,"lon",NF90_REAL,lon_dimid,lon_id))
      
      do jl=1,n_layers
      
        varname = "temperature_layer_" // trim(int2string(jl))
        call nc_check(nf90_def_var(ncid,varname,NF90_REAL,lat_lon_dimids,temperature_ids(jl)))
        call nc_check(nf90_put_att(ncid,temperature_ids(jl),"units","K"))
        
        varname = "pressure_layer_" // trim(int2string(jl))
        call nc_check(nf90_def_var(ncid,varname,NF90_REAL,lat_lon_dimids,pressure_ids(jl)))
        call nc_check(nf90_put_att(ncid,temperature_ids(jl),"units","Pa"))
        
        varname = "rel_hum_layer_" // trim(int2string(jl))
        call nc_check(nf90_def_var(ncid,varname,NF90_REAL,lat_lon_dimids,rel_hum_ids(jl)))
        call nc_check(nf90_put_att(ncid,temperature_ids(jl),"units","%"))
        
        varname = "wind_u_layer_" // trim(int2string(jl))
        call nc_check(nf90_def_var(ncid,varname,NF90_REAL,lat_lon_dimids,wind_u_ids(jl)))
        call nc_check(nf90_put_att(ncid,wind_u_ids(jl),"units","m/s"))
        
        varname = "wind_v_layer_" // trim(int2string(jl))
        call nc_check(nf90_def_var(ncid,varname,NF90_REAL,lat_lon_dimids,wind_v_ids(jl)))
        call nc_check(nf90_put_att(ncid,wind_v_ids(jl),"units","m/s"))
        
        varname = "rel_vort_layer_" // trim(int2string(jl))
        call nc_check(nf90_def_var(ncid,varname,NF90_REAL,lat_lon_dimids,rel_vort_ids(jl)))
        call nc_check(nf90_put_att(ncid,rel_vort_ids(jl),"units","1/s"))
        
        varname = "div_h_layer_" // trim(int2string(jl))
        call nc_check(nf90_def_var(ncid,varname,NF90_REAL,lat_lon_dimids,div_h_ids(jl)))
        call nc_check(nf90_put_att(ncid,div_h_ids(jl),"units","1/s"))
        
      enddo
      
      do jl=1,n_levels
        varname = "wind_w_layer_" // trim(int2string(jl))
        call nc_check(nf90_def_var(ncid,varname,NF90_REAL,lat_lon_dimids,wind_w_ids(jl)))
        call nc_check(nf90_put_att(ncid,wind_w_ids(jl),"units","m/s"))
      enddo
      
      call nc_check(nf90_enddef(ncid))
      
      ! writing the variables
      call nc_check(nf90_put_var(ncid,start_date_id,start_date))
      call nc_check(nf90_put_var(ncid,start_hour_id,start_hour))
      call nc_check(nf90_put_var(ncid,lat_id,lat_vector))
      call nc_check(nf90_put_var(ncid,lon_id,lon_vector))
      do jl=1,n_layers
      
        call interpolate_to_ll(diag%temperature(((jl-1)*n_scalars_h+1):(jl*n_scalars_h)),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,temperature_ids(jl),lat_lon_output_field))
        
        call interpolate_to_ll(pressure(((jl-1)*n_scalars_h+1):(jl*n_scalars_h)),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,pressure_ids(jl),lat_lon_output_field))
        
        call interpolate_to_ll(rh(((jl-1)*n_scalars_h+1):(jl*n_scalars_h)),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,rel_hum_ids(jl),lat_lon_output_field))
        
        call interpolate_to_ll(u_at_cell(((jl-1)*n_scalars_h+1):(jl*n_scalars_h)),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,wind_u_ids(jl),lat_lon_output_field))
        
        call interpolate_to_ll(v_at_cell(((jl-1)*n_scalars_h+1):(jl*n_scalars_h)),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,wind_v_ids(jl),lat_lon_output_field))
        
        call interpolate_to_ll(rel_vort_scalar_field(((jl-1)*n_scalars_h+1):(jl*n_scalars_h)),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,rel_vort_ids(jl),lat_lon_output_field))
        
        call interpolate_to_ll(div_h_all_layers(((jl-1)*n_scalars_h+1):(jl*n_scalars_h)),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,div_h_ids(jl),lat_lon_output_field))
        
      enddo
      
      do jl=1,n_levels
        call interpolate_to_ll(state%wind(((jl-1)*n_vectors_per_layer+1):((jl-1)*n_vectors_per_layer+n_scalars_h)), &
                               lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,wind_w_ids(jl),lat_lon_output_field))
      enddo
      
      ! closing the netcdf file
      call nc_check(nf90_close(ncid))
    endif
    
    deallocate(u_at_cell)
    deallocate(v_at_cell)
    ! output of the whole model state for data assimilation
    if ((ideal_input_id==-1 .or. ltotally_first_step) .and. time_since_init_min==time_to_next_analysis_min) then
    
      output_file = trim(run_id) // "+" // trim(int2string(time_since_init_min)) // "min_hex.nc"
      
      call nc_check(nf90_create(output_file,NF90_CLOBBER,ncid))
      call nc_check(nf90_def_dim(ncid,"single_int_index",1,single_int_dimid))
      call nc_check(nf90_def_dim(ncid,"scalar_index",n_scalars,scalar_dimid))
      call nc_check(nf90_def_dim(ncid,"soil_index",nsoillays*n_scalars_h,soil_dimid))
      call nc_check(nf90_def_dim(ncid,"vector_index",n_vectors,vector_dimid))
      call nc_check(nf90_def_dim(ncid,"densities_index",n_constituents*n_scalars,densities_dimid))
      
      ! Defining the variables.
      call nc_check(nf90_def_var(ncid,"start_date",NF90_INT,single_int_dimid,start_date_id))
      call nc_check(nf90_def_var(ncid,"start_hour",NF90_INT,single_int_dimid,start_hour_id))
      call nc_check(nf90_def_var(ncid,"densities",NF90_REAL,densities_dimid,densities_id))
      call nc_check(nf90_put_att(ncid,densities_id,"units","kg/m^3"))
      call nc_check(nf90_def_var(ncid,"temperature",NF90_REAL,scalar_dimid,temperature_id))
      call nc_check(nf90_put_att(ncid,temperature_id,"units","K"))
      call nc_check(nf90_def_var(ncid,"wind",NF90_REAL,vector_dimid,wind_id))
      call nc_check(nf90_put_att(ncid,wind_id,"units","m/s"))
      call nc_check(nf90_def_var(ncid,"tke",NF90_REAL,scalar_dimid,tke_id))
      call nc_check(nf90_put_att(ncid,tke_id,"units","J/kg"))
      call nc_check(nf90_def_var(ncid,"t_soil",NF90_REAL,soil_dimid,soil_id))
      call nc_check(nf90_put_att(ncid,soil_id,"units","K"))
      call nc_check(nf90_enddef(ncid))
      
      ! setting the variables
      call nc_check(nf90_put_var(ncid,start_date_id,start_date))
      call nc_check(nf90_put_var(ncid,start_hour_id,start_hour))
      call nc_check(nf90_put_var(ncid,densities_id,state%rho))
      call nc_check(nf90_put_var(ncid,temperature_id,diag%temperature))
      call nc_check(nf90_put_var(ncid,wind_id,state%wind))
      call nc_check(nf90_put_var(ncid,tke_id,diag%tke))
      call nc_check(nf90_put_var(ncid,soil_id,state%temperature_soil))
      
      ! closing the netcdf file
      call nc_check(nf90_close(ncid))
    endif
    
    deallocate(lat_lon_output_field)
    deallocate(div_h_all_layers)
    deallocate(rel_vort_scalar_field)
    deallocate(rh)
    deallocate(epv)
    deallocate(pressure)
    write(*,*) "Output written."
    
  end subroutine write_out

  function pseudopotential_temperature(state,diag,scalar_index,grid)
    
    ! This function returns the pseudopotential temperature,which is needed for diagnozing CAPE.
    
    type(t_state), intent(in) :: state                       ! state variables
    type(t_diag),  intent(in) :: diag                        ! diagnostic quantities
    integer,       intent(in) :: scalar_index                ! scalar index at which to compute the pseudopotential temperature
    type(t_grid),  intent(in) :: grid                        ! grid properties
    real(wp)                  :: pseudopotential_temperature ! the result
    
    ! local variables
    real(wp) :: r,alpha_1,alpha_2,alpha_3,pressure,t_lcl,vapour_pressure,saturation_pressure,rel_hum
    
    pseudopotential_temperature = 0._wp
    ! the dry case
    if (.not. lmoist) then
      pseudopotential_temperature = grid%theta_v_bg(scalar_index) + state%theta_v_pert(scalar_index)
    ! This is the moist case,based on
    ! Bolton,D. (1980). The Computation of Equivalent Potential Temperature,Monthly Weather Review,108(7),1046-1053.
    else
      
      ! the mixing ratio
      r = state%rho((n_condensed_constituents + 1)*n_scalars + scalar_index) &
      /(state%rho(n_condensed_constituents*n_scalars + scalar_index) &
      - state%rho((n_condensed_constituents + 1)*n_scalars + scalar_index))
      
      ! now,the first two required parameters can already be computed
      alpha_1 = 0.2854_wp*(1._wp - 0.28e-3_wp*r)
      alpha_3 = r*(1._wp + 0.81e-3_wp*r)
      
      ! calculating the pressure
      pressure = p_0*(grid%exner_bg(scalar_index) + state%exner_pert(scalar_index))**(c_d_p/r_d)
      
      ! computing the temperature t_lcl of the air parcel after raising it to the lifted condensation level (LCL)
      ! therefore we firstly compute the saturation pressure,the vapour pressure and the relative humidity
      if (diag%temperature(scalar_index)>=t_0) then
        saturation_pressure = saturation_pressure_over_water(diag%temperature(scalar_index))
      else
        saturation_pressure = saturation_pressure_over_ice(diag%temperature(scalar_index))
      endif
      vapour_pressure = state%rho((n_condensed_constituents + 1)*n_scalars + scalar_index)*r_v*diag%temperature(scalar_index)
      rel_hum = vapour_pressure/saturation_pressure
      ! we compute t_lcl using Eq. (22) of Bolton (1980)
      t_lcl = 1._wp/(1._wp/(diag%temperature(scalar_index) - 55._wp) - log(rel_hum)/2840._wp) + 55._wp
      
      ! the last remaining parameter can be computed now
      alpha_2 = 3.376_wp/t_lcl - 0.00254_wp
      
      ! the final formula by Bolton
      pseudopotential_temperature = diag%temperature(scalar_index)*(p_0/pressure)**alpha_1*exp(alpha_2*alpha_3)
    endif
    
  end function pseudopotential_temperature
  
end module mo_write_output









