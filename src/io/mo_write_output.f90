! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_write_output

  ! In this module,the output is written to netCDF files and integrals are written to text files if configured that way.
  ! In addition to that,some postprocessing diagnostics are also calculated here.
  
  use iso_c_binding
  use netcdf
  use mo_various_helpers, only: nc_check
  
  implicit none
  
  contains

  ! the number of pressure levels for the pressure level output
  const int N_PRESSURE_LEVELS = 6

  subroutine interpolate_to_ll(double in_field(),double out_field()(N_LON_IO_POINTS),Grid *grid)
  
    ! This subroutine interpolates a single-layer scalar field to a lat-lon grid.
      
    ! loop over all output points
    !$omp parallel do
    do (int i = 0 i < N_LAT_IO_POINTS ++i)
      do (int j = 0 j < N_LON_IO_POINTS ++j)
        ! initializing the result with zero
        out_field(i)(j) = 0._wp
        ! 1/r-average
        do (int k = 0 k < 5 ++k)
          if (in_field(grid -> latlon_interpol_indices(5*(j + N_LON_IO_POINTS*i) + k))/=9999) then
            out_field(i)(j) += grid -> latlon_interpol_weights(5*(j + N_LON_IO_POINTS*i) + k)*in_field(grid -> latlon_interpol_indices(5*(j + N_LON_IO_POINTS*i) + k))
          else
            out_field(i)(j) = 9999
            break
          endif
        enddo
      enddo
    enddo
    !$omp end parallel do
  
  end subroutine interpolate_to_ll

  subroutine write_out_integral(State *state_write_out,double time_since_init,Grid *grid,Dualgrid *dualgrid,Diagnostics *diagnostics,int integral_id)
  
    
    ! integral_id:
    ! 0: dry mass
    ! 1: entropy
    ! 2: energy
    
    if (integral_id < 0 .or. integral_id > 2) then
      printf("integral_id can only be 0,1 or 2.\n")
      printf("Aborting.\n")
      exit(1)
    endif
    
    double global_integral = 0.0
    FILE *global_integral_file
    
    char integral_file_pre(200)
    
    if (integral_id == 0) then
      sprintf(integral_file_pre,"%s","masses")
    endif
    if (integral_id == 1) then
      sprintf(integral_file_pre,"%s","potential_temperature_density")
    endif
    if (integral_id == 2) then
      sprintf(integral_file_pre,"%s","energy")
    endif
    char integral_file(strlen(integral_file_pre) + 1)
    strcpy(integral_file,integral_file_pre)
    sprintf(integral_file,"%s",integral_file_pre)
    if (integral_id == 0) then
      ! masses
      global_integral_file = fopen(integral_file,"a")
      fprintf(global_integral_file,"%lf\t",time_since_init)
      do (int const_id = 0 const_id < N_CONSTITUENTS ++const_id)
        !$omp parallel do
        do (int i = 0 i < N_SCALARS ++i)
          diagnostics -> scalar_field_placeholder(i) = state_write_out -> rho(const_id*N_SCALARS + i)
        enddo
        global_integral = 0._wp
        do (int i = 0 i < N_SCALARS ++i)
          global_integral += diagnostics -> scalar_field_placeholder(i)*grid -> volume(i)
        enddo
        if (const_id == N_CONSTITUENTS - 1) 
          fprintf(global_integral_file,"%lf\n",global_integral)
        else
          fprintf(global_integral_file,"%lf\t",global_integral)
        endif
      enddo
      fclose(global_integral_file)
    endif
    if (integral_id == 1) then
      ! density times virtual potential temperature
      global_integral_file = fopen(integral_file,"a")
      global_integral = 0._wp
      do (int i = 0 i < N_SCALARS ++i)
        global_integral += state_write_out -> rhotheta_v(i)*grid -> volume(i)
      enddo
      fprintf(global_integral_file,"%lf\t%lf\n",time_since_init,global_integral)
      fclose(global_integral_file)
    endif
    if (integral_id == 2) then
      double kinetic_integral,potential_integral,internal_integral
      global_integral_file = fopen(integral_file,"a")
      Scalar_field *e_kin_density = malloc(sizeof(Scalar_field))
      inner_product(state_write_out -> wind,state_write_out -> wind,*e_kin_density,grid -> adjacent_vector_indices_h,grid -> inner_product_weights)
      !$omp parallel do
      do (int i = 0 i < N_SCALARS ++i)
        diagnostics -> scalar_field_placeholder(i) = state_write_out -> rho(N_CONDENSED_CONSTITUENTS*N_SCALARS + i)
      enddo
      !$omp end parallel do
      !$omp parallel do
      do (int i = 0 i < N_SCALARS ++i)
        (*e_kin_density)(i) = diagnostics -> scalar_field_placeholder(i)*(*e_kin_density)(i)
      enddo
      !$omp end parallel do
      kinetic_integral = 0._wp
      do (int i = 0 i < N_SCALARS ++i)
        kinetic_integral += (*e_kin_density)(i)*grid -> volume(i)
      enddo
      deallocate(e_kin_density)
      Scalar_field *pot_energy_density = malloc(sizeof(Scalar_field))
      !$omp parallel do
      do (int i = 0 i < N_SCALARS ++i)
        (*pot_energy_density)(i) = diagnostics -> scalar_field_placeholder(i)*grid -> gravity_potential(i)
      enddo
      !$omp end parallel do
      potential_integral = 0._wp
      do (int i = 0 i < N_SCALARS ++i)
        potential_integral += (*pot_energy_density)(i)*grid -> volume(i)
      enddo
      deallocate(pot_energy_density)
      Scalar_field *int_energy_density = malloc(sizeof(Scalar_field))
      !$omp parallel do
      do (int i = 0 i < N_SCALARS ++i)
        (*int_energy_density)(i) = diagnostics -> scalar_field_placeholder(i)*diagnostics -> temperature(i)
      enddo
      !$omp end parallel do
      internal_integral = 0._wp
      do (int i = 0 i < N_SCALARS ++i)
        internal_integral += (*int_energy_density)(i)*grid -> volume(i)
      enddo
      fprintf(global_integral_file,"%lf\t%lf\t%lf\t%lf\n",time_since_init,0.5*kinetic_integral,potential_integral,C_D_V*internal_integral)
      deallocate(int_energy_density)
      fclose(global_integral_file)
    endif
    
  end subroutine write_out_integral

  function pseudopotential_temperature(State *state,Diagnostics *diagnostics,Grid *grid,int scalar_index)
    
    ! This function returns the pseudopotential temperature,which is needed for diagnozing CAPE.
    
    double result = 0._wp
    ! the dry case
    if (MOISTURE_ON == 0) then
      result = grid -> theta_v_bg(scalar_index) + state -> theta_v_pert(scalar_index)
    endif
    ! This is the moist case,based on
    ! Bolton,D. (1980). The Computation of Equivalent Potential Temperature,Monthly Weather Review,108(7),1046-1053.
    else
      ! some parameters we need to compute
      double alpha_1,alpha_2,alpha_3,r,pressure,t_lcl
      ! some helper variables
      double vapour_pressure,saturation_pressure,rel_hum
      
      ! the mixing ratio
      r = state -> rho((N_CONDENSED_CONSTITUENTS + 1)*N_SCALARS + scalar_index)
      /(state -> rho(N_CONDENSED_CONSTITUENTS*N_SCALARS + scalar_index)
      - state -> rho((N_CONDENSED_CONSTITUENTS + 1)*N_SCALARS + scalar_index))
      
      ! now,the first two required parameters can already be computed
      alpha_1 = 0.2854_wp*(1._wp - 0.28e-3_wp*r)
      alpha_3 = r*(1._wp + 0.81e-3_wp*r)
      
      ! calculating the pressure
      pressure = P_0*pow(grid -> exner_bg(scalar_index) + state -> exner_pert(scalar_index),C_D_P/R_D)
      
      ! computing the temperature t_lcl of the air parcel after raising it to the lifted condensation level (LCL)
      ! therefore we firstly compute the saturation pressure,the vapour pressure and the relative humidity
      if (diagnostics -> temperature(scalar_index) >= T_0) then
        saturation_pressure = saturation_pressure_over_water(diagnostics -> temperature(scalar_index))
      else
        saturation_pressure = saturation_pressure_over_ice(diagnostics -> temperature(scalar_index))
      endif
      vapour_pressure = state -> rho((N_CONDENSED_CONSTITUENTS + 1)*N_SCALARS + scalar_index)*R_V*diagnostics -> temperature(scalar_index)
      rel_hum = vapour_pressure/saturation_pressure
      ! we compute t_lcl using Eq. (22) of Bolton (1980)
      t_lcl = 1._wp/(1._wp/(diagnostics -> temperature(scalar_index) - 55._wp) - log(rel_hum)/2840._wp) + 55._wp
      
      ! the last remaining parameter can be computed now
      alpha_2 = 3.376_wp/t_lcl - 0.00254_wp
      
      ! the final formula by Bolton
      result = diagnostics -> temperature(scalar_index)*pow(P_0/pressure,alpha_1)*exp(alpha_2*alpha_3)
    enddo
    
  end function pseudopotential_temperature

  subroutine write_out()
  
    write(*,*) "Writing output ..."
    
    int no_of_layers = n_layers
    
    ! latitude resolution of the grid
    double delta_latitude = M_PI/N_LAT_IO_POINTS
    ! longitude resolution of the grid
    double delta_longitude = 2._wp*M_PI/N_LON_IO_POINTS
    
    double lat_vector(N_LAT_IO_POINTS)
    do (int i = 0 i < N_LAT_IO_POINTS ++i)
      lat_vector(i) = M_PI/2._wp - 0.5_wp*delta_latitude - i*delta_latitude
    enddo
    double lon_vector(N_LON_IO_POINTS)
    do (int i = 0 i < N_LON_IO_POINTS ++i)
      lon_vector(i) = i*delta_longitude
    enddo
    
    ! Time stuff.
    time_t t_init_t = (time_t) t_init
    ! t_init is in UTC
    struct tm *p_init_time = gmtime(t_init_t)
    int init_year = p_init_time -> tm_year
    int init_month = p_init_time -> tm_mon
    int init_day = p_init_time -> tm_mday
    int init_hour = p_init_time -> tm_hour
    int init_date = 10000*(init_year + 1900) + 100*(init_month + 1) + init_day
    int init_time = 100*init_hour
    
    ! precipitation rates smaller than this value are set to zero to not confuse users
    double min_precip_rate_mmh = 0.01_wp
    double min_precip_rate = min_precip_rate_mmh/(1000._wp*3600._wp/1024._wp)
    ! this heuristic coefficient converts the cloud water content to cloud cover
    double cloud_water2cloudiness = 10._wp
    
    int layer_index,closest_index,second_closest_index
    double cloud_water_content
    double vector_to_minimize(n_layers)
    
    double (*lat_lon_output_field)(N_LON_IO_POINTS) = malloc(sizeof(double(N_LAT_IO_POINTS)(N_LON_IO_POINTS)))
    
    ! diagnosing the temperature
    temperature_diagnostics(diagnostics -> temperature,grid -> theta_v_bg,state_write_out -> theta_v_pert,grid -> exner_bg,state_write_out -> exner_pert,state_write_out -> rho)
    
    int time_since_init_min = (int) (t_write - t_init)
    time_since_init_min = time_since_init_min/60._wp
    
    ! needed for netcdf
    int ncid,single_int_dimid,lat_dimid,lon_dimid,start_day_id,start_hour_id,lat_id,lon_id
    
    ! Surface output including diagnostics.
    ! -------------------------------------
    
    if (config_io -> surface_output_switch == 1) then
    
      double *mslp = malloc(n_scalars_h*sizeof(double))
      double *sp = malloc(n_scalars_h*sizeof(double))
      double *t2 = malloc(n_scalars_h*sizeof(double))
      double *tcc = malloc(n_scalars_h*sizeof(double))
      double *rprate = malloc(n_scalars_h*sizeof(double))
      double *sprate = malloc(n_scalars_h*sizeof(double))
      double *cape = malloc(n_scalars_h*sizeof(double))
      double *sfc_sw_down = malloc(n_scalars_h*sizeof(double))
      double temp_lowest_layer,pressure_value,mslp_factor,sp_factor,temp_mslp,temp_surface,z_height,theta_v,
      cape_integrand,delta_z,temp_closest,temp_second_closest,delta_z_temp,temperature_gradient,theta_e
      double z_tropopause = 12e3_wp
      double standard_vert_lapse_rate = 0.0065_wp
      !$omp parallel do private(temp_lowest_layer,pressure_value,mslp_factor,sp_factor,temp_mslp,temp_surface, &
      !$omp z_height,theta_v,cape_integrand,delta_z,temp_closest,temp_second_closest,delta_z_temp,temperature_gradient, &
      !$omp theta_e,layer_index,closest_index,second_closest_index,cloud_water_content,vector_to_minimize)
      do (int i = 0 i < n_scalars_h ++i)
        ! Now the aim is to determine the value of the mslp.
        temp_lowest_layer = diagnostics -> temperature((n_layers - 1)*n_scalars_h + i)
        int index = (n_layers - 1)*n_scalars_h + i
        pressure_value = state_write_out -> rho(N_CONDENSED_CONSTITUENTS*N_SCALARS + (n_layers - 1)*n_scalars_h + i)
        *gas_constant_diagnostics(state_write_out -> rho,index)
        *temp_lowest_layer
        temp_mslp = temp_lowest_layer + standard_vert_lapse_rate*grid -> z_scalar(i + (n_layers - 1)*n_scalars_h)
        mslp_factor = pow(1 - (temp_mslp - temp_lowest_layer)/temp_mslp,grid -> gravity_m((n_layers - 1)*n_vectors_per_layer + i)/
        (gas_constant_diagnostics(state_write_out -> rho,index)*standard_vert_lapse_rate))
        mslp(i) = pressure_value/mslp_factor
        
        ! Now the aim is to determine the value of the surface pressure.
        temp_surface = temp_lowest_layer + standard_vert_lapse_rate*(grid -> z_scalar(i + (n_layers - 1)*n_scalars_h) - grid -> z_vector(N_VECTORS - n_scalars_h + i))
        sp_factor = pow(1.0 - (temp_surface - temp_lowest_layer)/temp_surface,grid -> gravity_m((n_layers - 1)*n_vectors_per_layer + i)/
        (gas_constant_diagnostics(state_write_out -> rho,index)*standard_vert_lapse_rate))
        sp(i) = pressure_value/sp_factor
        
        ! Now the aim is to calculate the 2 m temperature.
        do (int j = 0 j < n_layers ++j)
          vector_to_minimize(j) = abs(grid -> z_vector(n_layers*n_vectors_per_layer + i) + 2 - grid -> z_scalar(i + j*n_scalars_h))
        enddo
        closest_index = find_min_index(vector_to_minimize,no_of_layers)
        temp_closest = diagnostics -> temperature(closest_index*n_scalars_h + i)
        delta_z_temp = grid -> z_vector(n_layers*n_vectors_per_layer + i) + 2 - grid -> z_scalar(i + closest_index*n_scalars_h)
        ! real radiation
        if (prog_soil_temp == 1) then
          temperature_gradient = (temp_closest - state_write_out -> temperature_soil(i))
          /(grid -> z_scalar(i + closest_index*n_scalars_h) - grid -> z_vector(n_layers*n_vectors_per_layer + i))
        ! no real radiation
        else
          second_closest_index = closest_index - 1
          if (grid -> z_scalar(i + closest_index*n_scalars_h) > grid -> z_vector(n_layers*n_vectors_per_layer + i) + 2 .and. closest_index < n_layers - 1) then
            second_closest_index = closest_index + 1
          endif
          temp_second_closest = diagnostics -> temperature(second_closest_index*n_scalars_h + i)
          ! calculating the vertical temperature gradient that will be used for the extrapolation
          temperature_gradient = (temp_closest - temp_second_closest)/(grid -> z_scalar(i + closest_index*n_scalars_h) - grid -> z_scalar(i + second_closest_index*n_scalars_h))
        endif
        ! performing the interpolation / extrapolation to two meters above the surface
        t2(i) = temp_closest + delta_z_temp*temperature_gradient
        
        ! diagnozing CAPE
        ! initializing CAPE with zero
        cape(i) = 0._wp
        layer_index = n_layers - 1
        z_height = grid -> z_scalar(layer_index*n_scalars_h + i)
        ! pseduovirtual potential temperature of the particle in the lowest layer
        theta_e = pseudopotential_temperature(state_write_out,diagnostics,grid,layer_index*n_scalars_h + i)
        do while (z_height<z_tropopause)
          ! full virtual potential temperature in the grid box
          theta_v = grid -> theta_v_bg(layer_index*n_scalars_h + i) + state_write_out -> theta_v_pert(layer_index*n_scalars_h + i)
          ! thickness of the gridbox
          delta_z = grid -> layer_thickness(layer_index*n_scalars_h + i)
          ! this is the candidate that we might want to add to the integral
          cape_integrand
          = grid -> gravity_m(layer_index*n_vectors_per_layer + i)*(theta_e - theta_v)/theta_v
          ! we do not add negative values to CAPE (see the definition of CAPE)
          if (cape_integrand>0._wp) then
            cape(i) += cape_integrand*delta_z
          endif
          --layer_index
          z_height = grid -> z_scalar(layer_index*n_scalars_h + i)
        enddo
        
        sfc_sw_down(i) = diagnostics -> sfc_sw_in(i)/(1._wp - grid -> sfc_albedo(i) + EPSILON_SECURITY)
        
        ! Now come the hydrometeors.
        ! Calculation of the total cloud cover
        if (N_CONDENSED_CONSTITUENTS == 4) then
          ! calculating the cloud water content in this column
          cloud_water_content = 0._wp
          do (int k = 0 k < n_layers ++k)
            if (grid -> z_scalar(k*n_scalars_h + i)<z_tropopause) then
              cloud_water_content += (state_write_out -> rho(2*N_SCALARS + k*n_scalars_h + i)
              + state_write_out -> rho(3*N_SCALARS + k*n_scalars_h + i))
              *(grid -> z_vector(i + k*n_vectors_per_layer) - grid -> z_vector(i + (k + 1)*n_vectors_per_layer))
            endif
          enddo
          ! some heuristic ansatz for the total cloud cover
          tcc(i) = min(cloud_water2cloudiness*cloud_water_content,1._wp)
          ! conversion of the total cloud cover into a percentage
          tcc(i) = 100._wp*tcc(i)
          ! setting too small values to zero to not confuse users
          if (tcc(i)<0._wp) then
            tcc(i) = 0._wp
          else
            tcc(i) = 0._wp
          endif
          ! solid precipitation rate
          sprate(i) = 0._wp
          if (N_CONDENSED_CONSTITUENTS==4) then
            sprate(i) = snow_velocity*state_write_out -> rho((n_layers - 1)*n_scalars_h + i)
          endif
          ! liquid precipitation rate
          rprate(i) = 0._wp
          if (N_CONDENSED_CONSTITUENTS==4) then
            rprate(i) = rain_velocity*state_write_out -> rho(N_SCALARS + (n_layers - 1)*n_scalars_h + i)
          endif
          ! setting very small values to zero
          if (rprate(i) < min_precip_rate) then
            rprate(i) = 0._wp
          endif
          ! setting very small values to zero
          if (sprate(i)<min_precip_rate) then
            sprate(i) = 0._wp
          endif
        endif
      enddo
      !$omp end parallel do
      
      ! 10 m wind diagnostics
      ! ---------------------
      
      double wind_tangential,wind_u_value,wind_v_value
      int j
      double *wind_10_m_mean_u = malloc(N_VECS_H*sizeof(double))
      double *wind_10_m_mean_v = malloc(N_VECS_H*sizeof(double))
      ! temporal average over the ten minutes output interval
      !$omp parallel do private(j,wind_tangential,wind_u_value,wind_v_value)
      do (int h_index = 0 h_index < N_VECS_H ++h_index)
        ! initializing the means with zero
        wind_10_m_mean_u(h_index) = 0._wp
        wind_10_m_mean_v(h_index) = 0._wp
        ! loop over the time steps
        do (int time_step_10_m_wind = 0 time_step_10_m_wind < min_no_of_output_steps ++time_step_10_m_wind)
          j = time_step_10_m_wind*N_VECS_H + h_index
          wind_tangential = 0._wp
          do (int i = 0 i < 10 ++i)
            wind_tangential += grid -> trsk_weights(10*h_index + i)*wind_h_lowest_layer_array(time_step_10_m_wind*N_VECS_H + grid -> trsk_indices(10*h_index + i))
          enddo
          wind_10_m_mean_u(h_index) = wind_10_m_mean_u(h_index) + 1._wp/min_no_of_output_steps*wind_h_lowest_layer_array(j)
          wind_10_m_mean_v(h_index) = wind_10_m_mean_v(h_index) + 1._wp/min_no_of_output_steps*wind_tangential
        enddo
        ! passive turn to obtain the u- and v-components of the wind
        double m_direction = -grid -> direction(h_index)
        passive_turn(wind_10_m_mean_u(h_index),wind_10_m_mean_v(h_index),m_direction,wind_u_value,wind_v_value)
        wind_10_m_mean_u(h_index) = wind_u_value
        wind_10_m_mean_v(h_index) = wind_v_value
      enddo
      !$omp end parallel do
      ! vertically extrapolating to ten meters above the surface
      double roughness_length_extrapolation,actual_roughness_length,z_sfc,z_agl,rescale_factor
      !$omp parallel do private(roughness_length_extrapolation,actual_roughness_length,z_sfc,z_agl,rescale_factor)
      do (int i = 0 i < N_VECS_H ++i)
        actual_roughness_length = 0.5_wp*(grid -> roughness_length(grid -> from_index(i)) + grid -> roughness_length(grid -> to_index(i)))
        ! roughness length of grass according to WMO
        roughness_length_extrapolation = 0.02_wp
        if (grid -> is_land(grid -> from_index(i)) == 0) then
          roughness_length_extrapolation = actual_roughness_length
        endif
        z_sfc = 0.5_wp*(grid -> z_vector(N_VECTORS - n_scalars_h + grid -> from_index(i)) + grid -> z_vector(N_VECTORS - n_scalars_h + grid -> to_index(i)))
        z_agl = grid -> z_vector(N_VECTORS - n_vectors_per_layer + i) - z_sfc
        
        ! rescale factor for computing the wind in a height of 10 m
        rescale_factor = log(10._wp/roughness_length_extrapolation)/log(z_agl/actual_roughness_length)
        
        wind_10_m_mean_u(i) = rescale_factor*wind_10_m_mean_u(i)
        wind_10_m_mean_v(i) = rescale_factor*wind_10_m_mean_v(i)
      enddo
      !$omp end parallel do
      
      ! averaging the wind quantities to cell centers for output
      double *wind_10_m_mean_u_at_cell = malloc(n_scalars_h*sizeof(double))
      edges_to_cells_lowest_layer(wind_10_m_mean_u,wind_10_m_mean_u_at_cell,grid -> adjacent_vector_indices_h,grid -> inner_product_weights)
      deallocate(wind_10_m_mean_u)
      double *wind_10_m_mean_v_at_cell = malloc(n_scalars_h*sizeof(double))
      edges_to_cells_lowest_layer(wind_10_m_mean_v,wind_10_m_mean_v_at_cell,grid -> adjacent_vector_indices_h,grid -> inner_product_weights)
      deallocate(wind_10_m_mean_v)
      
      ! gust diagnostics
      double u_850_surrogate,u_950_surrogate
      double u_850_proxy_height = 8000._wp*log(1000._wp/850._wp)
      double u_950_proxy_height = 8000._wp*log(1000._wp/950._wp)
      double *wind_10_m_gusts_speed_at_cell = malloc(n_scalars_h*sizeof(double))
      !$omp parallel do private(closest_index,second_closest_index,u_850_surrogate,u_950_surrogate)
      do (int i = 0 i < n_scalars_h ++i)
      
        ! This is the normal case.
        if ((sfc_sensible_heat_flux == 1 .or. sfc_phase_trans == 1 .or. pbl_scheme == 1)
        .and. abs(diagnostics -> monin_obukhov_length(i)) > EPSILON_SECURITY) then
          ! This follows IFS DOCUMENTATION – Cy43r1 - Operational implementation 22 Nov 2016 - PART IV: PHYSICAL PROCESSES.
          wind_10_m_gusts_speed_at_cell(i) = pow(pow(wind_10_m_mean_u_at_cell(i),2) + pow(wind_10_m_mean_v_at_cell(i),2),0.5_wp)
          + 7.71*diagnostics -> roughness_velocity(i)*pow(max(1._wp - 0.5_wp/12._wp*1000._wp/diagnostics -> monin_obukhov_length(i),0._wp),1._wp/3._wp)
          ! calculating the wind speed in a height representing 850 hPa
          do (int j = 0 j < n_layers ++j)
            vector_to_minimize(j) = abs(grid -> z_scalar(j*n_scalars_h + i) - (grid -> z_vector(N_VECTORS - n_scalars_h + i) + u_850_proxy_height))
          enddo
          closest_index = find_min_index(vector_to_minimize,no_of_layers)
          second_closest_index = closest_index - 1
          if (closest_index < n_layers - 1
          .and. grid -> z_scalar(closest_index*n_scalars_h + i) - grid -> z_vector(N_VECTORS - n_scalars_h + i) > u_850_proxy_height) then
            second_closest_index = closest_index + 1
          endif
          u_850_surrogate = pow(diagnostics -> v_squared(i + closest_index*n_scalars_h),0.5_wp)
          + (pow(diagnostics -> v_squared(i + closest_index*n_scalars_h),0.5_wp) - pow(diagnostics -> v_squared(i + second_closest_index*n_scalars_h),0.5_wp))
          /(grid -> z_scalar(i + closest_index*n_scalars_h) - grid -> z_scalar(i + second_closest_index*n_scalars_h))
          *(grid -> z_vector(N_VECTORS - n_scalars_h + i) + u_850_proxy_height - grid -> z_scalar(i + closest_index*n_scalars_h))
          ! calculating the wind speed in a height representing 950 hPa
          do (int j = 0 j < n_layers ++j)
            vector_to_minimize(j) = abs(grid -> z_scalar(j*n_scalars_h + i) - (grid -> z_vector(N_VECTORS - n_scalars_h + i) + u_950_proxy_height))
          enddo
          closest_index = find_min_index(vector_to_minimize,no_of_layers)
          second_closest_index = closest_index - 1
          if (closest_index < n_layers - 1
          .and. grid -> z_scalar(closest_index*n_scalars_h + i) - grid -> z_vector(N_VECTORS - n_scalars_h + i) > u_950_proxy_height) then
            second_closest_index = closest_index + 1
          endif
          u_950_surrogate = pow(diagnostics -> v_squared(i + closest_index*n_scalars_h),0.5_wp)
          + (pow(diagnostics -> v_squared(i + closest_index*n_scalars_h),0.5_wp) - pow(diagnostics -> v_squared(i + second_closest_index*n_scalars_h),0.5_wp))
          /(grid -> z_scalar(i + closest_index*n_scalars_h) - grid -> z_scalar(i + second_closest_index*n_scalars_h))
          *(grid -> z_vector(N_VECTORS - n_scalars_h + i) + u_950_proxy_height - grid -> z_scalar(i + closest_index*n_scalars_h))
          ! adding the baroclinic and convective component to the gusts
          wind_10_m_gusts_speed_at_cell(i) += 0.6_wp*max(0._wp,u_850_surrogate - u_950_surrogate)
          wind_10_m_gusts_speed_at_cell(i) = min(wind_10_m_gusts_speed_at_cell(i),3._wp*pow(pow(wind_10_m_mean_u_at_cell(i),2) + pow(wind_10_m_mean_v_at_cell(i),2),0.5_wp))
        ! This is used if the turbulence quantities are not populated.
        else
          wind_10_m_gusts_speed_at_cell(i) = 1.67_wp*pow(pow(wind_10_m_mean_u_at_cell(i),2) + pow(wind_10_m_mean_v_at_cell(i),2),0.5_wp)
        endif

      enddo
      !$omp end parallel do
      
      char OUTPUT_FILE_PRE(300)
      sprintf(OUTPUT_FILE_PRE,"%s+%dmin_surface.nc",config_io -> run_id,time_since_init_min)
      char OUTPUT_FILE(strlen(OUTPUT_FILE_PRE) + 1)
      sprintf(OUTPUT_FILE,"%s+%dmin_surface.nc",config_io -> run_id,time_since_init_min)
      int mslp_id,sp_id,rprate_id,sprate_id,
      cape_id,tcc_id,t2_id,u10_id,v10_id,gusts_id,sfc_sw_down_id,start_day_id,start_hour_id
      
      call nc_check(nf90_create(OUTPUT_FILE,NC_CLOBBER,ncid))
      call nc_check(nf90_def_dim(ncid,"single_int_index",1,single_int_dimid))
      call nc_check(nf90_def_dim(ncid,"lat_index",N_LAT_IO_POINTS,lat_dimid))
      call nc_check(nf90_def_dim(ncid,"lon_index",N_LON_IO_POINTS,lon_dimid))
        
      int lat_lon_dimids(2)
      lat_lon_dimids(1) = lat_dimid
      lat_lon_dimids(2) = lon_dimid
      
      ! defining the variables
      call nc_check(nf90_def_var(ncid,"start_day",NC_INT,1,single_int_dimid,start_day_id))
      call nc_check(nf90_def_var(ncid,"start_hour",NC_INT,1,single_int_dimid,start_hour_id))
      call nc_check(nf90_def_var(ncid,"lat",NC_DOUBLE,1,lat_dimid,lat_id))
      call nc_check(nf90_def_var(ncid,"lon",NC_DOUBLE,1,lon_dimid,lon_id))
      call nc_check(nf90_def_var(ncid,"mslp",NC_DOUBLE,2,lat_lon_dimids,mslp_id))
      call nc_check(nf90_put_att(ncid,mslp_id,"units",strlen("Pa"),"Pa"))
      call nc_check(nf90_def_var(ncid,"sp",NC_DOUBLE,2,lat_lon_dimids,sp_id))
      call nc_check(nf90_put_att(ncid,sp_id,"units",strlen("Pa"),"Pa"))
      call nc_check(nf90_def_var(ncid,"t2",NC_DOUBLE,2,lat_lon_dimids,t2_id))
      call nc_check(nf90_put_att(ncid,t2_id,"units",strlen("K"),"K"))
      call nc_check(nf90_def_var(ncid,"tcc",NC_DOUBLE,2,lat_lon_dimids,tcc_id))
      call nc_check(nf90_put_att(ncid,tcc_id,"units",strlen("%"),"%"))
      call nc_check(nf90_def_var(ncid,"rprate",NC_DOUBLE,2,lat_lon_dimids,rprate_id))
      call nc_check(nf90_put_att(ncid,rprate_id,"units",strlen("kg/(m^2s)"),"kg/(m^2s)"))
      call nc_check(nf90_def_var(ncid,"sprate",NC_DOUBLE,2,lat_lon_dimids,sprate_id))
      call nc_check(nf90_put_att(ncid,sprate_id,"units",strlen("kg/(m^2s)"),"kg/(m^2s)"))
      call nc_check(nf90_def_var(ncid,"cape",NC_DOUBLE,2,lat_lon_dimids,cape_id))
      call nc_check(nf90_put_att(ncid,cape_id,"units",strlen("J/kg"),"J/kg"))
      call nc_check(nf90_def_var(ncid,"sfc_sw_down",NC_DOUBLE,2,lat_lon_dimids,sfc_sw_down_id))
      call nc_check(nf90_put_att(ncid,sfc_sw_down_id,"units",strlen("W/m^2"),"W/m^2"))
      call nc_check(nf90_def_var(ncid,"u10",NC_DOUBLE,2,lat_lon_dimids,u10_id))
      call nc_check(nf90_put_att(ncid,u10_id,"units",strlen("m/s"),"m/s"))
      call nc_check(nf90_def_var(ncid,"v10",NC_DOUBLE,2,lat_lon_dimids,v10_id))
      call nc_check(nf90_put_att(ncid,v10_id,"units",strlen("m/s"),"m/s"))
      call nc_check(nf90_def_var(ncid,"gusts10",NC_DOUBLE,2,lat_lon_dimids,gusts_id))
      call nc_check(nf90_put_att(ncid,gusts_id,"units",strlen("m/s"),"m/s"))
      call nc_check(nc_enddef(ncid))
      
      ! writing the variables
      call nc_check(nf90_put_var(ncid,start_day_id,init_date))
      call nc_check(nf90_put_var(ncid,start_hour_id,init_time))
      call nc_check(nf90_put_var(ncid,lat_id,lat_vector))
      call nc_check(nf90_put_var(ncid,lon_id,lon_vector))
      call interpolate_to_ll(mslp,lat_lon_output_field,grid)
      call nc_check(nf90_put_var(ncid,mslp_id,lat_lon_output_field))
      call interpolate_to_ll(sp,lat_lon_output_field,grid)
      call nc_check(nf90_put_var(ncid,sp_id,lat_lon_output_field))
      call interpolate_to_ll(t2,lat_lon_output_field,grid)
      call nc_check(nf90_put_var(ncid,t2_id,lat_lon_output_field))
      call interpolate_to_ll(tcc,lat_lon_output_field,grid)
      call nc_check(nf90_put_var(ncid,tcc_id,lat_lon_output_field)
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
      call interpolate_to_ll(wind_10_m_gusts_speed_at_cell,lat_lon_output_field,grid)
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
    Scalar_field *div_h_all_layers = calloc(1,sizeof(Scalar_field))
    div_h(state_write_out -> wind,*div_h_all_layers,
    grid -> adjacent_signs_h,grid -> adjacent_vector_indices_h,grid -> inner_product_weights,grid -> slope,grid -> area,grid -> volume)
    calc_rel_vort(state_write_out->wind,diagnostics->rel_vort_on_triangles,grid->z_vector,dualgrid->z_vector,diagnostics->rel_vort,
                  dualgrid->vorticity_indices_triangles,dualgrid->vorticity_signs_triangles,grid->normal_distance,
                  dualgrid->area,grid->from_index,grid->to_index,dualgrid->from_index,dualgrid->to_index,grid->inner_product_weights,
                  grid->slope)
    Scalar_field *rel_vort = calloc(1,sizeof(Scalar_field))
    curl_field_to_cells(diagnostics -> rel_vort,*rel_vort,grid -> adjacent_vector_indices_h,grid -> inner_product_weights)
    
    ! Diagnozing the u and v wind components at the vector points.
    calc_uv_at_edge(state_write_out -> wind,diagnostics -> u_at_edge,diagnostics -> v_at_edge,grid -> trsk_indices,grid -> trsk_weights,grid -> direction)
    ! Averaging to cell centers for output.
    edges_to_cells(diagnostics -> u_at_edge,diagnostics -> u_at_cell,grid -> adjacent_vector_indices_h,grid -> inner_product_weights)
    edges_to_cells(diagnostics -> v_at_edge,diagnostics -> v_at_cell,grid -> adjacent_vector_indices_h,grid -> inner_product_weights)
    Scalar_field *rh = calloc(1,sizeof(Scalar_field))
    Scalar_field *epv = calloc(1,sizeof(Scalar_field))
    Scalar_field *pressure = calloc(1,sizeof(Scalar_field))
    !$omp parallel do
    do (int i = 0 i < N_SCALARS ++i)
      if (N_CONSTITUENTS >= 4) then
        (*rh)(i) = 100._wp*rel_humidity(state_write_out -> rho((N_CONDENSED_CONSTITUENTS + 1)*N_SCALARS + i),diagnostics -> temperature(i))
      endif
      (*pressure)(i) = state_write_out -> rho(N_CONDENSED_CONSTITUENTS*N_SCALARS + i)*gas_constant_diagnostics(state_write_out -> rho,i)*diagnostics -> temperature(i)
    enddo
    !$omp end parallel do
    
    !$omp parallel do
    do (int i = 0 i < N_SCALARS ++i)
      diagnostics -> scalar_field_placeholder(i) = state_write_out -> rho(N_CONDENSED_CONSTITUENTS*N_SCALARS + i)
    enddo
    !$omp end parallel do
    calc_pot_vort(state_write_out->wind,diagnostics->rel_vort_on_triangles,grid->z_vector,dualgrid->z_vector,diagnostics->rel_vort,
                  dualgrid->vorticity_indices_triangles,dualgrid->vorticity_signs_triangles,grid->normal_distance,
                  dualgrid->area,grid->from_index,grid->to_index,dualgrid->from_index,dualgrid->to_index,grid->inner_product_weights,
                  grid->slope,dualgrid->f_vec,diagnostics -> pot_vort,grid -> density_to_rhombi_indices,grid -> density_to_rhombi_weights,
                  diagnostics -> scalar_field_placeholder)
    epv_diagnostics(epv,grid -> from_index,grid -> to_index,grid -> inner_product_weights,diagnostics -> pot_vort,grid -> trsk_indices,grid -> trsk_weights,
                    grid -> adjacent_vector_indices_h,grid -> slope,grid -> normal_distance,grid -> theta_v_bg,state_write_out -> theta_v_pert,grid -> z_vector)
    
    ! pressure level output
    double closest_weight
    if (config_io -> pressure_level_output_switch == 1) then
      double *pressure_levels = malloc(sizeof(double)*N_PRESSURE_LEVELS)
      pressure_levels(1) = 20000._wp
      pressure_levels(2) = 30000._wp
      pressure_levels(3) = 50000._wp
      pressure_levels(4) = 70000._wp
      pressure_levels(5) = 85000._wp
      pressure_levels(6) = 92500._wp
      ! allocating memory for the variables on pressure levels
      double (*geopotential_height)(n_scalars_h) = malloc(sizeof(double(N_PRESSURE_LEVELS)(n_scalars_h)))
      double (*t_on_pressure_levels)(n_scalars_h) = malloc(sizeof(double(N_PRESSURE_LEVELS)(n_scalars_h)))
      double (*rh_on_pressure_levels)(n_scalars_h) = malloc(sizeof(double(N_PRESSURE_LEVELS)(n_scalars_h)))
      double (*epv_on_pressure_levels)(n_scalars_h) = malloc(sizeof(double(N_PRESSURE_LEVELS)(n_scalars_h)))
      double (*u_on_pressure_levels)(n_scalars_h) = malloc(sizeof(double(N_PRESSURE_LEVELS)(n_scalars_h)))
      double (*v_on_pressure_levels)(n_scalars_h) = malloc(sizeof(double(N_PRESSURE_LEVELS)(n_scalars_h)))
      double (*rel_vort_on_pressure_levels)(n_scalars_h) = malloc(sizeof(double(N_PRESSURE_LEVELS)(n_scalars_h)))
      
      ! vertical interpolation to the pressure levels
      !$omp parallel do private(vector_to_minimize,closest_index,second_closest_index,closest_weight)
      do (int j = 0 j < N_PRESSURE_LEVELS ++j)
        do (int i = 0 i < n_scalars_h ++i)
          do (int k = 0 k < n_layers ++k)
            
            ! It is approx. p = p_0exp(-z/H) => log(p) = log(p_0) - z/H => z/H = log(p_0) - log(p) = log(p_0/p) => z = H*log(p_0/p).
            ! This leads to abs(z_2 - z_1) = abs(H*log(p_2/p) - H*log(p_1/p)) = H*abs(log(p_2/p) - log(p_1/p)) = H*abs(log(p_2/p_1))
            ! propto abs(log(p_2/p_1)).
            
            vector_to_minimize(k) = abs(log(pressure_levels(j)/(*pressure)(k*n_scalars_h + i)))
          enddo
          ! finding the model layer that is the closest to the desired pressure level
          closest_index = find_min_index(vector_to_minimize,no_of_layers)
          ! first guess for the other layer that will be used for the interpolation
          second_closest_index = closest_index + 1
          ! in this case,the layer above the closest layer will be used for the interpolation
          if (pressure_levels(j) < (*pressure)(closest_index*n_scalars_h + i)) then
            second_closest_index = closest_index - 1
          endif
          ! in this case,a missing value will be written
          if ((closest_index == n_layers - 1 .and. second_closest_index == n_layers) .or. (closest_index < 0 .or. second_closest_index < 0)) then
            geopotential_height(j)(i) = 9999
            t_on_pressure_levels(j)(i) = 9999
            rh_on_pressure_levels(j)(i) = 9999
            epv_on_pressure_levels(j)(i) = 9999
            rel_vort_on_pressure_levels(j)(i) = 9999
            u_on_pressure_levels(j)(i) = 9999
            v_on_pressure_levels(j)(i) = 9999
          else
            ! this is the interpolation weight:
            ! closest_weight = 1 - abs((delta z)_{closest})/(abs(z_{closest} - z_{other}))
            
            closest_weight = 1._wp - vector_to_minimize(closest_index)/
            (abs(log((*pressure)(closest_index*n_scalars_h + i)/(*pressure)(second_closest_index*n_scalars_h + i))) + EPSILON_SECURITY)
            geopotential_height(j)(i) = closest_weight*grid -> gravity_potential(closest_index*n_scalars_h + i)
            + (1._wp - closest_weight)*grid -> gravity_potential(second_closest_index*n_scalars_h + i)
            geopotential_height(j)(i) = geopotential_height(j)(i)/G_MEAN_SFC_ABS
            t_on_pressure_levels(j)(i) = closest_weight*diagnostics -> temperature(closest_index*n_scalars_h + i)
            + (1._wp - closest_weight)*diagnostics -> temperature(second_closest_index*n_scalars_h + i)
            rh_on_pressure_levels(j)(i) = closest_weight*(*rh)(closest_index*n_scalars_h + i)
            + (1._wp - closest_weight)*(*rh)(second_closest_index*n_scalars_h + i)
            epv_on_pressure_levels(j)(i) = closest_weight*(*epv)(closest_index*n_scalars_h + i)
            + (1._wp - closest_weight)*(*epv)(second_closest_index*n_scalars_h + i)
            rel_vort_on_pressure_levels(j)(i) = closest_weight*(*rel_vort)(closest_index*n_scalars_h + i)
            + (1._wp - closest_weight)*(*rel_vort)(second_closest_index*n_scalars_h + i)
            u_on_pressure_levels(j)(i) = closest_weight*diagnostics-> u_at_cell(closest_index*n_scalars_h + i)
            + (1._wp - closest_weight)*diagnostics-> u_at_cell(second_closest_index*n_scalars_h + i)
            v_on_pressure_levels(j)(i) = closest_weight*diagnostics-> v_at_cell(closest_index*n_scalars_h + i)
            + (1._wp - closest_weight)*diagnostics-> v_at_cell(second_closest_index*n_scalars_h + i)
          endif
        enddo
      enddo
      !$omp end parallel do
      
      int OUTPUT_FILE_PRESSURE_LEVEL_LENGTH = 300
      char *OUTPUT_FILE_PRESSURE_LEVEL_PRE = malloc((OUTPUT_FILE_PRESSURE_LEVEL_LENGTH + 1)*sizeof(char))
      sprintf(OUTPUT_FILE_PRESSURE_LEVEL_PRE,"%s+%dmin_pressure_levels.nc",config_io -> run_id,time_since_init_min)
      OUTPUT_FILE_PRESSURE_LEVEL_LENGTH = strlen(OUTPUT_FILE_PRESSURE_LEVEL_PRE)
      deallocate(OUTPUT_FILE_PRESSURE_LEVEL_PRE)
      char *OUTPUT_FILE_PRESSURE_LEVEL = malloc((OUTPUT_FILE_PRESSURE_LEVEL_LENGTH + 1)*sizeof(char))
      sprintf(OUTPUT_FILE_PRESSURE_LEVEL,"%s+%dmin_pressure_levels.nc",config_io -> run_id,time_since_init_min)
      
      int gh_ids(N_PRESSURE_LEVELS),temp_p_ids(N_PRESSURE_LEVELS),rh_p_ids(N_PRESSURE_LEVELS),
      wind_u_p_ids(N_PRESSURE_LEVELS),wind_v_p_ids(N_PRESSURE_LEVELS),
      epv_p_ids(N_PRESSURE_LEVELS),rel_vort_p_ids(N_PRESSURE_LEVELS)
      
      
      call nc_check(nf90_create(OUTPUT_FILE_PRESSURE_LEVEL,NC_CLOBBER,ncid))
      call nc_check(nf90_def_dim(ncid,"single_int_index",1,single_int_dimid))
      call nc_check(nf90_def_dim(ncid,"lat_index",N_LAT_IO_POINTS,lat_dimid))
      call nc_check(nf90_def_dim(ncid,"lon_index",N_LON_IO_POINTS,lon_dimid))
        
      int lat_lon_dimids(2)
      lat_lon_dimids(1) = lat_dimid
      lat_lon_dimids(2) = lon_dimid
      
      ! defining the variables
      call nc_check(nf90_def_var(ncid,"start_day",NC_INT,1,single_int_dimid,start_day_id))
      call nc_check(nf90_def_var(ncid,"start_hour",NC_INT,1,single_int_dimid,start_hour_id))
      call nc_check(nf90_def_var(ncid,"lat",NC_DOUBLE,1,lat_dimid,lat_id))
      call nc_check(nf90_def_var(ncid,"lon",NC_DOUBLE,1,lon_dimid,lon_id))
      
      char varname(100)
      do (int i = 0 i < N_PRESSURE_LEVELS ++i)
        int pressure_level_hpa = (int) pressure_levels(i)/100._wp
        sprintf(varname,"geopot_layer_%d",pressure_level_hpa)
        call nc_check(nf90_def_var(ncid,varname,NC_DOUBLE,2,lat_lon_dimids,gh_ids(i)))
        call nc_check(nf90_put_att(ncid,gh_ids(i),"units",strlen("gpm"),"gpm"))
        sprintf(varname,"temperature_layer_%d",pressure_level_hpa)
        call nc_check(nf90_def_var(ncid,varname,NC_DOUBLE,2,lat_lon_dimids,temp_p_ids(i)))
        call nc_check(nf90_put_att(ncid,temp_p_ids(i),"units",strlen("K"),"K"))
        sprintf(varname,"rel_hum_layer_%d",pressure_level_hpa)
        call nc_check(nf90_def_var(ncid,varname,NC_DOUBLE,2,lat_lon_dimids,rh_p_ids(i)))
        call nc_check(nf90_put_att(ncid,rh_p_ids(i),"units",strlen("%"),"%"))
        sprintf(varname,"wind_u_layer_%d",pressure_level_hpa)
        call nc_check(nf90_def_var(ncid,varname,NC_DOUBLE,2,lat_lon_dimids,wind_u_p_ids(i)))
        call nc_check(nf90_put_att(ncid,wind_u_p_ids(i),"units",strlen("m/s"),"m/s"))
        sprintf(varname,"wind_v_layer_%d",pressure_level_hpa)
        call nc_check(nf90_def_var(ncid,varname,NC_DOUBLE,2,lat_lon_dimids,wind_v_p_ids(i)))
        call nc_check(nf90_put_att(ncid,wind_v_p_ids(i),"units",strlen("m/s"),"m/s"))
        sprintf(varname,"rel_vort_layer_%d",pressure_level_hpa)
        call nc_check(nf90_def_var(ncid,varname,NC_DOUBLE,2,lat_lon_dimids,epv_p_ids(i)))
        call nc_check(nf90_put_att(ncid,epv_p_ids(i),"units",strlen("PVU"),"PVU"))
        sprintf(varname,"epv_layer_%d",pressure_level_hpa)
        call nc_check(nf90_def_var(ncid,varname,NC_DOUBLE,2,lat_lon_dimids,rel_vort_p_ids(i)))
        call nc_check(nf90_put_att(ncid,rel_vort_p_ids(i),"units",strlen("K*m^2/(ks*s)"),"K*m^2/(ks*s)"))
      enddo
      
      call nc_check(nc_enddef(ncid))
      
      ! writing the variables
      call nc_check(nf90_put_var(ncid,start_day_id,init_date))
      call nc_check(nf90_put_var(ncid,start_hour_id,init_time))
      call nc_check(nf90_put_var(ncid,lat_id,lat_vector(0)))
      call nc_check(nf90_put_var(ncid,lon_id,lon_vector(0)))
      do (int i = 0 i < N_PRESSURE_LEVELS ++i)
        
        call interpolate_to_ll(geopotential_height(i)(0),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,gh_ids(i),lat_lon_output_field(0)(0)))
        call interpolate_to_ll(t_on_pressure_levels(i)(0),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,temp_p_ids(i),lat_lon_output_field(0)(0)))
        call interpolate_to_ll(rh_on_pressure_levels(i)(0),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,rh_p_ids(i),lat_lon_output_field(0)(0)))
        call interpolate_to_ll(u_on_pressure_levels(i)(0),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,wind_u_p_ids(i),lat_lon_output_field(0)(0)))
        call interpolate_to_ll(v_on_pressure_levels(i)(0),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,wind_v_p_ids(i),lat_lon_output_field(0)(0)))
        call interpolate_to_ll(epv_on_pressure_levels(i)(0),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,epv_p_ids(i),lat_lon_output_field(0)(0)))
        call interpolate_to_ll(rel_vort_on_pressure_levels(i)(0),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,rel_vort_p_ids(i),lat_lon_output_field(0)(0)))
      
      enddo
      
      ! closing the netcdf file
      call nc_check(nf90_close(ncid))
      
      deallocate(OUTPUT_FILE_PRESSURE_LEVEL)
      deallocate(geopotential_height)
      deallocate(t_on_pressure_levels)
      deallocate(rh_on_pressure_levels)
      deallocate(u_on_pressure_levels)
      deallocate(v_on_pressure_levels)
      deallocate(epv_on_pressure_levels)
      deallocate(pressure_levels)
    endif

    ! model level output
    if (config_io -> model_level_output_switch == 1) then
    
      char OUTPUT_FILE_PRE(300)
      sprintf(OUTPUT_FILE_PRE,"%s+%dmin.nc",config_io -> run_id,time_since_init_min)
      char OUTPUT_FILE(strlen(OUTPUT_FILE_PRE) + 1)
      sprintf(OUTPUT_FILE,"%s+%dmin.nc",config_io -> run_id,time_since_init_min)
      
      int temperature_ids(n_layers),pressure_ids(n_layers),rel_hum_ids(n_layers),
      wind_u_ids(n_layers),wind_v_ids(n_layers),
      rel_vort_ids(n_layers),div_h_ids(n_layers),wind_w_ids(N_LEVELS)
      
      call nc_check(nf90_create(OUTPUT_FILE,NC_CLOBBER,ncid))
      call nc_check(nf90_def_dim(ncid,"single_int_index",1,single_int_dimid))
      call nc_check(nf90_def_dim(ncid,"lat_index",N_LAT_IO_POINTS,lat_dimid))
      call nc_check(nf90_def_dim(ncid,"lon_index",N_LON_IO_POINTS,lon_dimid))
        
      int lat_lon_dimids(2)
      lat_lon_dimids(0) = lat_dimid
      lat_lon_dimids(1) = lon_dimid
      
      ! defining the variables
      call nc_check(nf90_def_var(ncid,"start_day",NC_INT,1,single_int_dimid,start_day_id))
      call nc_check(nf90_def_var(ncid,"start_hour",NC_INT,1,single_int_dimid,start_hour_id))
      call nc_check(nf90_def_var(ncid,"lat",NC_DOUBLE,1,lat_dimid,lat_id))
      call nc_check(nf90_def_var(ncid,"lon",NC_DOUBLE,1,lon_dimid,lon_id))
      
      char varname(100)
      do (int i = 0 i < n_layers ++i)
        sprintf(varname,"temperature_layer_%d",i)
        call nc_check(nf90_def_var(ncid,varname,NC_DOUBLE,2,lat_lon_dimids,temperature_ids(i)))
        call nc_check(nf90_put_att(ncid,temperature_ids(i),"units",strlen("K"),"K"))
        sprintf(varname,"pressure_layer_%d",i)
        call nc_check(nf90_def_var(ncid,varname,NC_DOUBLE,2,lat_lon_dimids,pressure_ids(i)))
        call nc_check(nf90_put_att(ncid,temperature_ids(i),"units",strlen("Pa"),"Pa"))
        sprintf(varname,"rel_hum_layer_%d",i)
        call nc_check(nf90_def_var(ncid,varname,NC_DOUBLE,2,lat_lon_dimids,rel_hum_ids(i)))
        call nc_check(nf90_put_att(ncid,temperature_ids(i),"units",strlen("%"),"%"))
        sprintf(varname,"wind_u_layer_%d",i)
        call nc_check(nf90_def_var(ncid,varname,NC_DOUBLE,2,lat_lon_dimids,wind_u_ids(i)))
        call nc_check(nf90_put_att(ncid,wind_u_ids(i),"units",strlen("m/s"),"m/s"))
        sprintf(varname,"wind_v_layer_%d",i)
        call nc_check(nf90_def_var(ncid,varname,NC_DOUBLE,2,lat_lon_dimids,wind_v_ids(i)))
        call nc_check(nf90_put_att(ncid,wind_v_ids(i),"units",strlen("m/s"),"m/s"))
        sprintf(varname,"rel_vort_layer_%d",i)
        call nc_check(nf90_def_var(ncid,varname,NC_DOUBLE,2,lat_lon_dimids,rel_vort_ids(i)))
        call nc_check(nf90_put_att(ncid,rel_vort_ids(i),"units",strlen("1/s"),"1/s"))
        sprintf(varname,"div_h_layer_%d",i)
        call nc_check(nf90_def_var(ncid,varname,NC_DOUBLE,2,lat_lon_dimids,div_h_ids(i)))
        call nc_check(nf90_put_att(ncid,div_h_ids(i),"units",strlen("1/s"),"1/s"))
      enddo
      do (int i = 0 i < N_LEVELS ++i)
        sprintf(varname,"wind_w_layer_%d",i)
        call nc_check(nf90_def_var(ncid,varname,NC_DOUBLE,2,lat_lon_dimids,wind_w_ids(i)))
        call nc_check(nf90_put_att(ncid,wind_w_ids(i),"units",strlen("m/s"),"m/s"))
      enddo
      call nc_check(nc_enddef(ncid))
      
      ! writing the variables
      call nc_check(nf90_put_var(ncid,start_day_id,init_date))
      call nc_check(nf90_put_var(ncid,start_hour_id,init_time))
      call nc_check(nf90_put_var(ncid,lat_id,lat_vector(0)))
      call nc_check(nf90_put_var(ncid,lon_id,lon_vector(0)))
      do (int i = 0 i < n_layers ++i)
        interpolate_to_ll(diagnostics -> temperature(i*n_scalars_h),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,temperature_ids(i),lat_lon_output_field(0)(0)))
        interpolate_to_ll((*pressure)(i*n_scalars_h),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,pressure_ids(i),lat_lon_output_field(0)(0)))
        interpolate_to_ll((*rh)(i*n_scalars_h),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,rel_hum_ids(i),lat_lon_output_field(0)(0)))
        interpolate_to_ll(diagnostics-> u_at_cell(i*n_scalars_h),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,wind_u_ids(i),lat_lon_output_field(0)(0)))
        interpolate_to_ll(diagnostics-> v_at_cell(i*n_scalars_h),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,wind_v_ids(i),lat_lon_output_field(0)(0)))
        interpolate_to_ll((*rel_vort)(i*n_scalars_h),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,rel_vort_ids(i),lat_lon_output_field(0)(0)))
        interpolate_to_ll((*div_h_all_layers)(i*n_scalars_h),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,div_h_ids(i),lat_lon_output_field(0)(0)))
      enddo
      
      do (int i = 0 i < N_LEVELS ++i)
        interpolate_to_ll(state_write_out -> wind(i*n_vectors_per_layer),lat_lon_output_field,grid)
        call nc_check(nf90_put_var(ncid,wind_w_ids(i),lat_lon_output_field(0)(0)))
      enddo
      
      ! closing the netcdf file
      call nc_check(nf90_close(ncid))
    endif
    
    ! output of the whole model state for data assimilation
    if ((config_io -> ideal_input_id == -1 .or. totally_first_step_bool == 1)
    .and. time_since_init_min == time_to_next_analysis_min) then
    
      char OUTPUT_FILE_PRE(300)
      sprintf(OUTPUT_FILE_PRE,"%s+%dmin_hex.nc",config_io -> run_id,time_since_init_min)
      char OUTPUT_FILE(strlen(OUTPUT_FILE_PRE) + 1)
      sprintf(OUTPUT_FILE,"%s+%dmin_hex.nc",config_io -> run_id,time_since_init_min)
      int scalar_dimid,soil_dimid,vector_dimid,densities_dimid,densities_id,temperature_id,wind_id,
      tke_id,soil_id,single_int_dimid
      
      call nc_check(nf90_create(OUTPUT_FILE,NC_CLOBBER,ncid))
      call nc_check(nf90_def_dim(ncid,"single_int_index",1,single_int_dimid))
      call nc_check(nf90_def_dim(ncid,"scalar_index",N_SCALARS,scalar_dimid))
      call nc_check(nf90_def_dim(ncid,"soil_index",N_SOIL_LAYERS*n_scalars_h,soil_dimid))
      call nc_check(nf90_def_dim(ncid,"vector_index",N_VECTORS,vector_dimid))
      call nc_check(nf90_def_dim(ncid,"densities_index",N_CONSTITUENTS*N_SCALARS,densities_dimid))
      
      ! Defining the variables.
      call nc_check(nf90_def_var(ncid,"start_day",NC_INT,1,single_int_dimid,start_day_id))
      call nc_check(nf90_def_var(ncid,"start_hour",NC_INT,1,single_int_dimid,start_hour_id))
      call nc_check(nf90_def_var(ncid,"densities",NC_DOUBLE,1,densities_dimid,densities_id))
      call nc_check(nf90_put_att(ncid,densities_id,"units",strlen("kg/m^3"),"kg/m^3"))
      call nc_check(nf90_def_var(ncid,"temperature",NC_DOUBLE,1,scalar_dimid,temperature_id))
      call nc_check(nf90_put_att(ncid,temperature_id,"units",strlen("K"),"K"))
      call nc_check(nf90_def_var(ncid,"wind",NC_DOUBLE,1,vector_dimid,wind_id))
      call nc_check(nf90_put_att(ncid,wind_id,"units",strlen("m/s"),"m/s"))
      call nc_check(nf90_def_var(ncid,"tke",NC_DOUBLE,1,scalar_dimid,tke_id))
      call nc_check(nf90_put_att(ncid,tke_id,"units",strlen("J/kg"),"J/kg"))
      call nc_check(nf90_def_var(ncid,"t_soil",NC_DOUBLE,1,soil_dimid,soil_id))
      call nc_check(nf90_put_att(ncid,soil_id,"units",strlen("K"),"K"))
      call nc_check(nc_enddef(ncid))
      
      ! setting the variables
      call nc_check(nf90_put_var(ncid,start_day_id,init_date))
      call nc_check(nf90_put_var(ncid,start_hour_id,init_time))
      call nc_check(nf90_put_var(ncid,densities_id,state_write_out -> rho(0)))
      call nc_check(nf90_put_var(ncid,temperature_id,diagnostics -> temperature(0)))
      call nc_check(nf90_put_var(ncid,wind_id,state_write_out -> wind(0)))
      call nc_check(nf90_put_var(ncid,tke_id,diagnostics -> tke(0)))
      call nc_check(nf90_put_var(ncid,soil_id,state_write_out -> temperature_soil(0)))
      
      ! closing the netcdf file
      call nc_check(nf90_close(ncid))
    endif
    
    deallocate(lat_lon_output_field)
    deallocate(div_h_all_layers)
    deallocate(rel_vort)
    deallocate(rh)
    deallocate(epv)
    deallocate(pressure)
    write(*,*) "Output written."
    
  end subroutine write_out


end module mo_write_output









