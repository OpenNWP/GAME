! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

program control

  ! This program organizes the model, manages the time stepping, calls model output, collects the lowest model layer wind for 10 m wind mean and so on.
  ! All the memory needed for the integration is allocated and freed here.

  use mo_definitions,        only: t_grid,t_state
  use mo_gradient_operators, only: grad

  implicit none
  
  type(t_state) :: state_1,state_2,state_tendency,state_write
  type(t_grid)  :: grid

  ! taking the timestamp to measure the performance
  ! clock_t begin = clock()
  
  ! console output
  write(*,*) "*****************************************************************"
  write(*,*) "*\t\t\t\t\t\t\t\t\t\t*\n"
  write(*,*) "*\t\t\t\tThis is the GAME\t\t\t\t*\n"
  write(*,*) "*\t\t\tGeophysical Fluids Modeling Framework\t\t\t*\n"
  write(*,*) "*\t\t\t\t\t\t\t\t\t\t*\n"
  write(*,*) "*****************************************************************"
  write(*,*) "Released under the MIT license,visit https://github.com/OpenNWP/GAME for more information."
  write(*,*) "*****************************************************************"
  
  ! Allocating memory
  ! ------------------
  
  allocate(grid%normal_distance(n_vectors))
  allocate(grid%volume(n_scalars))
  allocate(grid%z_scalar(n_scalars))
  allocate(grid%from_index(n_vectors_h))
  allocate(grid%to_index(n_vectors_h))
  allocate(state_1%rho(n_constituents*n_scalars))
  allocate(state_1%rhotheta_v(n_scalars))
  allocate(state_1%theta_v_pert(n_scalars))
  allocate(state_1%exner_pert(n_scalars))
  allocate(state_1%wind(n_vectors))
  allocate(state_1%temperature_soil(nsoillays*n_scalars_h))
  allocate(state_2%rho(n_constituents*n_scalars))
  allocate(state_2%rhotheta_v(n_scalars))
  allocate(state_2%theta_v_pert(n_scalars))
  allocate(state_2%exner_pert(n_scalars))
  allocate(state_2%wind(n_vectors))
  allocate(state_2%temperature_soil(nsoillays*n_scalars_h))
  allocate(state_tendency%rho(n_constituents*n_scalars))
  allocate(state_tendency%rhotheta_v(n_scalars))
  allocate(state_tendency%theta_v_pert(n_scalars))
  allocate(state_tendency%exner_pert(n_scalars))
  allocate(state_tendency%wind(n_vectors))
  allocate(state_tendency%temperature_soil(nsoillays*n_scalars_h))
  allocate(diagnostics%flux_density(n_vectors))
  allocate(diagnostics%flux_density_div(n_scalars))
  allocate(diagnostics%rel_vort_on_triangles(n_vectors))
  allocate(diagnostics%rel_vort((2*n_layers+1)*n_vectors_h))
  allocate(diagnostics%pot_vort((2*n_layers+1)*n_vectors_h))
  allocate(diagnostics%temperature(n_scalars))
  allocate(diagnostics%c_g_p_field(n_scalars))
  allocate(diagnostics%v_squared(n_scalars))
  allocate(diagnostics%wind_div(n_scalars))
  allocate(diagnostics%curl_of_vorticity(n_vectors))
  allocate(diagnostics%scalar_field_placeholder(n_scalars))
  allocate(diagnostics%vector_field_placeholder(n_vectors))
  allocate(diagnostics%u_at_edge(n_vectors))
  allocate(diagnostics%v_at_edge(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(diagnostics%rel_vort(n_vectors))
  allocate(state_write%rho(n_constituents*n_scalars))
  allocate(state_write%rhotheta_v(n_scalars))
  allocate(state_write%theta_v_pert(n_scalars))
  allocate(state_write%exner_pert(n_scalars))
  allocate(state_write%wind(n_vectors))
  allocate(state_write%temperature_soil(nsoillays*n_scalars_h))

  ! setting up the namelists (model configuration)
  call grid_nml_setup()
  call run_nml_setup()
  call constituents_nml_setup()
  call diff_nml_setup()
  call surface_nml_setup()
  call io_nml_setup()
  
  ! reading the grid
  write(*,*) "Reading grid data ..."
  call set_grid_properties(grid%normal_distance,grid%volume,grid%area,grid%z_scalar,grid%z_vector, &
                           grid%gravity_potential,grid%theta_v_bg,grid%exner_bg, &
                           grid%layer_thickness,grid%trsk_indices,grid%trsk_modified_curl_indices,grid%from_index, &
                           grid%to_index,grid%adjacent_vector_indices_h,grid%adjacent_signs_h,grid%density_to_rhombi_indices, &
                           grid%latitude_scalar,grid%longitude_scalar,grid%inner_product_weights,grid%direction, &
                           grid%density_to_rhombi_weights,grid%trsk_weights,grid%sfc_albedo,grid%sfc_rho_c, &
                           grid%t_conduc_soil,grid%roughness_length,grid%is_land,grid%latlon_interpol_indices, &
                           grid%latlon_interpol_weights,grid%z_soil_interface,grid%z_soil_center, &
                           grid%t_const_soil,&grid%z_t_const,&grid%toa,&grid%stretching_parameter,&grid%radius, &
                           grid%area_dual,grid%z_vector_dual,grid%normal_distance_dual,grid%from_index_dual, &
                           grid%to_index_dual,grid%vorticity_indices_triangles,grid%vorticity_signs_triangles,grid%f_vec)
  
  call grad_hor_cov(grid%z_scalar,grid%slope,grid%from_index,grid%to_index,grid%normal_distance)
  call grad(grid%gravity_potential,grid%gravity_m,grid%from_index,grid%to_index,grid%normal_distance,grid%inner_product_weights,grid%slope)
  call grad(grid%exner_bg,grid%exner_bg_grad,grid%from_index,grid%to_index,grid%normal_distance,grid%inner_product_weights,grid%slope)
  write(*,*) "Grid loaded successfully."
  
  rad_nml_setup()
  
  ! rescaling times for small Earth experiments
  double radius_rescale = grid%radius/RADIUS
  total_run_span_min = radius_rescale*total_run_span_min
  write_out_interval_min = radius_rescale*write_out_interval_min
  
  ! Reading and processing user input finished.
  
  write(*,*) "Setting initial state ..."
  
  ! ideal test case
  if (ideal_input_id/=-1) then
    call set_ideal_init(state_1%exner_pert,state_1%theta_v_pert,diagnostics%scalar_field_placeholder,grid%exner_bg, &
                        grid%theta_v_bg,grid%adjacent_vector_indices_h,grid%area_dual,grid%density_to_rhombi_indices, &
                        grid%density_to_rhombi_weights,grid%f_vec,diagnostics%flux_density,grid%from_index,grid%to_index, &
                        grid%from_index_dual,grid%to_index_dual,state_1%rho,grid%inner_product_weights,grid%normal_distance, &
                        diagnostics%pot_vort_tend,grid%z_scalar,state_1%rhotheta_v,state_1%wind,diagnostics%v_squared, &
                        grid%direction,grid%latitude_scalar,grid%longitude_scalar,grid%z_vector,grid%slope, &
                        grid%gravity_potential,diagnostics%pot_vort,diagnostics%rel_vort,diagnostics%rel_vort_on_triangles, &
                        grid%trsk_indices,grid%trsk_weights,grid%trsk_modified_curl_indices,grid%z_vector_dual, &
                        grid%vorticity_indices_triangles,grid%vorticity_signs_triangles,grid%t_const_soil, &
                        grid%is_land,state_1%temperature_soil,grid%z_t_const)
  ! NWP mode
  else
    call read_init_data()
  endif
  
  write(*,*) "Initial state set successfully."
  write(*,*) "%s",stars)
  
  if (rad_config>0) then
    write(*,*) "Radiation time step: %lf s\n", radiation_delta_t
  endif
  
  ! finding the minimum horizontal grid distance
  double normal_dist_min_hor = grid%eff_hor_res
  do ji=1,n_vectors_h
    if (grid%normal_distance(n_vectors - n_vectors_per_layer + ji) < normal_dist_min_hor) then
      normal_dist_min_hor = grid%normal_distance(n_vectors - n_vectors_per_layer + ji)
    endif
  enddo
  ! finding the minimum vertical grid distance
  double normal_dist_min_vert = grid%z_vector(1)/n_layers
  do ji=1,n_scalars_h
    if (grid%normal_distance(n_vectors - n_vectors_per_layer - n_scalars_h + ji) < normal_dist_min_vert) then
      normal_dist_min_vert = grid%normal_distance(n_vectors - n_vectors_per_layer - n_scalars_h + ji)
    endif
  enddo
  
  ! setting the hydrometeor falling velcotities
  cloud_droplets_velocity = 0._wp1
  rain_velocity = 10._wp
  snow_velocity = 5._wp
  write(*,*) "Cloud droplets falling velocity set to %lf m/s.\n",cloud_droplets_velocity
  write(*,*) "Rain falling velocity set to %lf m/s.\n",rain_velocity
  write(*,*) "Snow falling velocity set to %lf m/s.\n",snow_velocity
  
  write(*,*) "Effective horizontal resolution: %lf km\n",1e-3*grid%eff_hor_res
  write(*,*) "Minimum horizontal normal distance: %lf km\n",1e-3*normal_dist_min_hor
  double max_speed_hor = 100._wp
  write(*,*) "Horizontal advective Courant number: %lf\n",delta_t/normal_dist_min_hor*max_speed_hor
  double max_speed_vert = 0.1_wp
  write(*,*) "Vertical advective Courant number: %lf\n",delta_t/normal_dist_min_vert*max_speed_vert
  write(*,*) "%s",stars)
  write(*,*) "It begins."
  write(*,*) "%s",stars)
  
  double *wind_h_lowest_layer = calloc(1,min_no_of_10m_wind_avg_steps*n_vectors_h*sizeof(double))
  double t_write = t_init
  #pragma omp parallel for
  do h_index=1,n_vectors_h
    ! here,for all output time steps,the initial value is used
    do time_step_10_m_wind=1,min_no_of_10m_wind_avg_steps
      wind_h_lowest_layer(time_step_10_m_wind*n_vectors_h + h_index) = state_1%wind(n_vectors - n_vectors_per_layer + h_index)
    enddo
  enddo
  call temperature_diagnostics(diagnostics%temperature,grid%theta_v_bg,state_1%theta_v_pert,grid%exner_bg,state_1%exner_pert,state_1%rho)
  call inner_product(state_1%wind,state_1%wind,diagnostics%v_squared,grid%adjacent_vector_indices_h,grid%inner_product_weights)
  
  ! time coordinate of the old RK step
  double t_0
  t_0 = t_init
  
  ! configuring radiation and calculating radiative fluxes for the first time
  rad_update = 1
  double t_rad_update = t_0
  if (rad_config>0) then
    if (rad_config==1) then
      call radiation_init()
    endif
    call call_radiation(grid%latitude_scalar,grid%longitude_scalar,state_1%temperature_soil,grid%sfc_albedo,grid%z_scalar, &
                        grid%z_vector,state_1%rho,diagnostics%temperature,diagnostics%radiation_tendency, &
                        diagnostics%sfc_sw_in,diagnostics%sfc_lw_out,&t_0)
    rad_update = 1
    t_rad_update = t_rad_update + radiation_delta_t
  endif
  
  ! This is necessary because at the very first step of the model integration,some things are handled differently
  ! in the time stepping and in writing the output.
  totally_first_step_bool = 1
  ! writing out the initial state of the model run
  call write_out(diagnostics%scalar_field_placeholder,state_1%wind,grid%latlon_interpol_indices,grid%latlon_interpol_weights,grid%exner_bg,
                 grid%inner_product_weights,grid%volume,grid%gravity_potential,grid%from_index,grid%to_index,grid%z_vector,grid%f_vec,diagnostics%temperature,
                 state_1%temperature_soil,grid%area,state_1%rho,grid%z_scalar,grid%slope,grid%gravity_m,grid%adjacent_signs_h,grid%adjacent_vector_indices_h,
                 grid%area_dual,grid%density_to_rhombi_indices,grid%density_to_rhombi_weights,state_1%exner_pert,diagnostics%tke,t_init,t_write,
                 grid%from_index_dual,grid%to_index_dual,diagnostics%v_squared,grid%is_land,diagnostics%monin_obukhov_length,diagnostics%roughness_velocity,
                 grid%roughness_length,grid%direction,grid%trsk_indices,grid%sfc_albedo,diagnostics%sfc_sw_in,grid%layer_thickness,state_1%theta_v_pert,
                 grid%theta_v_bg,grid%z_vector_dual,grid%vorticity_indices_triangles,grid%vorticity_signs_triangles,grid%trsk_weights,
                 totally_first_step_bool,wind_h_lowest_layer,diagnostics%rel_vort_on_triangles,diagnostics%rel_vort,diagnostics%pot_vort,
                 grid%normal_distance)
  
  t_write = t_write + 60._wp*write_out_interval_min
  write(*,*) "Run progress: %f h\n",(t_init - t_init)/3600)
  time_step_counter = 0
  clock_t first_time,second_time
  first_time = clock()
  double time = 0._wp
  if (write_out_integrals == 1) then
    call write_out_integral(state_1%wind,state_1%rhotheta_v,diagnostics%temperature,state_1%rho, &
                            grid%volume,grid%inner_product_weights,grid%gravity_potential,grid%adjacent_vector_indices_h,&time)
  endif
  
  ! Preparation of the actual integration.
  ! --------------------------------------
  wind_lowest_layer_step_counter = 0
  call linear_combine_two_states(state_1%rho,state_1%rhotheta_v,state_1%exner_pert,state_1%wind,state_1%temperature_soil, &
                                 state_1%rho,state_1%rhotheta_v,state_1%exner_pert,state_1%wind,state_1%temperature_soil, &
                                 state_2%rho,state_2%rhotheta_v,state_2%theta_v_pert,state_2%exner_pert,state_2%wind, &
                                 state_2%temperature_soil,1._wp,0._wp,grid%theta_v_bg)
  
  ! This is the loop over the time steps.
  ! -------------------------------------
  ! this is to store the speed of the model integration
  speed
  while (t_0<t_init + 60*total_run_span_min + radius_rescale*300)
    
    ! Checking if the radiative fluxes need to be updated:
    ! ----------------------------------------------------
    
    if (t_0 <= t_rad_update .and. t_0 + delta_t >= t_rad_update .and. totally_first_step_bool/=1) then
      rad_update = 1
      t_rad_update += radiation_delta_t
    else
      rad_update = 0
    endif
    
    ! Time step integration.
    if (mod(time_step_counter,2)==0) then
      call manage_pchevi(grid%adjacent_signs_h,grid%adjacent_vector_indices_h,grid%area,grid%layer_thickness, &
                         grid%z_scalar,grid%z_vector,grid%volume,grid%vorticity_indices_triangles,grid%vorticity_signs_triangles, &
                         grid%z_t_const,grid%z_soil_center,grid%z_soil_interface,diagnostics%v_squared,grid%trsk_weights, &
                         grid%from_index,grid%to_index,grid%from_index_dual,grid%to_index_dual,grid%trsk_modified_curl_indices, &
                         grid%area_dual,grid%z_vector_dual,state_1%wind,state_tendency%wind,state_2%wind,grid%trsk_indices, &
                         diagnostics%temperature,diagnostics%wind_div,diagnostics%viscosity_triangles,diagnostics%viscosity, &
                         diagnostics%viscosity_rhombi, &
                         diagnostics%condensates_sediment_heat,diagnostics%molecular_diffusion_coeff,diagnostics%v_squared_grad, &
                         t_0,diagnostics%vert_hor_viscosity,diagnostics%vector_field_placeholder,diagnostics%sfc_sw_in, &
                         totally_first_step_bool,grid%gravity_m,diagnostics%curl_of_vorticity,diagnostics%tke,grid%t_const_soil, &
                         state_1%theta_v_pert,state_2%theta_v_pert,grid%theta_v_bg,state_2%temperature_soil,diagnostics%sfc_lw_out, &
                         state_1%temperature_soil,diagnostics%temperature_diffusion_heating,grid%slope,grid%t_conduc_soil, &
                         diagnostics%temp_diffusion_coeff_numerical_h,diagnostics%temp_diffusion_coeff_numerical_v,diagnostics%friction_acc, &
                         grid%sfc_rho_c,grid%sfc_albedo,diagnostics%scalar_flux_resistance,diagnostics%scalar_field_placeholder, &
                         diagnostics%roughness_velocity,grid%roughness_length,state_tendency%rhotheta_v,state_1%rhotheta_v, &
                         state_2%rhotheta_v,state_tendency%rho,state_1%rho,state_2%rho,diagnostics%radiation_tendency,grid%exner_bg, &
                         diagnostics%pot_vort_tend,grid%normal_distance,diagnostics%n_squared,grid%inner_product_weights,grid%exner_bg_grad, &
                         grid%normal_distance_dual,diagnostics%power_flux_density_latent,diagnostics%power_flux_density_sensible, &
                         grid%density_to_rhombi_weights,grid%density_to_rhombi_indices,diagnostics%rel_vort_on_triangles, &
                         diagnostics%phase_trans_heating_rate,state_2%exner_pert,state_1%exner_pert,diagnostics%rel_vort, &
                         grid%f_vec,&rad_update, &
                         diagnostics%pressure_gradient_decel_factor,diagnostics%pressure_gradient_acc_neg_l, &
                         diagnostics%pressure_gradient_acc_neg_nl, &
                         diagnostics%pot_vort,diagnostics%flux_density,diagnostics%pressure_grad_condensates_v, &
                         diagnostics%flux_density_div,diagnostics%dv_hdz, &
                         diagnostics%monin_obukhov_length,diagnostics%heating_diss,grid%is_land,diagnostics%pgrad_acc_old, &
                         diagnostics%mass_diff_tendency, &
                         grid%latitude_scalar,grid%longitude_scalar,diagnostics%mass_diffusion_coeff_numerical_h, &
                         diagnostics%mass_diffusion_coeff_numerical_v,diagnostics%phase_trans_rates)
    else
      call manage_pchevi(grid%adjacent_signs_h,grid%adjacent_vector_indices_h,grid%area,grid%layer_thickness, &
                         grid%z_scalar,grid%z_vector,grid%volume,grid%vorticity_indices_triangles,grid%vorticity_signs_triangles, &
                         grid%z_t_const,grid%z_soil_center,grid%z_soil_interface,diagnostics%v_squared,grid%trsk_weights, &
                         grid%from_index,grid%to_index,grid%from_index_dual,grid%to_index_dual,grid%trsk_modified_curl_indices, &
                         grid%area_dual,grid%z_vector_dual,state_2%wind,state_tendency%wind,state_1%wind,grid%trsk_indices, &
                         diagnostics%temperature,diagnostics%wind_div,diagnostics%viscosity_triangles,diagnostics%viscosity, &
                         diagnostics%viscosity_rhombi, &
                         diagnostics%condensates_sediment_heat,diagnostics%molecular_diffusion_coeff,diagnostics%v_squared_grad, &
                         t_0,diagnostics%vert_hor_viscosity,diagnostics%vector_field_placeholder,diagnostics%sfc_sw_in, &
                         totally_first_step_bool,grid%gravity_m,diagnostics%curl_of_vorticity,diagnostics%tke,grid%t_const_soil, &
                         state_2%theta_v_pert,state_1%theta_v_pert,grid%theta_v_bg,state_1%temperature_soil,diagnostics%sfc_lw_out, &
                         state_2%temperature_soil,diagnostics%temperature_diffusion_heating,grid%slope,grid%t_conduc_soil, &
                         diagnostics%temp_diffusion_coeff_numerical_h,diagnostics%temp_diffusion_coeff_numerical_v,diagnostics%friction_acc, &
                         grid%sfc_rho_c,grid%sfc_albedo,diagnostics%scalar_flux_resistance,diagnostics%scalar_field_placeholder, &
                         diagnostics%roughness_velocity,grid%roughness_length,state_tendency%rhotheta_v,state_2%rhotheta_v, &
                         state_1%rhotheta_v,state_tendency%rho,state_2%rho,state_1%rho,diagnostics%radiation_tendency,grid%exner_bg, &
                         diagnostics%pot_vort_tend,grid%normal_distance,diagnostics%n_squared,grid%inner_product_weights,grid%exner_bg_grad, &
                         grid%normal_distance_dual,diagnostics%power_flux_density_latent,diagnostics%power_flux_density_sensible, &
                         grid%density_to_rhombi_weights,grid%density_to_rhombi_indices,diagnostics%rel_vort_on_triangles, &
                         diagnostics%phase_trans_heating_rate,state_1%exner_pert,state_2%exner_pert,diagnostics%rel_vort, &
                         grid%f_vec,&rad_update, &
                         diagnostics%pressure_gradient_decel_factor,diagnostics%pressure_gradient_acc_neg_l, &
                         diagnostics%pressure_gradient_acc_neg_nl, &
                         diagnostics%pot_vort,diagnostics%flux_density,diagnostics%pressure_grad_condensates_v, &
                         diagnostics%flux_density_div,diagnostics%dv_hdz, &
                         diagnostics%monin_obukhov_length,diagnostics%heating_diss,grid%is_land,diagnostics%pgrad_acc_old, &
                         diagnostics%mass_diff_tendency, &
                         grid%latitude_scalar,grid%longitude_scalar,diagnostics%mass_diffusion_coeff_numerical_h, &
                         diagnostics%mass_diffusion_coeff_numerical_v,diagnostics%phase_trans_rates)
    endif
  
    ! Writing out integrals over the model domain if requested by the user.
    ! ---------------------------------------------------------------------
  
    if (write_out_integrals == 1) then
      if (mod(time_step_counter%2)==0) then
        call write_out_integral(state_2%wind,state_2%rhotheta_v,diagnostics%temperature,state_2%rho,
                                grid%volume,grid%inner_product_weights,grid%gravity_potential,grid%adjacent_vector_indices_h,t_0+delta_t-t_init)
      else
        call write_out_integral(state_1%wind,state_1%rhotheta_v,diagnostics%temperature,state_1%rho,
                                grid%volume,grid%inner_product_weights,grid%gravity_potential,grid%adjacent_vector_indices_h,t_0+delta_t-t_init)
      endif
    endif
  
    ! Writing the actual output.
    ! --------------------------
    
    ! interpolating to the output time
    double new_weight,old_weight
    if(t_0 + delta_t >= t_write .and. t_0 <= t_write) then
      if (fmod(time_step_counter,2) == 0) then
        new_weight = (t_write - t_0)/delta_t
        old_weight = 1._wp - new_weight
        call linear_combine_two_states(state_1%rho,state_1%rhotheta_v,state_1%exner_pert,state_1%wind,state_1%temperature_soil, &
                                       state_2%rho,state_2%rhotheta_v,state_2%exner_pert,state_2%wind,state_2%temperature_soil, &
                                       state_write%rho,state_write%rhotheta_v,state_write%theta_v_pert,state_write%exner_pert, &
                                       state_write%wind,state_write%temperature_soil,old_weight,new_weight,grid%theta_v_bg)
      else
        new_weight = (t_write - t_0)/delta_t
        old_weight = 1._wp - new_weight
        call linear_combine_two_states(state_2%rho,state_2%rhotheta_v,state_2%exner_pert,state_2%wind,state_2%temperature_soil, &
                                       state_1%rho,state_1%rhotheta_v,state_1%exner_pert,state_1%wind,state_1%temperature_soil, &
                                       state_write%rho,state_write%rhotheta_v,state_write%theta_v_pert,state_write%exner_pert, &
                                       state_write%wind,state_write%temperature_soil,old_weight,new_weight,grid%theta_v_bg)
      endif
    endif
  
    ! 5 minutes before the output time,the wind in the lowest layer needs to be collected for 10 m wind diagnostics.
    if (t_0 >= t_write - radius_rescale*300) then
      if (wind_lowest_layer_step_counter < min_no_of_10m_wind_avg_steps) then
        if (mod(time_step_counter%2)==0) then
          !$omp parallel private(h_index)
          do h_index=1,n_vectors_h
            wind_h_lowest_layer(wind_lowest_layer_step_counter*n_vectors_h + h_index) = state_1%wind(N_VECTORS - n_vectors_per_layer + h_index)
          enddo
          !$omp end parallel do
        else
          !$omp parallel do private(h_index)
          do h_index=1,n_vectors_h
            wind_h_lowest_layer(wind_lowest_layer_step_counter*n_vectors_h + h_index) = state_2%wind(N_VECTORS - n_vectors_per_layer + h_index)
          enddo
          !$omp end parallel do
        endif
        wind_lowest_layer_step_counter += 1
      endif
    endif
  
    ! 5 minutes after the output time,the 10 m wind diagnostics can be executed,so output can actually be written
    if(t_0+delta_t>=t_write+radius_rescale*300 .and. t_0<=t_write + radius_rescale*300) then
      ! here,output is actually written
      call write_out(diagnostics%scalar_field_placeholder,state_write%wind,grid%latlon_interpol_indices, &
                     grid%latlon_interpol_weights,grid%exner_bg, &
                     grid%inner_product_weights,grid%volume,grid%gravity_potential,grid%from_index,grid%to_index, &
                     grid%z_vector,grid%f_vec,diagnostics%temperature, &
                     state_write%temperature_soil,grid%area,state_write%rho,grid%z_scalar,grid%slope,grid%gravity_m, &
                     grid%adjacent_signs_h,grid%adjacent_vector_indices_h, &
                     grid%area_dual,grid%density_to_rhombi_indices,grid%density_to_rhombi_weights,state_write%exner_pert, &
                     diagnostics%tke,&t_init,&t_write, &
                     grid%from_index_dual,grid%to_index_dual,diagnostics%v_squared,grid%is_land, &
                     diagnostics%monin_obukhov_length,diagnostics%roughness_velocity, &
                     grid%roughness_length,grid%direction,grid%trsk_indices,grid%sfc_albedo,diagnostics%sfc_sw_in, &
                     grid%layer_thickness,state_write%theta_v_pert, &
                     grid%theta_v_bg,grid%z_vector_dual,grid%vorticity_indices_triangles,grid%vorticity_signs_triangles, &
                     grid%trsk_weights, &
                     totally_first_step_bool,wind_h_lowest_layer,diagnostics%rel_vort_on_triangles,diagnostics%rel_vort, &
                     diagnostics%pot_vort,grid%normal_distance)
      ! setting the next output time
      t_write = t_write + 60._wp*write_out_interval_min
      
      ! Calculating the speed of the model.
      second_time = clock()
      speed = CLOCKS_PER_SEC*60*write_out_interval_min/((double) second_time - first_time)
      write(*,*) "Current speed: %lf\n",speed)
      first_time = clock()
      write(*,*) "Run progress: %f h\n",(t_0 + delta_t - t_init)/3600)
      
      ! resetting the wind in the lowest layer to zero
      !$omp parallel do private(ji)
      do ji=1,min_no_of_10m_wind_avg_steps*n_vectors_h
        wind_h_lowest_layer(ji) = 0._wp
      enddo
      !$omp end parallel do
      wind_lowest_layer_step_counter = 0
    
    endif
    
    
    
    ! This switch can be set to zero now and remains there.
    totally_first_step_bool = 0

    time_step_counter += 1
  
    ! giving the user information on the run progress
    write(*,*) "Step %d completed.\n",time_step_counter)
    
    ! updating the model time
    t_0 = t_0 + delta_t
  
  enddo
  
  ! Clean-up.
  ! ---------
  deallocate(config_io)
  deallocate(diagnostics)
  deallocate(state_tendency)
  deallocate(grid)
  deallocate(state_1)
  deallocate(state_2)
  deallocate(state_write)
  write(*,*) "%s",stars
  deallocate(stars)
  clock_t end = clock()
  speed = CLOCKS_PER_SEC*(60*total_run_span_min + 300)/((double) end - begin)
  deallocate(config)
  write(*,*) "Average speed: %lf\n",speed)
  write(*,*) "GAME over.\n")

end program control










