! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

program control

  ! This program organizes the model, manages the time stepping, calls model output, collects the lowest model layer wind for 10 m wind mean and so on.
  ! All the memory needed for the integration is allocated and freed here.

  use mo_definitions,            only: t_grid,t_state,t_diag,wp
  use mo_grid_nml,               only: n_scalars,n_layers,n_scalars_h,n_vectors,n_vectors_h,n_dual_vectors, &
                                       n_dual_scalars_h,n_dual_v_vectors,n_h_vectors,n_latlon_io_points, &
                                       n_vectors_per_layer,grid_nml_setup
  use mo_gradient_operators,     only: grad_hor_cov,grad
  use mo_inner_product,          only: inner_product
  use mo_constituents_nml,       only: cloud_droplets_velocity,rain_velocity,snow_velocity,n_constituents, &
                                       n_condensed_constituents,constituents_nml_setup
  use mo_run_nml,                only: dtime,ideal_input_id,total_run_span_min,run_nml_setup
  use mo_grid_setup,             only: eff_hor_res,radius_rescale,set_grid_properties
  use mo_io_nml,                 only: n_output_steps_10m_wind,lwrite_integrals,write_out_interval_min, &
                                       io_nml_setup
  use mo_surface_nml,            only: nsoillays,surface_nml_setup
  use mo_rad_nml,                only: rad_config,radiation_dtime,rad_nml_setup
  use mo_diff_nml,               only: diff_nml_setup
  use mo_set_initial_state,      only: set_ideal_init,read_init_data
  use mo_derived,                only: temperature_diagnostics
  use mo_write_output,           only: write_out,write_out_integrals
  use mo_linear_combination,     only: linear_combine_two_states
  use mo_manage_radiation_calls, only: call_radiation
  use mo_rrtmgp_coupler,         only: radiation_init
  use mo_manage_pchevi,          only: manage_pchevi

  implicit none
  
  type(t_state)         :: state_1,state_2,state_tendency,state_write
  type(t_diag)          :: diag
  type(t_grid)          :: grid
  integer               :: ji,time_step_counter,h_index,rad_update,wind_lowest_layer_step_counter,time_step_10_m_wind, &
                           totally_first_step_bool
  real(wp)              :: t_init,t_0,t_write,t_rad_update,new_weight,old_weight,max_speed_hor,max_speed_ver, &
                           normal_dist_min_hor,normal_dist_min_ver
  real(wp), allocatable :: wind_h_lowest_layer(:)
  character(len=128)    :: init_state_file

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
  allocate(grid%area(n_vectors))
  allocate(grid%z_scalar(n_scalars))
  allocate(grid%z_vector(n_vectors))
  allocate(grid%gravity_potential(n_vectors))
  allocate(grid%gravity_m(n_vectors))
  allocate(grid%slope(n_vectors))
  allocate(grid%theta_v_bg(n_scalars))
  allocate(grid%exner_bg(n_scalars))
  allocate(grid%exner_bg_grad(n_vectors))
  allocate(grid%layer_thickness(n_scalars))
  allocate(grid%area_dual(n_dual_vectors))
  allocate(grid%z_vector_dual(n_dual_vectors))
  allocate(grid%normal_distance_dual(n_dual_vectors))
  allocate(grid%vorticity_indices_triangles(3*n_dual_scalars_h))
  allocate(grid%vorticity_signs_triangles(3*n_dual_scalars_h))
  allocate(grid%f_vec(2*n_vectors_h))
  allocate(grid%trsk_indices(10*n_vectors_h))
  allocate(grid%trsk_modified_curl_indices(10*n_vectors_h))
  allocate(grid%from_index(n_vectors_h))
  allocate(grid%to_index(n_vectors_h))
  allocate(grid%from_index_dual(n_vectors_h))
  allocate(grid%to_index_dual(n_vectors_h))
  allocate(grid%adjacent_vector_indices_h(6*n_scalars_h))
  allocate(grid%adjacent_signs_h(6*n_scalars_h))
  allocate(grid%density_to_rhombi_indices(4*n_vectors_h))
  allocate(grid%latitude_scalar(n_scalars_h))
  allocate(grid%longitude_scalar(n_scalars_h))
  allocate(grid%inner_product_weights(8*n_scalars))
  allocate(grid%direction(n_vectors_h))
  allocate(grid%density_to_rhombi_weights(4*n_vectors_h))
  allocate(grid%trsk_weights(10*n_vectors_h))
  allocate(grid%sfc_albedo(n_scalars_h))
  allocate(grid%sfc_rho_c(n_scalars_h))
  allocate(grid%t_conduc_soil(n_scalars_h))
  allocate(grid%roughness_length(n_scalars_h))
  allocate(grid%is_land(n_scalars_h))
  allocate(grid%latlon_interpol_indices(5*n_latlon_io_points))
  allocate(grid%latlon_interpol_weights(5*n_latlon_io_points))
  allocate(grid%z_soil_interface(nsoillays+1))
  allocate(grid%z_soil_center(nsoillays))
  allocate(grid%t_const_soil(n_scalars_h))
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
  allocate(diag%flux_density(n_vectors))
  allocate(diag%flux_density_div(n_scalars))
  allocate(diag%rel_vort_on_triangles(n_vectors))
  allocate(diag%rel_vort((2*n_layers+1)*n_vectors_h))
  allocate(diag%pot_vort((2*n_layers+1)*n_vectors_h))
  allocate(diag%temperature(n_scalars))
  allocate(diag%c_g_p_field(n_scalars))
  allocate(diag%v_squared(n_scalars))
  allocate(diag%wind_div(n_scalars))
  allocate(diag%curl_of_vorticity(n_vectors))
  allocate(diag%scalar_field_placeholder(n_scalars))
  allocate(diag%vector_field_placeholder(n_vectors))
  allocate(diag%u_at_edge(n_vectors))
  allocate(diag%v_at_edge(n_vectors))
  allocate(diag%u_at_cell(n_scalars))
  allocate(diag%v_at_cell(n_scalars))
  allocate(diag%n_squared(n_scalars))
  allocate(diag%dv_hdz(n_h_vectors+n_vectors_h))
  allocate(diag%scalar_flux_resistance(n_scalars_h))
  allocate(diag%power_flux_density_sensible(n_scalars_h))
  allocate(diag%power_flux_density_latent(n_scalars_h))
  allocate(diag%roughness_velocity(n_scalars_h))
  allocate(diag%monin_obukhov_length(n_scalars_h))
  allocate(diag%temperature_diffusion_heating(n_scalars))
  allocate(diag%friction_acc(n_vectors))
  allocate(diag%heating_diss(n_scalars))
  allocate(diag%molecular_diffusion_coeff(n_scalars))
  allocate(diag%mass_diffusion_coeff_numerical_h(n_scalars))
  allocate(diag%mass_diffusion_coeff_numerical_v(n_scalars))
  allocate(diag%temp_diffusion_coeff_numerical_h(n_scalars))
  allocate(diag%temp_diffusion_coeff_numerical_v(n_scalars))
  allocate(diag%pressure_gradient_decel_factor(n_scalars))
  allocate(diag%condensates_sediment_heat(n_scalars))
  allocate(diag%mass_diff_tendency(n_constituents*n_scalars))
  allocate(diag%phase_trans_rates((n_condensed_constituents+1)*n_scalars))
  allocate(diag%phase_trans_heating_rate(n_scalars))
  allocate(diag%viscosity(n_scalars))
  allocate(diag%viscosity_rhombi(n_vectors))
  allocate(diag%viscosity_triangles(n_dual_v_vectors))
  allocate(diag%vert_hor_viscosity(n_h_vectors+n_vectors_h))
  allocate(diag%tke(n_scalars))
  allocate(diag%pgrad_acc_old(n_vectors))
  allocate(diag%pressure_gradient_acc_neg_nl(n_vectors))
  allocate(diag%pressure_gradient_acc_neg_l(n_vectors))
  allocate(diag%pressure_grad_condensates_v(n_vectors))
  allocate(diag%v_squared_grad(n_vectors))
  allocate(diag%pot_vort_tend(n_vectors))
  allocate(diag%sfc_sw_in(n_scalars_h))
  allocate(diag%sfc_lw_out(n_scalars_h))
  allocate(diag%radiation_tendency(n_scalars))
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
  call rad_nml_setup()
  
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
                           grid%t_const_soil, &
                           grid%area_dual,grid%z_vector_dual,grid%normal_distance_dual,grid%from_index_dual, &
                           grid%to_index_dual,grid%vorticity_indices_triangles,grid%vorticity_signs_triangles,grid%f_vec)
  
  call grad_hor_cov(grid%z_scalar,grid%slope,grid%from_index,grid%to_index,grid%normal_distance)
  call grad(grid%gravity_potential,grid%gravity_m,grid%from_index,grid%to_index, &
            grid%normal_distance,grid%inner_product_weights,grid%slope)
  call grad(grid%exner_bg,grid%exner_bg_grad,grid%from_index,grid%to_index, &
            grid%normal_distance,grid%inner_product_weights,grid%slope)
  write(*,*) "Grid loaded successfully."
  
  ! rescaling times for small Earth experiments
  total_run_span_min = radius_rescale*total_run_span_min
  write_out_interval_min = radius_rescale*write_out_interval_min
  
  ! Reading and processing user input finished.
  
  write(*,*) "Setting initial state ..."
  
  ! ideal test case
  if (ideal_input_id/=-1) then
    call set_ideal_init(state_1%exner_pert,state_1%theta_v_pert,diag%scalar_field_placeholder,grid%exner_bg, &
                        grid%theta_v_bg,grid%adjacent_vector_indices_h,grid%area_dual,grid%density_to_rhombi_indices, &
                        grid%density_to_rhombi_weights,grid%f_vec,diag%flux_density,grid%from_index,grid%to_index, &
                        grid%from_index_dual,grid%to_index_dual,state_1%rho,grid%inner_product_weights,grid%normal_distance, &
                        diag%pot_vort_tend,grid%z_scalar,state_1%rhotheta_v,state_1%wind,diag%v_squared, &
                        grid%direction,grid%latitude_scalar,grid%longitude_scalar,grid%z_vector,grid%slope, &
                        grid%gravity_potential,diag%pot_vort,diag%rel_vort,diag%rel_vort_on_triangles, &
                        grid%trsk_indices,grid%trsk_weights,grid%trsk_modified_curl_indices,grid%z_vector_dual, &
                        grid%vorticity_indices_triangles,grid%vorticity_signs_triangles,grid%t_const_soil, &
                        grid%is_land,state_1%temperature_soil)
  ! NWP mode
  else
    init_state_file = "a"
    call read_init_data(init_state_file,state_1%rho,state_1%wind,state_1%rhotheta_v,state_1%theta_v_pert, &
                        state_1%exner_pert,grid%t_const_soil,diag%tke, &
                        state_1%temperature_soil,grid%is_land,grid%theta_v_bg,grid%exner_bg)
  endif
  
  write(*,*) "Initial state set successfully."
  write(*,*) "**"
  
  if (rad_config>0) then
    write(*,*) "Radiation time step: %lf s\n",radiation_dtime
  endif
  
  ! finding the minimum horizontal grid distance
  normal_dist_min_hor = eff_hor_res
  do ji=1,n_vectors_h
    if (grid%normal_distance(n_vectors - n_vectors_per_layer + ji)<normal_dist_min_hor) then
      normal_dist_min_hor = grid%normal_distance(n_vectors - n_vectors_per_layer + ji)
    endif
  enddo
  ! finding the minimum vertical grid distance
  normal_dist_min_ver = grid%z_vector(1)/n_layers
  do ji=1,n_scalars_h
    if (grid%normal_distance(n_vectors - n_vectors_per_layer - n_scalars_h + ji)<normal_dist_min_ver) then
      normal_dist_min_ver = grid%normal_distance(n_vectors - n_vectors_per_layer - n_scalars_h + ji)
    endif
  enddo
  
  write(*,*) "Cloud droplets falling velocity set to %lf m/s.\n",cloud_droplets_velocity
  write(*,*) "Rain falling velocity set to %lf m/s.\n",rain_velocity
  write(*,*) "Snow falling velocity set to %lf m/s.\n",snow_velocity
  
  write(*,*) "Effective horizontal resolution: %lf km\n",1e-3*eff_hor_res
  write(*,*) "Minimum horizontal normal distance: %lf km\n",1e-3*normal_dist_min_hor
  max_speed_hor = 100._wp
  write(*,*) "Horizontal advective Courant number: %lf\n",dtime/normal_dist_min_hor*max_speed_hor
  max_speed_ver = 0.1_wp
  write(*,*) "Vertical advective Courant number: %lf\n",dtime/normal_dist_min_ver*max_speed_ver
  write(*,*) "**"
  write(*,*) "It begins."
  write(*,*) "**"
  
  allocate(wind_h_lowest_layer(n_output_steps_10m_wind*n_vectors_h))
  !$omp parallel do private(h_index)
  do h_index=1,n_vectors_h
    ! here,for all output time steps,the initial value is used
    do time_step_10_m_wind=1,n_output_steps_10m_wind
      wind_h_lowest_layer(time_step_10_m_wind*n_vectors_h + h_index) = state_1%wind(n_vectors - n_vectors_per_layer + h_index)
    enddo
  enddo
  !$omp end parallel do
  call temperature_diagnostics(diag%temperature,grid%theta_v_bg,state_1%theta_v_pert,grid%exner_bg,state_1%exner_pert,state_1%rho)
  call inner_product(state_1%wind,state_1%wind,diag%v_squared,grid%adjacent_vector_indices_h,grid%inner_product_weights)
  
  ! time coordinate of the old RK step
  t_0 = t_init
  
  ! configuring radiation and calculating radiative fluxes for the first time
  rad_update = 1
  t_rad_update = t_0
  if (rad_config>0) then
    if (rad_config==1) then
      call radiation_init()
    endif
    call call_radiation(grid%latitude_scalar,grid%longitude_scalar,state_1%temperature_soil,grid%sfc_albedo,grid%z_scalar, &
                        grid%z_vector,state_1%rho,diag%temperature,diag%radiation_tendency, &
                        diag%sfc_sw_in,diag%sfc_lw_out,t_0)
    rad_update = 1
    t_rad_update = t_rad_update+radiation_dtime
  endif
  
  ! This is necessary because at the very first step of the model integration,some things are handled differently
  ! in the time stepping and in writing the output.
  totally_first_step_bool = 1
  ! writing out the initial state of the model run
  call write_out(diag%scalar_field_placeholder,state_1%wind,grid%latlon_interpol_indices, &
  grid%latlon_interpol_weights,grid%exner_bg, &
                 grid%inner_product_weights,grid%volume,grid%gravity_potential, &
                 grid%from_index,grid%to_index,grid%z_vector,grid%f_vec,diag%temperature, &
                 state_1%temperature_soil,grid%area,state_1%rho,grid%z_scalar,grid%slope, &
                 grid%gravity_m,grid%adjacent_signs_h,grid%adjacent_vector_indices_h, &
                 grid%area_dual,grid%density_to_rhombi_indices,grid%density_to_rhombi_weights, &
                 state_1%exner_pert,diag%tke,t_init,t_write, &
                 grid%from_index_dual,grid%to_index_dual,diag%v_squared,grid%is_land, &
                 diag%monin_obukhov_length,diag%roughness_velocity, &
                 grid%roughness_length,grid%direction,grid%trsk_indices,grid%sfc_albedo, &
                 diag%sfc_sw_in,grid%layer_thickness,state_1%theta_v_pert, &
                 grid%theta_v_bg,grid%z_vector_dual,grid%vorticity_indices_triangles, &
                 grid%vorticity_signs_triangles,grid%trsk_weights, &
                 totally_first_step_bool,wind_h_lowest_layer,diag%rel_vort_on_triangles,diag%rel_vort,diag%pot_vort, &
                 grid%normal_distance)
  
  t_write = t_write + 60._wp*write_out_interval_min
  write(*,*) "Run progress:", (t_init - t_init)/3600._wp, "h"
  ! clock_t first_time,second_time
  ! first_time = clock()
  if (lwrite_integrals) then
    call write_out_integrals(state_1%wind,state_1%rhotheta_v,diag%temperature,state_1%rho, &
                             grid%volume,grid%inner_product_weights,grid%gravity_potential,grid%adjacent_vector_indices_h,0._wp)
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
  do while (t_0<t_init+60._wp*total_run_span_min+radius_rescale*300._wp)
    
    ! Checking if the radiative fluxes need to be updated:
    ! ----------------------------------------------------
    
    if (t_0<=t_rad_update .and. t_0+dtime>=t_rad_update .and. totally_first_step_bool/=1) then
      rad_update = 1
      t_rad_update = t_rad_update+radiation_dtime
    else
      rad_update = 0
    endif
    
    ! time step integration
    if (mod(time_step_counter,2)==0) then
      call manage_pchevi(grid%adjacent_signs_h,grid%adjacent_vector_indices_h,grid%area,grid%layer_thickness, &
                         grid%z_scalar,grid%z_vector,grid%volume,grid%vorticity_indices_triangles,grid%vorticity_signs_triangles, &
                         grid%z_soil_center,grid%z_soil_interface,diag%v_squared,grid%trsk_weights, &
                         grid%from_index,grid%to_index,grid%from_index_dual,grid%to_index_dual,grid%trsk_modified_curl_indices, &
                         grid%area_dual,grid%z_vector_dual,state_1%wind,state_tendency%wind,state_2%wind,grid%trsk_indices, &
                         diag%temperature,diag%wind_div,diag%viscosity_triangles,diag%viscosity, &
                         diag%viscosity_rhombi, &
                         diag%condensates_sediment_heat,diag%molecular_diffusion_coeff,diag%v_squared_grad, &
                         t_0,diag%vert_hor_viscosity,diag%vector_field_placeholder,diag%sfc_sw_in, &
                         totally_first_step_bool,grid%gravity_m,diag%curl_of_vorticity,diag%tke,grid%t_const_soil, &
                         state_1%theta_v_pert,state_2%theta_v_pert,grid%theta_v_bg,state_2%temperature_soil,diag%sfc_lw_out, &
                         state_1%temperature_soil,diag%temperature_diffusion_heating,grid%slope,grid%t_conduc_soil, &
                         diag%temp_diffusion_coeff_numerical_h,diag%temp_diffusion_coeff_numerical_v,diag%friction_acc, &
                         grid%sfc_rho_c,grid%sfc_albedo,diag%scalar_flux_resistance,diag%scalar_field_placeholder, &
                         diag%roughness_velocity,grid%roughness_length,state_tendency%rhotheta_v,state_1%rhotheta_v, &
                         state_2%rhotheta_v,state_tendency%rho,state_1%rho,state_2%rho,diag%radiation_tendency,grid%exner_bg, &
                         diag%pot_vort_tend,grid%normal_distance,diag%n_squared,grid%inner_product_weights,grid%exner_bg_grad, &
                         grid%normal_distance_dual,diag%power_flux_density_latent,diag%power_flux_density_sensible, &
                         grid%density_to_rhombi_weights,grid%density_to_rhombi_indices,diag%rel_vort_on_triangles, &
                         diag%phase_trans_heating_rate,state_2%exner_pert,state_1%exner_pert,diag%rel_vort, &
                         grid%f_vec,rad_update, &
                         diag%pressure_gradient_decel_factor,diag%pressure_gradient_acc_neg_l, &
                         diag%pressure_gradient_acc_neg_nl, &
                         diag%pot_vort,diag%flux_density,diag%pressure_grad_condensates_v, &
                         diag%flux_density_div,diag%dv_hdz, &
                         diag%monin_obukhov_length,diag%heating_diss,grid%is_land,diag%pgrad_acc_old, &
                         diag%mass_diff_tendency, &
                         grid%latitude_scalar,grid%longitude_scalar,diag%mass_diffusion_coeff_numerical_h, &
                         diag%mass_diffusion_coeff_numerical_v,diag%phase_trans_rates)
    else
      call manage_pchevi(grid%adjacent_signs_h,grid%adjacent_vector_indices_h,grid%area,grid%layer_thickness, &
                         grid%z_scalar,grid%z_vector,grid%volume,grid%vorticity_indices_triangles,grid%vorticity_signs_triangles, &
                         grid%z_soil_center,grid%z_soil_interface,diag%v_squared,grid%trsk_weights, &
                         grid%from_index,grid%to_index,grid%from_index_dual,grid%to_index_dual,grid%trsk_modified_curl_indices, &
                         grid%area_dual,grid%z_vector_dual,state_2%wind,state_tendency%wind,state_1%wind,grid%trsk_indices, &
                         diag%temperature,diag%wind_div,diag%viscosity_triangles,diag%viscosity, &
                         diag%viscosity_rhombi, &
                         diag%condensates_sediment_heat,diag%molecular_diffusion_coeff,diag%v_squared_grad, &
                         t_0,diag%vert_hor_viscosity,diag%vector_field_placeholder,diag%sfc_sw_in, &
                         totally_first_step_bool,grid%gravity_m,diag%curl_of_vorticity,diag%tke,grid%t_const_soil, &
                         state_2%theta_v_pert,state_1%theta_v_pert,grid%theta_v_bg,state_1%temperature_soil,diag%sfc_lw_out, &
                         state_2%temperature_soil,diag%temperature_diffusion_heating,grid%slope,grid%t_conduc_soil, &
                         diag%temp_diffusion_coeff_numerical_h,diag%temp_diffusion_coeff_numerical_v,diag%friction_acc, &
                         grid%sfc_rho_c,grid%sfc_albedo,diag%scalar_flux_resistance,diag%scalar_field_placeholder, &
                         diag%roughness_velocity,grid%roughness_length,state_tendency%rhotheta_v,state_2%rhotheta_v, &
                         state_1%rhotheta_v,state_tendency%rho,state_2%rho,state_1%rho,diag%radiation_tendency,grid%exner_bg, &
                         diag%pot_vort_tend,grid%normal_distance,diag%n_squared,grid%inner_product_weights,grid%exner_bg_grad, &
                         grid%normal_distance_dual,diag%power_flux_density_latent,diag%power_flux_density_sensible, &
                         grid%density_to_rhombi_weights,grid%density_to_rhombi_indices,diag%rel_vort_on_triangles, &
                         diag%phase_trans_heating_rate,state_1%exner_pert,state_2%exner_pert,diag%rel_vort, &
                         grid%f_vec,rad_update, &
                         diag%pressure_gradient_decel_factor,diag%pressure_gradient_acc_neg_l, &
                         diag%pressure_gradient_acc_neg_nl, &
                         diag%pot_vort,diag%flux_density,diag%pressure_grad_condensates_v, &
                         diag%flux_density_div,diag%dv_hdz, &
                         diag%monin_obukhov_length,diag%heating_diss,grid%is_land,diag%pgrad_acc_old, &
                         diag%mass_diff_tendency, &
                         grid%latitude_scalar,grid%longitude_scalar,diag%mass_diffusion_coeff_numerical_h, &
                         diag%mass_diffusion_coeff_numerical_v,diag%phase_trans_rates)
    endif
  
    ! Writing out integrals over the model domain if requested by the user.
    ! ---------------------------------------------------------------------
  
    if (lwrite_integrals) then
      if (mod(time_step_counter,2)==0) then
        call write_out_integrals(state_2%wind,state_2%rhotheta_v,diag%temperature,state_2%rho, &
                                 grid%volume,grid%inner_product_weights,grid%gravity_potential, &
                                 grid%adjacent_vector_indices_h,t_0+dtime-t_init)
      else
        call write_out_integrals(state_1%wind,state_1%rhotheta_v,diag%temperature,state_1%rho, &
                                 grid%volume,grid%inner_product_weights,grid%gravity_potential, &
                                 grid%adjacent_vector_indices_h,t_0+dtime-t_init)
      endif
    endif
  
    ! Writing the actual output.
    ! --------------------------
    
    ! interpolating to the output time
    if(t_0+dtime>=t_write .and. t_0<=t_write) then
      if (mod(time_step_counter,2) == 0) then
        new_weight = (t_write-t_0)/dtime
        old_weight = 1._wp-new_weight
        call linear_combine_two_states(state_1%rho,state_1%rhotheta_v,state_1%exner_pert,state_1%wind,state_1%temperature_soil, &
                                       state_2%rho,state_2%rhotheta_v,state_2%exner_pert,state_2%wind,state_2%temperature_soil, &
                                       state_write%rho,state_write%rhotheta_v,state_write%theta_v_pert,state_write%exner_pert, &
                                       state_write%wind,state_write%temperature_soil,old_weight,new_weight,grid%theta_v_bg)
      else
        new_weight = (t_write-t_0)/dtime
        old_weight = 1._wp-new_weight
        call linear_combine_two_states(state_2%rho,state_2%rhotheta_v,state_2%exner_pert,state_2%wind,state_2%temperature_soil, &
                                       state_1%rho,state_1%rhotheta_v,state_1%exner_pert,state_1%wind,state_1%temperature_soil, &
                                       state_write%rho,state_write%rhotheta_v,state_write%theta_v_pert,state_write%exner_pert, &
                                       state_write%wind,state_write%temperature_soil,old_weight,new_weight,grid%theta_v_bg)
      endif
    endif
  
    ! 5 minutes before the output time,the wind in the lowest layer needs to be collected for 10 m wind diag.
    if (t_0>=t_write-radius_rescale*300._wp) then
      if (wind_lowest_layer_step_counter<n_output_steps_10m_wind) then
        if (mod(time_step_counter,2)==0) then
          !$omp parallel do private(h_index)
          do h_index=1,n_vectors_h
            wind_h_lowest_layer(wind_lowest_layer_step_counter*n_vectors_h + h_index) &
            = state_1%wind(N_VECTORS - n_vectors_per_layer + h_index)
          enddo
          !$omp end parallel do
        else
          !$omp parallel do private(h_index)
          do h_index=1,n_vectors_h
            wind_h_lowest_layer(wind_lowest_layer_step_counter*n_vectors_h + h_index) &
            = state_2%wind(N_VECTORS - n_vectors_per_layer + h_index)
          enddo
          !$omp end parallel do
        endif
        wind_lowest_layer_step_counter = wind_lowest_layer_step_counter+1
      endif
    endif
  
    ! 5 minutes after the output time,the 10 m wind diag can be executed,so output can actually be written
    if(t_0+dtime>=t_write+radius_rescale*300._wp .and. t_0<=t_write+radius_rescale*300._wp) then
      ! here,output is actually written
      call write_out(diag%scalar_field_placeholder,state_write%wind,grid%latlon_interpol_indices, &
                     grid%latlon_interpol_weights,grid%exner_bg, &
                     grid%inner_product_weights,grid%volume,grid%gravity_potential,grid%from_index,grid%to_index, &
                     grid%z_vector,grid%f_vec,diag%temperature, &
                     state_write%temperature_soil,grid%area,state_write%rho,grid%z_scalar,grid%slope,grid%gravity_m, &
                     grid%adjacent_signs_h,grid%adjacent_vector_indices_h, &
                     grid%area_dual,grid%density_to_rhombi_indices,grid%density_to_rhombi_weights,state_write%exner_pert, &
                     diag%tke,t_init,t_write, &
                     grid%from_index_dual,grid%to_index_dual,diag%v_squared,grid%is_land, &
                     diag%monin_obukhov_length,diag%roughness_velocity, &
                     grid%roughness_length,grid%direction,grid%trsk_indices,grid%sfc_albedo,diag%sfc_sw_in, &
                     grid%layer_thickness,state_write%theta_v_pert, &
                     grid%theta_v_bg,grid%z_vector_dual,grid%vorticity_indices_triangles,grid%vorticity_signs_triangles, &
                     grid%trsk_weights, &
                     totally_first_step_bool,wind_h_lowest_layer,diag%rel_vort_on_triangles,diag%rel_vort, &
                     diag%pot_vort,grid%normal_distance)
      ! setting the next output time
      t_write = t_write + 60._wp*write_out_interval_min
      
      ! Calculating the speed of the model.
      !second_time = clock()
      !speed = CLOCKS_PER_SEC*60*write_out_interval_min/((double) second_time - first_time)
      !write(*,*) "Current speed: %lf\n",speed
      !first_time = clock()
      write(*,*) "Run progress:",(t_0+dtime-t_init)/3600._wp,"h"
      
      ! resetting the wind in the lowest layer to zero
      !$omp parallel do private(ji)
      do ji=1,n_output_steps_10m_wind*n_vectors_h
        wind_h_lowest_layer(ji) = 0._wp
      enddo
      !$omp end parallel do
      wind_lowest_layer_step_counter = 0
    
    endif
    
    
    
    ! This switch can be set to zero now and remains there.
    totally_first_step_bool = 0

    time_step_counter = time_step_counter+1
  
    ! giving the user information on the run progress
    write(*,*) "Step",time_step_counter,"completed."
    
    ! updating the model time
    t_0 = t_0+dtime
  
  enddo
  
  ! Clean-up.
  ! ---------
  deallocate(grid%normal_distance)
  deallocate(grid%volume)
  deallocate(grid%area)
  deallocate(grid%z_scalar)
  deallocate(grid%z_vector)
  deallocate(grid%gravity_potential)
  deallocate(grid%gravity_m)
  deallocate(grid%slope)
  deallocate(grid%theta_v_bg)
  deallocate(grid%exner_bg)
  deallocate(grid%exner_bg_grad)
  deallocate(grid%layer_thickness)
  deallocate(grid%area_dual)
  deallocate(grid%z_vector_dual)
  deallocate(grid%normal_distance_dual)
  deallocate(grid%vorticity_indices_triangles)
  deallocate(grid%vorticity_signs_triangles)
  deallocate(grid%f_vec)
  deallocate(grid%trsk_indices)
  deallocate(grid%trsk_modified_curl_indices)
  deallocate(grid%from_index)
  deallocate(grid%to_index)
  deallocate(grid%from_index_dual)
  deallocate(grid%to_index_dual)
  deallocate(grid%adjacent_vector_indices_h)
  deallocate(grid%adjacent_signs_h)
  deallocate(grid%density_to_rhombi_indices)
  deallocate(grid%latitude_scalar)
  deallocate(grid%longitude_scalar)
  deallocate(grid%inner_product_weights)
  deallocate(grid%direction)
  deallocate(grid%density_to_rhombi_weights)
  deallocate(grid%trsk_weights)
  deallocate(grid%sfc_albedo)
  deallocate(grid%sfc_rho_c)
  deallocate(grid%t_conduc_soil)
  deallocate(grid%roughness_length)
  deallocate(grid%is_land)
  deallocate(grid%latlon_interpol_indices)
  deallocate(grid%latlon_interpol_weights)
  deallocate(grid%z_soil_interface)
  deallocate(grid%z_soil_center)
  deallocate(grid%t_const_soil)
  deallocate(state_1%rho)
  deallocate(state_1%rhotheta_v)
  deallocate(state_1%theta_v_pert)
  deallocate(state_1%exner_pert)
  deallocate(state_1%wind)
  deallocate(state_1%temperature_soil)
  deallocate(state_2%rho)
  deallocate(state_2%rhotheta_v)
  deallocate(state_2%theta_v_pert)
  deallocate(state_2%exner_pert)
  deallocate(state_2%wind)
  deallocate(state_2%temperature_soil)
  deallocate(state_tendency%rho)
  deallocate(state_tendency%rhotheta_v)
  deallocate(state_tendency%theta_v_pert)
  deallocate(state_tendency%exner_pert)
  deallocate(state_tendency%wind)
  deallocate(state_tendency%temperature_soil)
  deallocate(diag%flux_density)
  deallocate(diag%flux_density_div)
  deallocate(diag%rel_vort_on_triangles)
  deallocate(diag%rel_vort)
  deallocate(diag%pot_vort)
  deallocate(diag%temperature)
  deallocate(diag%c_g_p_field)
  deallocate(diag%v_squared)
  deallocate(diag%wind_div)
  deallocate(diag%curl_of_vorticity)
  deallocate(diag%scalar_field_placeholder)
  deallocate(diag%vector_field_placeholder)
  deallocate(diag%u_at_edge)
  deallocate(diag%v_at_edge)
  deallocate(diag%u_at_cell)
  deallocate(diag%v_at_cell)
  deallocate(diag%n_squared)
  deallocate(diag%dv_hdz)
  deallocate(diag%scalar_flux_resistance)
  deallocate(diag%power_flux_density_sensible)
  deallocate(diag%power_flux_density_latent)
  deallocate(diag%roughness_velocity)
  deallocate(diag%monin_obukhov_length)
  deallocate(diag%temperature_diffusion_heating)
  deallocate(diag%friction_acc)
  deallocate(diag%heating_diss)
  deallocate(diag%molecular_diffusion_coeff)
  deallocate(diag%mass_diffusion_coeff_numerical_h)
  deallocate(diag%mass_diffusion_coeff_numerical_v)
  deallocate(diag%temp_diffusion_coeff_numerical_h)
  deallocate(diag%temp_diffusion_coeff_numerical_v)
  deallocate(diag%pressure_gradient_decel_factor)
  deallocate(diag%condensates_sediment_heat)
  deallocate(diag%mass_diff_tendency)
  deallocate(diag%phase_trans_rates)
  deallocate(diag%phase_trans_heating_rate)
  deallocate(diag%viscosity)
  deallocate(diag%viscosity_rhombi)
  deallocate(diag%viscosity_triangles)
  deallocate(diag%vert_hor_viscosity)
  deallocate(diag%tke)
  deallocate(diag%pgrad_acc_old)
  deallocate(diag%pressure_gradient_acc_neg_nl)
  deallocate(diag%pressure_gradient_acc_neg_l)
  deallocate(diag%pressure_grad_condensates_v)
  deallocate(diag%v_squared_grad)
  deallocate(diag%pot_vort_tend)
  deallocate(diag%sfc_sw_in)
  deallocate(diag%sfc_lw_out)
  deallocate(diag%radiation_tendency)
  deallocate(state_write%rho)
  deallocate(state_write%rhotheta_v)
  deallocate(state_write%theta_v_pert)
  deallocate(state_write%exner_pert)
  deallocate(state_write%wind)
  deallocate(state_write%temperature_soil)
  write(*,*) "**"
  ! clock_t end = clock()
  ! speed = CLOCKS_PER_SEC*(60*total_run_span_min + 300)/((double) end - begin)
  ! write(*,*) "Average speed: %lf",speed
  write(*,*) "GAME over."

end program control










