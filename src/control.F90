! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

program control

  ! This program organizes the model, manages the time stepping, calls model output, collects the lowest model layer wind for 10 m wind mean and so on.
  ! All the memory needed for the integration is allocated and freed here.

  use mo_definitions,            only: wp,t_grid,t_state,t_diag
  use mo_grid_nml,               only: n_layers,n_cells,n_edges,n_triangles,n_lat_io_points,n_lon_io_points, &
                                       n_levels,grid_nml_setup
  use mo_constituents_nml,       only: n_constituents,n_condensed_constituents,constituents_nml_setup
  use mo_run_nml,                only: run_span_min,run_nml_setup,t_init,luse_bg_state
  use mo_grid_setup,             only: eff_hor_res,radius_rescale,set_grid_properties,set_background_state,dtime,toa
  use mo_io_nml,                 only: n_output_steps_10m_wind,lwrite_integrals,write_out_interval_min, &
                                       io_nml_setup,ideal_input_id
  use mo_surface_nml,            only: nsoillays,surface_nml_setup
  use mo_rad_nml,                only: rad_config,radiation_dtime,rad_nml_setup
  use mo_diff_nml,               only: diff_nml_setup
  use mo_set_initial_state,      only: set_ideal_init,read_init_data
  use mo_derived,                only: temperature_diagnostics
  use mo_write_output,           only: write_out,write_out_integrals
  use mo_manage_radiation_calls, only: update_rad_fluxes
  use mo_rrtmgp_coupler,         only: radiation_init
  use mo_manage_pchevi,          only: manage_pchevi
  use mo_linear_combination,     only: linear_combine_two_states
  use mo_gradient_operators,     only: grad_hor_cov,grad_hor,grad_vert
  use mo_inner_product,          only: inner_product
  use omp_lib,                   only: omp_get_num_threads

  implicit none
  
  type(t_state)         :: state_1                        ! states (prognostic variables)
  type(t_state)         :: state_2                        ! states (prognostic variables)
  type(t_state)         :: state_write                    ! states (prognostic variables) to be written out
  type(t_state)         :: state_tend                     ! state containing the explicit tendencies of the prognostic variables
  type(t_diag)          :: diag                           ! diagnostic quantities
  type(t_grid)          :: grid                           ! grid quantities
  logical               :: ltotally_first_step            ! boolean that is only true at the very first step
  logical               :: lrad_update                    ! boolean indicating that radiative fluxes need to be updated
  integer               :: omp_num_threads                ! number of OMP threads
  integer               :: ji                             ! horizontal index
  integer               :: time_step_counter              ! total time step counter
  integer               :: time_step_10_m_wind            ! time step of the 10 m wind averaging
  real(wp)              :: t_0                            ! Unix timestamp of the old time step
  real(wp)              :: t_write                        ! Unix timestamp of the next output time
  real(wp)              :: t_rad_update                   ! Unix timestamp of the next radiation update
  real(wp)              :: new_weight                     ! time interpolation weight for computing the output state
  real(wp)              :: old_weight                     ! time interpolation weight for computing the output state
  real(wp)              :: max_speed_hor                  ! maximum horizontal wind speed used for computing the horizontal advective Courant number
  real(wp)              :: max_speed_ver                  ! maximum vertical wind speed used for computing the vertical advective Courant number
  real(wp)              :: normal_dist_min_hor            ! minimum horizontal normal distance
  real(wp)              :: normal_dist_min_ver            ! minimum vertical normal distance
  real(wp)              :: init_timestamp                 ! CPU time of the model start
  real(wp)              :: begin_timestamp                ! CPU time of the beginning of a runtime measurement interval
  real(wp)              :: end_timestamp                  ! CPU time of the end of a runtime measurement interval
  real(wp), allocatable :: wind_h_lowest_layer(:,:)       ! horizontal wind in the lowest layer (collected for computing the 10-m wind)
  character(len=82)     :: stars                          ! string containing stars for console output

  ! taking the timestamp to measure the performance
  call cpu_time(init_timestamp)
  begin_timestamp = init_timestamp
  
  ! console output
  stars = "**********************************************************************************"
  write(*,*) stars
  write(*,*) "*                                                                                *"
  write(*,*) "*                                  This is GAME                                  *"
  write(*,*) "*                      Geophysical Fluids Modeling Framework                     *"
  write(*,*) "*                                                                                *"
  write(*,*) "*                         Released under the MIT license.                        *"
  write(*,*) "*          Visit https://github.com/OpenNWP/GAME for more information.           *"
  write(*,*) "*                                                                                *"
  write(*,*) stars
  
  ! setting up the namelists (model configuration)
  call grid_nml_setup()
  call run_nml_setup()
  call constituents_nml_setup()
  call diff_nml_setup()
  call surface_nml_setup()
  
  ! Allocating memory
  ! -----------------
  
  allocate(grid%dx(n_edges,n_layers))
  allocate(grid%dz(n_cells,n_levels))
  allocate(grid%volume(n_cells,n_layers))
  allocate(grid%area_h(n_edges,n_layers))
  allocate(grid%area_v(n_cells,n_levels))
  allocate(grid%z_scalar(n_cells,n_layers))
  allocate(grid%z_vector_h(n_edges,n_layers))
  allocate(grid%z_vector_v(n_cells,n_levels))
  allocate(grid%gravity_potential(n_cells,n_layers))
  allocate(grid%gravity_m_v(n_cells,n_levels))
  allocate(grid%slope(n_edges,n_layers))
  allocate(grid%theta_v_bg(n_cells,n_layers))
  allocate(grid%exner_bg(n_cells,n_layers))
  allocate(grid%exner_bg_grad_h(n_edges,n_layers))
  allocate(grid%exner_bg_grad_v(n_cells,n_levels))
  allocate(grid%layer_thickness(n_cells,n_layers))
  allocate(grid%area_dual_h(n_edges,n_levels))
  allocate(grid%area_dual_v(n_triangles,n_layers))
  allocate(grid%z_vector_dual_h(n_edges,n_levels))
  allocate(grid%z_vector_dual_v(n_triangles,n_layers))
  allocate(grid%dy(n_edges,n_levels))
  allocate(grid%dz_dual(n_triangles,n_layers))
  allocate(grid%vorticity_indices_triangles(3,n_triangles))
  allocate(grid%vorticity_signs_triangles(3,n_triangles))
  allocate(grid%f_vec_h(n_edges))
  allocate(grid%f_vec_v(n_edges))
  allocate(grid%trsk_indices(10,n_edges))
  allocate(grid%trsk_modified_curl_indices(10,n_edges))
  allocate(grid%from_cell(n_edges))
  allocate(grid%to_cell(n_edges))
  allocate(grid%from_cell_dual(n_edges))
  allocate(grid%to_cell_dual(n_edges))
  allocate(grid%adjacent_edges(6,n_cells))
  allocate(grid%adjacent_signs(6,n_cells))
  allocate(grid%density_to_rhombi_indices(4,n_edges))
  allocate(grid%lat_c(n_cells))
  allocate(grid%lon_c(n_cells))
  allocate(grid%inner_product_weights(8,n_cells,n_layers))
  allocate(grid%direction(n_edges))
  allocate(grid%density_to_rhombi_weights(4,n_edges))
  allocate(grid%trsk_weights(10,n_edges))
  allocate(grid%sfc_albedo(n_cells))
  allocate(grid%sfc_rho_c(n_cells))
  allocate(grid%t_conduc_soil(n_cells))
  allocate(grid%roughness_length(n_cells))
  allocate(grid%land_fraction(n_cells))
  allocate(grid%lake_fraction(n_cells))
  allocate(grid%latlon_interpol_indices(5,n_lat_io_points,n_lon_io_points))
  allocate(grid%latlon_interpol_weights(5,n_lat_io_points,n_lon_io_points))
  allocate(grid%lat_output_vector(n_lat_io_points))
  allocate(grid%lon_output_vector(n_lon_io_points))
  allocate(grid%z_soil_interface(nsoillays+1))
  allocate(grid%z_soil_center(nsoillays))
  allocate(grid%t_const_soil(n_cells))
  allocate(state_1%rho(n_cells,n_layers,n_constituents))
  allocate(state_1%rhotheta_v(n_cells,n_layers))
  allocate(state_1%theta_v_pert(n_cells,n_layers))
  allocate(state_1%exner_pert(n_cells,n_layers))
  allocate(state_1%wind_h(n_edges,n_layers))
  allocate(state_1%wind_v(n_cells,n_levels))
  allocate(state_1%temperature_soil(n_cells,nsoillays))
  allocate(state_2%rho(n_cells,n_layers,n_constituents))
  allocate(state_2%rhotheta_v(n_cells,n_layers))
  allocate(state_2%theta_v_pert(n_cells,n_layers))
  allocate(state_2%exner_pert(n_cells,n_layers))
  allocate(state_2%wind_h(n_edges,n_layers))
  allocate(state_2%wind_v(n_cells,n_levels))
  allocate(state_2%temperature_soil(n_cells,nsoillays))
  allocate(state_tend%rho(n_cells,n_layers,n_constituents))
  allocate(state_tend%rhotheta_v(n_cells,n_layers))
  allocate(state_tend%theta_v_pert(n_cells,n_layers))
  allocate(state_tend%exner_pert(n_cells,n_layers))
  allocate(state_tend%wind_h(n_edges,n_layers))
  allocate(state_tend%wind_v(n_cells,n_levels))
  allocate(state_tend%temperature_soil(n_cells,nsoillays))
  allocate(diag%flux_density_h(n_edges,n_layers))
  allocate(diag%flux_density_v(n_cells,n_levels))
  allocate(diag%flux_density_div(n_cells,n_layers))
  allocate(diag%zeta_on_triangles(n_triangles,n_layers))
  allocate(diag%zeta_h(n_edges,n_levels))
  allocate(diag%zeta_v(n_edges,n_layers))
  allocate(diag%eta_h(n_edges,n_levels))
  allocate(diag%eta_v(n_edges,n_layers))
  allocate(diag%temperature(n_cells,n_layers))
  allocate(diag%v_squared(n_cells,n_layers))
  allocate(diag%wind_div(n_cells,n_layers))
  allocate(diag%curl_of_vorticity_h(n_edges,n_layers))
  allocate(diag%scalar_placeholder(n_cells,n_layers))
  allocate(diag%vector_placeholder_h(n_edges,n_layers))
  allocate(diag%vector_placeholder_v(n_cells,n_levels))
  allocate(diag%n_squared(n_cells,n_layers))
  allocate(diag%dv_hdz(n_edges,n_levels))
  allocate(diag%scalar_flux_resistance(n_cells))
  allocate(diag%power_flux_density_sens_sea(n_cells))
  allocate(diag%power_flux_density_sens_soil(n_cells))
  allocate(diag%power_flux_density_lat_sea(n_cells))
  allocate(diag%power_flux_density_lat_lake(n_cells))
  allocate(diag%roughness_length(n_cells))
  allocate(diag%roughness_velocity(n_cells))
  allocate(diag%monin_obukhov_length(n_cells))
  allocate(diag%temperature_diff_heating(n_cells,n_layers))
  allocate(diag%friction_acc_h(n_edges,n_layers))
  allocate(diag%friction_acc_v(n_cells,n_levels))
  allocate(diag%heating_diss(n_cells,n_layers))
  allocate(diag%molecular_diff_coeff(n_cells,n_layers))
  allocate(diag%mass_diff_coeff_eff_h(n_cells,n_layers))
  allocate(diag%mass_diff_coeff_eff_v(n_cells,n_layers))
  allocate(diag%temp_diff_coeff_eff_h(n_cells,n_layers))
  allocate(diag%temp_diff_coeff_eff_v(n_cells,n_layers))
  allocate(diag%p_grad_decel_factor(n_cells,n_layers))
  allocate(diag%condensates_sediment_heat(n_cells,n_layers))
  allocate(diag%mass_diff_tendency(n_cells,n_layers,n_constituents))
  allocate(diag%phase_trans_rates(n_cells,n_layers,n_condensed_constituents+1))
  allocate(diag%phase_trans_heating_rate(n_cells,n_layers))
  allocate(diag%viscosity(n_cells,n_layers))
  allocate(diag%viscosity_rhombi(n_edges,n_layers))
  allocate(diag%viscosity_triangles(n_triangles,n_layers))
  allocate(diag%vert_hor_viscosity(n_edges,n_levels))
  allocate(diag%tke(n_cells,n_layers))
  allocate(diag%sst(n_cells))
  allocate(diag%p_grad_acc_old_h(n_edges,n_layers))
  allocate(diag%p_grad_acc_neg_nl_h(n_edges,n_layers))
  allocate(diag%p_grad_acc_neg_nl_v(n_cells,n_levels))
  allocate(diag%p_grad_acc_neg_l_h(n_edges,n_layers))
  allocate(diag%p_grad_acc_neg_l_v(n_cells,n_levels))
  allocate(diag%p_grad_condensates_v(n_cells,n_levels))
  allocate(diag%v_squared_grad_h(n_edges,n_layers))
  allocate(diag%v_squared_grad_v(n_cells,n_levels))
  allocate(diag%pot_vort_tend_h(n_edges,n_layers))
  allocate(diag%pot_vort_tend_v(n_cells,n_levels))
  allocate(diag%sfc_sw_in(n_cells))
  allocate(diag%sfc_lw_out(n_cells))
  allocate(diag%radiation_tendency(n_cells,n_layers))
  allocate(diag%a_rain(n_cells,n_layers))
  allocate(state_write%rho(n_cells,n_layers,n_constituents))
  allocate(state_write%rhotheta_v(n_cells,n_layers))
  allocate(state_write%theta_v_pert(n_cells,n_layers))
  allocate(state_write%exner_pert(n_cells,n_layers))
  allocate(state_write%wind_h(n_edges,n_layers))
  allocate(state_write%wind_v(n_cells,n_levels))
  allocate(state_write%temperature_soil(n_cells,nsoillays))
  
  ! setting all arrays to zero
  !$omp parallel workshare
  grid%dx = 0._wp
  grid%dz = 0._wp
  grid%volume = 0._wp
  grid%area_h = 0._wp
  grid%area_v = 0._wp
  grid%z_scalar = 0._wp
  grid%z_vector_h = 0._wp
  grid%z_vector_v = 0._wp
  grid%gravity_potential = 0._wp
  grid%gravity_m_v = 0._wp
  grid%slope = 0._wp
  grid%theta_v_bg = 0._wp
  grid%exner_bg = 0._wp
  grid%exner_bg_grad_h = 0._wp
  grid%exner_bg_grad_v = 0._wp
  grid%layer_thickness = 0._wp
  grid%area_dual_h = 0._wp
  grid%area_dual_v = 0._wp
  grid%z_vector_dual_h = 0._wp
  grid%z_vector_dual_v = 0._wp
  grid%dy = 0._wp
  grid%dz_dual = 0._wp
  grid%vorticity_indices_triangles = 0
  grid%vorticity_signs_triangles = 0
  grid%f_vec_h = 0._wp
  grid%f_vec_v = 0._wp
  grid%trsk_indices = 0
  grid%trsk_modified_curl_indices = 0
  grid%from_cell = 0
  grid%to_cell = 0
  grid%from_cell_dual = 0
  grid%to_cell_dual = 0
  grid%adjacent_edges = 0
  grid%adjacent_signs = 0
  grid%density_to_rhombi_indices = 0
  grid%lat_c = 0._wp
  grid%lon_c = 0._wp
  grid%inner_product_weights = 0._wp
  grid%direction = 0._wp
  grid%density_to_rhombi_weights = 0._wp
  grid%trsk_weights = 0._wp
  grid%sfc_albedo = 0._wp
  grid%sfc_rho_c = 0._wp
  grid%t_conduc_soil = 0._wp
  grid%roughness_length = 0._wp
  grid%land_fraction = 0._wp
  grid%lake_fraction = 0._wp
  grid%latlon_interpol_indices = 0
  grid%latlon_interpol_weights = 0._wp
  grid%lat_output_vector = 0._wp
  grid%lon_output_vector = 0._wp
  grid%z_soil_interface = 0._wp
  grid%z_soil_center = 0._wp
  grid%t_const_soil = 0._wp
  state_1%rho = 0._wp
  state_1%rhotheta_v = 0._wp
  state_1%theta_v_pert = 0._wp
  state_1%exner_pert = 0._wp
  state_1%wind_h = 0._wp
  state_1%wind_v = 0._wp
  state_1%temperature_soil = 0._wp
  state_2%rho = 0._wp
  state_2%rhotheta_v = 0._wp
  state_2%theta_v_pert = 0._wp
  state_2%exner_pert = 0._wp
  state_2%wind_h = 0._wp
  state_2%wind_v = 0._wp
  state_2%temperature_soil = 0._wp
  state_tend%rho = 0._wp
  state_tend%rhotheta_v = 0._wp
  state_tend%theta_v_pert = 0._wp
  state_tend%exner_pert = 0._wp
  state_tend%wind_h = 0._wp
  state_tend%wind_v = 0._wp
  state_tend%temperature_soil = 0._wp
  diag%flux_density_h = 0._wp
  diag%flux_density_v = 0._wp
  diag%flux_density_div = 0._wp
  diag%zeta_on_triangles = 0._wp
  diag%zeta_h = 0._wp
  diag%zeta_v = 0._wp
  diag%eta_h = 0._wp
  diag%eta_v = 0._wp
  diag%temperature = 0._wp
  diag%v_squared = 0._wp
  diag%wind_div = 0._wp
  diag%curl_of_vorticity_h = 0._wp
  diag%scalar_placeholder = 0._wp
  diag%vector_placeholder_h = 0._wp
  diag%vector_placeholder_v = 0._wp
  diag%n_squared = 0._wp
  diag%dv_hdz = 0._wp
  diag%scalar_flux_resistance = 0._wp
  diag%power_flux_density_sens_sea = 0._wp
  diag%power_flux_density_sens_soil = 0._wp
  diag%power_flux_density_lat_sea = 0._wp
  diag%power_flux_density_lat_lake = 0._wp
  diag%roughness_length = 0._wp
  diag%roughness_velocity = 0._wp
  diag%monin_obukhov_length = 0._wp
  diag%temperature_diff_heating = 0._wp
  diag%friction_acc_h = 0._wp
  diag%friction_acc_v = 0._wp
  diag%heating_diss = 0._wp
  diag%molecular_diff_coeff = 0._wp
  diag%mass_diff_coeff_eff_h = 0._wp
  diag%mass_diff_coeff_eff_v = 0._wp
  diag%temp_diff_coeff_eff_h = 0._wp
  diag%temp_diff_coeff_eff_v = 0._wp
  diag%p_grad_decel_factor = 0._wp
  diag%condensates_sediment_heat = 0._wp
  diag%mass_diff_tendency = 0._wp
  diag%phase_trans_rates = 0._wp
  diag%phase_trans_heating_rate = 0._wp
  diag%viscosity = 0._wp
  diag%viscosity_rhombi = 0._wp
  diag%viscosity_triangles = 0._wp
  diag%vert_hor_viscosity = 0._wp
  diag%tke = 0._wp
  diag%sst = 0._wp
  diag%p_grad_acc_old_h = 0._wp
  diag%p_grad_acc_neg_nl_h = 0._wp
  diag%p_grad_acc_neg_nl_v = 0._wp
  diag%p_grad_acc_neg_l_h = 0._wp
  diag%p_grad_acc_neg_l_v = 0._wp
  diag%p_grad_condensates_v = 0._wp
  diag%v_squared_grad_h = 0._wp
  diag%v_squared_grad_v = 0._wp
  diag%pot_vort_tend_h = 0._wp
  diag%pot_vort_tend_v = 0._wp
  diag%sfc_sw_in = 0._wp
  diag%sfc_lw_out = 0._wp
  diag%radiation_tendency = 0._wp
  diag%a_rain = 0._wp
  state_write%rho = 0._wp
  state_write%rhotheta_v = 0._wp
  state_write%theta_v_pert = 0._wp
  state_write%exner_pert = 0._wp
  state_write%wind_h = 0._wp
  state_write%wind_v = 0._wp
  state_write%temperature_soil = 0._wp
  !$omp end parallel workshare
  
  ! reading the grid
  write(*,*) "Reading grid data ..."
  call set_grid_properties(grid)
  call set_background_state(grid)
  
  ! initializing the diagnostic roughness length
  !$omp parallel workshare
  diag%roughness_length = grid%roughness_length
  !$omp end parallel workshare
  
  call io_nml_setup()
  call rad_nml_setup(grid)
  
  ! some gradients have to be computed here
  call grad_hor_cov(grid%z_scalar,grid%slope,grid)
  call grad_vert(grid%gravity_potential,grid%gravity_m_v,grid)
  ! if we do not use a hydrostatic background state we write the negative acceleration due to gravity into the
  ! negative "linear" pressure gradient acceleration
  if (.not. luse_bg_state) then
    call grad_hor_cov(grid%gravity_potential,diag%p_grad_acc_neg_l_h,grid)
    !$omp parallel workshare
    diag%p_grad_acc_neg_l_h = diag%p_grad_acc_neg_l_h
    diag%p_grad_acc_neg_l_v = grid%gravity_m_v
    !$omp end parallel workshare
  endif
  call grad_vert(grid%exner_bg,grid%exner_bg_grad_v,grid)
  call grad_hor(grid%exner_bg,grid%exner_bg_grad_h,grid%exner_bg_grad_v,grid)
  write(*,*) "Grid loaded successfully."
  
  ! rescaling times for small Earth experiments
  run_span_min = radius_rescale*run_span_min
  write_out_interval_min = radius_rescale*write_out_interval_min
  
  ! Reading and processing user input finished.
  
  write(*,*) "Setting initial state ..."
  
  ! ideal test case
  if (ideal_input_id/=-1) then
    call set_ideal_init(state_1,diag,grid)
  ! NWP mode
  else
    call read_init_data(state_1,diag,grid)
  endif
  
  write(*,*) "Initial state set successfully."
  write(*,*) stars
  
  if (rad_config>0) then
    write(*,fmt="(A,F10.3,A2)") " Radiation time step:",radiation_dtime,"s"
  endif
  
  ! finding the minimum horizontal grid distance
  normal_dist_min_hor = eff_hor_res
  do ji=1,n_edges
    if (grid%dx(ji,n_layers)<normal_dist_min_hor) then
      normal_dist_min_hor = grid%dx(ji,n_layers)
    endif
  enddo
  ! finding the minimum vertical grid distance
  normal_dist_min_ver = toa/n_layers
  do ji=1,n_cells
    if (grid%dz(ji,n_layers)<normal_dist_min_ver) then
      normal_dist_min_ver = grid%dz(ji,n_layers)
    endif
  enddo
  
  write(*,fmt="(A,F8.3,A3)") " Effective horizontal resolution:",1e-3*eff_hor_res,"km"
  write(*,fmt="(A,F8.3,A3)") " Minimum horizontal normal distance:",1e-3*normal_dist_min_hor," km"
  max_speed_hor = 100._wp
  write(*,fmt="(A,F6.3)") " Horizontal advective Courant number:",dtime/normal_dist_min_hor*max_speed_hor
  max_speed_ver = 0.1_wp
  write(*,fmt="(A,F6.3)") " Vertical advective Courant number:",dtime/normal_dist_min_ver*max_speed_ver
  write(*,*) stars
  write(*,*) "It begins."
  write(*,*) stars
  
  allocate(wind_h_lowest_layer(n_edges,n_output_steps_10m_wind))
  ! here,for all output time steps,the initial value is used
  do time_step_10_m_wind=1,n_output_steps_10m_wind
    !$omp parallel workshare
    wind_h_lowest_layer(:,time_step_10_m_wind) = state_1%wind_h(:,n_layers)
    !$omp end parallel workshare
  enddo
  
  ! getting the number of OMP threads
  !$omp parallel
  omp_num_threads = omp_get_num_threads()
  !$omp end parallel
  
  ! time coordinate of the old RK step
  t_0 = t_init
  ! for radiation and writing output it is necessary that the temperature is set correctly
  call temperature_diagnostics(state_1,diag,grid)
  ! configuring radiation and calculating radiative fluxes for the first time
  t_rad_update = t_0
  if (rad_config>0) then
    if (rad_config==1) then
      call radiation_init()
    endif
    call update_rad_fluxes(state_1,diag,grid,t_0)
    t_rad_update = t_rad_update+radiation_dtime
  endif
  
  ! This is necessary because at the very first step of the model integration,some things are handled differently
  ! in the time stepping and in writing the output.
  ltotally_first_step = .true.
  ! writing out the initial state of the model run
  t_write = t_init
  call inner_product(state_1%wind_h,state_1%wind_v,state_1%wind_h,state_1%wind_v,diag%v_squared,grid)
  call write_out(state_1,diag,grid,wind_h_lowest_layer,t_init,t_write,ltotally_first_step)
  t_write = t_0 + 60._wp*write_out_interval_min
  
  write(*,"(A,F10.3,A2)") " Run progress:",(t_init-t_init)/3600._wp,"h"
  
  if (lwrite_integrals) then
    call write_out_integrals(state_1,diag,grid,0._wp)
  endif
  
  ! Preparation of the actual integration.
  ! --------------------------------------
  time_step_10_m_wind = 1
  !$omp parallel workshare
  state_2 = state_1
  !$omp end parallel workshare
  
  ! This is the loop over the time steps.
  ! -------------------------------------
  ! this is to store the speed of the model integration
  time_step_counter = 0
  do while (t_0<t_init+60._wp*run_span_min+radius_rescale*300._wp)
    
    ! Checking if the radiative fluxes need to be updated:
    ! ---------------------------------------------------
    
    if (t_0<=t_rad_update .and. t_0+dtime>=t_rad_update .and. .not.ltotally_first_step) then
      lrad_update = .true.
      t_rad_update = t_rad_update+radiation_dtime
    else
      lrad_update = .false.
    endif
    
    ! time step integration
    if (mod(time_step_counter,2)==0) then
      call manage_pchevi(state_1,state_2,state_tend,diag,grid,lrad_update,ltotally_first_step,t_0)
    else
      call manage_pchevi(state_2,state_1,state_tend,diag,grid,lrad_update,ltotally_first_step,t_0)
    endif
  
    ! Writing out integrals over the model domain if requested by the user.
    ! ---------------------------------------------------------------------
  
    if (lwrite_integrals) then
      if (mod(time_step_counter,2)==0) then
        call write_out_integrals(state_2,diag,grid,t_0+dtime-t_init)
      else
        call write_out_integrals(state_1,diag,grid,t_0+dtime-t_init)
      endif
    endif
  
    ! Writing the actual output.
    ! --------------------------
    
    ! interpolating to the output time
    if(t_0+dtime>=t_write .and. t_0<=t_write) then
      new_weight = (t_write-t_0)/dtime
      old_weight = 1._wp-new_weight
      if (mod(time_step_counter,2)==0) then
        call linear_combine_two_states(state_1,state_2,state_write,old_weight,new_weight,grid)
      else
        call linear_combine_two_states(state_2,state_1,state_write,old_weight,new_weight,grid)
      endif
    endif
  
    ! 5 minutes before the output time,the wind in the lowest layer needs to be collected for 10 m wind diag.
    if (t_0>=t_write-radius_rescale*300._wp) then
      if (time_step_10_m_wind<=n_output_steps_10m_wind) then
        if (mod(time_step_counter,2)==0) then
          !$omp parallel workshare
          wind_h_lowest_layer(:,time_step_10_m_wind) = state_1%wind_h(:,n_layers)
          !$omp end parallel workshare
        else
          !$omp parallel workshare
          wind_h_lowest_layer(:,time_step_10_m_wind) = state_2%wind_h(:,n_layers)
          !$omp end parallel workshare
        endif
        time_step_10_m_wind = time_step_10_m_wind+1
      endif
    endif
    
    ! 5 minutes after the output time,the 10 m wind diag can be executed,so output can actually be written
    if(t_0+dtime>=t_write+radius_rescale*300._wp .and. t_0<=t_write+radius_rescale*300._wp) then
      ! here,output is actually written
      call write_out(state_write,diag,grid,wind_h_lowest_layer,t_init,t_write,ltotally_first_step)
      ! setting the next output time
      t_write = t_write + 60._wp*write_out_interval_min
      
      ! calculating the speed of the model
      call cpu_time(end_timestamp)
      write(*,fmt="(A,F9.3)") " Current speed:",60._wp*write_out_interval_min/((end_timestamp-begin_timestamp)/omp_num_threads)
      call cpu_time(begin_timestamp)
      write(*,fmt="(A,F10.3,A2)") " Run progress:",(t_0+dtime-t_init)/3600._wp,"h"
      
      ! resetting the wind in the lowest layer to zero
      !$omp parallel workshare
      wind_h_lowest_layer = 0._wp
      !$omp end parallel workshare
      time_step_10_m_wind = 1
    
    endif
    
    ! This switch can be set to false now and remains there.
    ltotally_first_step = .false.

    time_step_counter = time_step_counter+1
  
    ! giving the user information on the run progress
    write(*,*) "Step",time_step_counter,"completed."
    
    ! updating the model time
    t_0 = t_0+dtime
  
  enddo
  
  ! Clean-up.
  ! ---------
  ! freeing the memory
  deallocate(wind_h_lowest_layer)
  deallocate(grid%dx)
  deallocate(grid%dz)
  deallocate(grid%volume)
  deallocate(grid%area_h)
  deallocate(grid%area_v)
  deallocate(grid%z_scalar)
  deallocate(grid%z_vector_h)
  deallocate(grid%z_vector_v)
  deallocate(grid%gravity_potential)
  deallocate(grid%gravity_m_v)
  deallocate(grid%slope)
  deallocate(grid%theta_v_bg)
  deallocate(grid%exner_bg)
  deallocate(grid%exner_bg_grad_h)
  deallocate(grid%exner_bg_grad_v)
  deallocate(grid%layer_thickness)
  deallocate(grid%area_dual_h)
  deallocate(grid%area_dual_v)
  deallocate(grid%z_vector_dual_h)
  deallocate(grid%z_vector_dual_v)
  deallocate(grid%dy)
  deallocate(grid%dz_dual)
  deallocate(grid%vorticity_indices_triangles)
  deallocate(grid%vorticity_signs_triangles)
  deallocate(grid%f_vec_h)
  deallocate(grid%f_vec_v)
  deallocate(grid%trsk_indices)
  deallocate(grid%trsk_modified_curl_indices)
  deallocate(grid%from_cell)
  deallocate(grid%to_cell)
  deallocate(grid%from_cell_dual)
  deallocate(grid%to_cell_dual)
  deallocate(grid%adjacent_edges)
  deallocate(grid%adjacent_signs)
  deallocate(grid%density_to_rhombi_indices)
  deallocate(grid%lat_c)
  deallocate(grid%lon_c)
  deallocate(grid%inner_product_weights)
  deallocate(grid%direction)
  deallocate(grid%density_to_rhombi_weights)
  deallocate(grid%trsk_weights)
  deallocate(grid%sfc_albedo)
  deallocate(grid%sfc_rho_c)
  deallocate(grid%t_conduc_soil)
  deallocate(grid%roughness_length)
  deallocate(grid%land_fraction)
  deallocate(grid%lake_fraction)
  deallocate(grid%latlon_interpol_indices)
  deallocate(grid%latlon_interpol_weights)
  deallocate(grid%lat_output_vector)
  deallocate(grid%lon_output_vector)
  deallocate(grid%z_soil_interface)
  deallocate(grid%z_soil_center)
  deallocate(grid%t_const_soil)
  deallocate(state_1%rho)
  deallocate(state_1%rhotheta_v)
  deallocate(state_1%theta_v_pert)
  deallocate(state_1%exner_pert)
  deallocate(state_1%wind_h)
  deallocate(state_1%wind_v)
  deallocate(state_1%temperature_soil)
  deallocate(state_2%rho)
  deallocate(state_2%rhotheta_v)
  deallocate(state_2%theta_v_pert)
  deallocate(state_2%exner_pert)
  deallocate(state_2%wind_h)
  deallocate(state_2%wind_v)
  deallocate(state_2%temperature_soil)
  deallocate(state_tend%rho)
  deallocate(state_tend%rhotheta_v)
  deallocate(state_tend%theta_v_pert)
  deallocate(state_tend%exner_pert)
  deallocate(state_tend%wind_h)
  deallocate(state_tend%wind_v)
  deallocate(state_tend%temperature_soil)
  deallocate(diag%flux_density_h)
  deallocate(diag%flux_density_v)
  deallocate(diag%flux_density_div)
  deallocate(diag%zeta_on_triangles)
  deallocate(diag%zeta_h)
  deallocate(diag%zeta_v)
  deallocate(diag%eta_h)
  deallocate(diag%eta_v)
  deallocate(diag%temperature)
  deallocate(diag%v_squared)
  deallocate(diag%wind_div)
  deallocate(diag%curl_of_vorticity_h)
  deallocate(diag%scalar_placeholder)
  deallocate(diag%vector_placeholder_h)
  deallocate(diag%vector_placeholder_v)
  deallocate(diag%n_squared)
  deallocate(diag%dv_hdz)
  deallocate(diag%scalar_flux_resistance)
  deallocate(diag%power_flux_density_sens_sea)
  deallocate(diag%power_flux_density_sens_soil)
  deallocate(diag%power_flux_density_lat_sea)
  deallocate(diag%power_flux_density_lat_lake)
  deallocate(diag%roughness_length)
  deallocate(diag%roughness_velocity)
  deallocate(diag%monin_obukhov_length)
  deallocate(diag%temperature_diff_heating)
  deallocate(diag%friction_acc_h)
  deallocate(diag%friction_acc_v)
  deallocate(diag%heating_diss)
  deallocate(diag%molecular_diff_coeff)
  deallocate(diag%mass_diff_coeff_eff_h)
  deallocate(diag%mass_diff_coeff_eff_v)
  deallocate(diag%temp_diff_coeff_eff_h)
  deallocate(diag%temp_diff_coeff_eff_v)
  deallocate(diag%p_grad_decel_factor)
  deallocate(diag%condensates_sediment_heat)
  deallocate(diag%mass_diff_tendency)
  deallocate(diag%phase_trans_rates)
  deallocate(diag%phase_trans_heating_rate)
  deallocate(diag%viscosity)
  deallocate(diag%viscosity_rhombi)
  deallocate(diag%viscosity_triangles)
  deallocate(diag%vert_hor_viscosity)
  deallocate(diag%tke)
  deallocate(diag%sst)
  deallocate(diag%p_grad_acc_old_h)
  deallocate(diag%p_grad_acc_neg_nl_h)
  deallocate(diag%p_grad_acc_neg_nl_v)
  deallocate(diag%p_grad_acc_neg_l_h)
  deallocate(diag%p_grad_acc_neg_l_v)
  deallocate(diag%p_grad_condensates_v)
  deallocate(diag%v_squared_grad_h)
  deallocate(diag%v_squared_grad_v)
  deallocate(diag%pot_vort_tend_h)
  deallocate(diag%pot_vort_tend_v)
  deallocate(diag%sfc_sw_in)
  deallocate(diag%sfc_lw_out)
  deallocate(diag%radiation_tendency)
  deallocate(diag%a_rain)
  deallocate(state_write%rho)
  deallocate(state_write%rhotheta_v)
  deallocate(state_write%theta_v_pert)
  deallocate(state_write%exner_pert)
  deallocate(state_write%wind_h)
  deallocate(state_write%wind_v)
  deallocate(state_write%temperature_soil)
  write(*,*) stars
  call cpu_time(end_timestamp)
  write(*,fmt="(A,F9.3)") " Average speed:",(60._wp*run_span_min+300._wp)/((end_timestamp-init_timestamp)/omp_num_threads)
  write(*,*) "GAME over."

end program control










