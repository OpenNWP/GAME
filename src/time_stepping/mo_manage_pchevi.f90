! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_manage_pchevi

  ! This module manages the predictor-corrector HEVI time stepping.
  
  use mo_definitions,            only: wp,t_grid,t_diag
  use mo_constituents_nml,       only: n_constituents,lmoist,n_condensed_constituents
  use mo_grid_nml,               only: n_layers,n_vectors_per_layer,n_scalars_h,n_vectors,n_dual_vectors,n_scalars, &
                                       n_dual_scalars_h,n_dual_v_vectors,n_h_vectors,n_vectors_h
  use mo_run_nml,                only: dtime
  use mo_column_solvers,         only: three_band_solver_ver_waves,three_band_solver_gen_densities
  use mo_surface_nml,            only: nsoillays,lsfc_sensible_heat_flux,lsfc_phase_trans,pbl_scheme
  use mo_rad_nml,                only: rad_config
  use mo_pgrad,                  only: manage_pressure_gradient,calc_pressure_grad_condensates_v
  use mo_derived,                only: temperature_diagnostics
  use mo_pbl,                    only: update_sfc_turb_quantities
  use mo_scalar_tend_expl,       only: scalar_tend_expl
  use mo_vector_tend_expl,       only: vector_tend_expl
  use mo_phase_trans,            only: calc_h2otracers_source_rates
  use mo_manage_radiation_calls, only: call_radiation
  
  implicit none
  
  contains
  
  subroutine manage_pchevi(wind_old,wind_tend,wind_new, &
                           time_coordinate, &
                           totally_first_step_bool, &
                           theta_v_pert_old,theta_v_pert_new,temperature_soil_new, &
                           temperature_soil_old, &
                           rhotheta_v_tend,rhotheta_v_old, &
                           rhotheta_v_new,rho_tend,rho_old,rho_new, &
                           exner_pert_new,exner_pert_old,rad_update,diag,grid)
    
    real(wp), intent(out) :: wind_new(n_vectors),wind_tend(n_vectors), &
                             temperature_soil_new(nsoillays*n_scalars_h), &
                             theta_v_pert_new(n_scalars), &
                             rhotheta_v_tend(n_scalars),rhotheta_v_new(n_scalars), &
                             rho_tend(n_scalars*n_constituents),rho_new(n_scalars*n_constituents), &
                             exner_pert_new(n_scalars)
    integer,  intent(in)  :: totally_first_step_bool,rad_update
    real(wp), intent(in)  :: time_coordinate,theta_v_pert_old(n_scalars), &
                             wind_old(n_vectors),temperature_soil_old(nsoillays*n_scalars_h), &
                             rhotheta_v_old(n_scalars),rho_old(n_scalars*n_constituents), &
                             exner_pert_old(n_scalars)
    type(t_grid), intent(inout) :: grid
    type(t_diag), intent(inout) :: diag
    
    ! local variabels
    integer :: h_index,layer_index,vector_index,rk_step
    
    ! Preparations
    ! ------------
    
    ! diagnosing the temperature
    call  temperature_diagnostics(diag%temperature,grid%theta_v_bg,theta_v_pert_old,grid%exner_bg,exner_pert_old,rho_old)
    
    ! updating surface-related turbulence quantities if it is necessary
    if (lsfc_sensible_heat_flux .or. lsfc_phase_trans .or. pbl_scheme==1) then
      call update_sfc_turb_quantities(grid%is_land,grid%roughness_length,diag%monin_obukhov_length,grid%z_scalar,grid%z_vector, &
                                      grid%theta_v_bg,theta_v_pert_old,diag%v_squared,diag%roughness_velocity, &
                                      diag%scalar_flux_resistance)
    endif
    
    ! cloud microphysics
    if (lmoist) then
      call calc_h2otracers_source_rates(rho_old,diag%temperature,grid%layer_thickness,temperature_soil_old, &
                                        diag%phase_trans_rates,diag%phase_trans_heating_rate, &
                                        diag%scalar_flux_resistance,grid%is_land,diag%power_flux_density_latent)
    endif
    
    ! Radiation is updated here.
    if (rad_config>0 .and. rad_update==1) then
      call call_radiation(grid%latitude_scalar,grid%longitude_scalar,temperature_soil_old,grid%sfc_albedo,grid%z_scalar, &
                          grid%z_vector,rho_old,diag%temperature,diag%radiation_tendency, &
                          diag%sfc_sw_in,diag%sfc_lw_out,time_coordinate)
    endif
      
    ! Loop over the RK substeps
    ! -------------------------
    
    do rk_step=0,1
    
      ! state_old remains unchanged the whole time.
      ! At rk_step == 0, state_new contains garbage.
    
      ! 1.) explicit component of the momentum equation
      ! -----------------------------------------------
      ! Update of the pressure gradient.
      if (rk_step==0) then
        call manage_pressure_gradient(diag%pressure_gradient_acc_neg_nl,diag%pressure_gradient_acc_neg_l, &
                                      diag%pressure_gradient_decel_factor,diag%scalar_field_placeholder, &
                                      exner_pert_old,grid%theta_v_bg, &
                                      rho_old,diag%pgrad_acc_old, &
                                      theta_v_pert_old,grid%from_index,grid%to_index,grid%normal_distance,grid%exner_bg_grad, &
                                      grid%inner_product_weights,grid%slope,totally_first_step_bool)
      endif
      
      if (rk_step==0) then
       call calc_pressure_grad_condensates_v(diag%pressure_gradient_decel_factor, &
                                             rho_old,grid%gravity_m,diag%pressure_grad_condensates_v)
        ! Only the horizontal momentum is a forward tendency.
       call  vector_tend_expl(wind_old,diag%rel_vort_on_triangles,grid%z_vector,grid%z_vector_dual,diag%rel_vort, &
                              grid%vorticity_indices_triangles,grid%vorticity_signs_triangles,grid%normal_distance, &
                              grid%area_dual,grid%from_index,grid%to_index,grid%from_index_dual, &
                              grid%to_index_dual,grid%inner_product_weights, &
                              grid%slope,diag%temperature,diag%friction_acc,grid%adjacent_signs_h, &
                              grid%adjacent_vector_indices_h,grid%area, &
                              diag%molecular_diffusion_coeff,grid%normal_distance_dual, &
                              rho_old,diag%tke,diag%viscosity,diag%viscosity_triangles, &
                              grid%volume,diag%wind_div,diag%viscosity_rhombi,diag%vector_field_placeholder, &
                              diag%curl_of_vorticity, &
                              grid%gravity_m,grid%theta_v_bg,theta_v_pert_old, &
                              diag%scalar_field_placeholder,diag%n_squared,wind_tend, &
                              grid%density_to_rhombi_indices,grid%density_to_rhombi_weights,diag%dv_hdz,grid%exner_bg, &
                              exner_pert_old,grid%f_vec,diag%flux_density,diag%heating_diss,grid%layer_thickness, &
                              diag%monin_obukhov_length, &
                              diag%pot_vort,grid%roughness_length,grid%trsk_indices,grid%trsk_modified_curl_indices, &
                              grid%trsk_weights, &
                              diag%v_squared,diag%vert_hor_viscosity,grid%z_scalar,diag%pot_vort_tend,diag%v_squared_grad, &
                              diag%pressure_gradient_acc_neg_nl,diag%pressure_gradient_acc_neg_l,diag%pgrad_acc_old, &
                              diag%pressure_grad_condensates_v,rk_step,totally_first_step_bool,grid)
      endif
      if (rk_step==1) then
        call calc_pressure_grad_condensates_v(diag%pressure_gradient_decel_factor, &
                                              rho_new,grid%gravity_m,diag%pressure_grad_condensates_v)
        ! Only the horizontal momentum is a forward tendency.
        call vector_tend_expl(wind_new,diag%rel_vort_on_triangles,grid%z_vector,grid%z_vector_dual,diag%rel_vort, &
                              grid%vorticity_indices_triangles,grid%vorticity_signs_triangles,grid%normal_distance, &
                              grid%area_dual,grid%from_index,grid%to_index,grid%from_index_dual, &
                              grid%to_index_dual,grid%inner_product_weights, &
                              grid%slope,diag%temperature,diag%friction_acc,grid%adjacent_signs_h, &
                              grid%adjacent_vector_indices_h,grid%area, &
                              diag%molecular_diffusion_coeff,grid%normal_distance_dual, &
                              rho_new,diag%tke,diag%viscosity,diag%viscosity_triangles, &
                              grid%volume,diag%wind_div,diag%viscosity_rhombi, &
                              diag%vector_field_placeholder,diag%curl_of_vorticity, &
                              grid%gravity_m,grid%theta_v_bg,theta_v_pert_new, &
                              diag%scalar_field_placeholder,diag%n_squared,wind_tend, &
                              grid%density_to_rhombi_indices,grid%density_to_rhombi_weights,diag%dv_hdz,grid%exner_bg, &
                              exner_pert_new,grid%f_vec,diag%flux_density,diag%heating_diss,grid%layer_thickness, &
                              diag%monin_obukhov_length, &
                              diag%pot_vort,grid%roughness_length,grid%trsk_indices,grid%trsk_modified_curl_indices, &
                              grid%trsk_weights, &
                              diag%v_squared,diag%vert_hor_viscosity,grid%z_scalar,diag%pot_vort_tend,diag%v_squared_grad, &
                              diag%pressure_gradient_acc_neg_nl,diag%pressure_gradient_acc_neg_l,diag%pgrad_acc_old, &
                              diag%pressure_grad_condensates_v,rk_step,totally_first_step_bool,grid)
      endif
      
      ! time stepping for the horizontal momentum can be directly executed
      !$omp parallel do private(h_index,layer_index,vector_index)
      do h_index=1,n_vectors_h
        do layer_index=0,n_layers-1
          vector_index = n_scalars_h + layer_index*n_vectors_per_layer + h_index
          wind_new(vector_index) = wind_old(vector_index) + dtime*wind_tend(vector_index)
        enddo
      enddo
      !$omp end parallel do
      ! Horizontal velocity can be considered to be updated from now on.

      ! 2.) explicit component of the generalized density equations
      ! -----------------------------------------------------------
      if (rk_step==0) then
        call scalar_tend_expl(rho_old,diag%mass_diff_tendency,diag%scalar_field_placeholder,grid%adjacent_vector_indices_h, &
                              rhotheta_v_tend,grid%adjacent_signs_h,grid%area,diag%flux_density,grid%from_index,grid%to_index, &
                              grid%inner_product_weights,grid%layer_thickness,diag%mass_diffusion_coeff_numerical_h, &
                              diag%mass_diffusion_coeff_numerical_v,diag%molecular_diffusion_coeff,diag%n_squared, &
                              grid%normal_distance,rhotheta_v_old,grid%slope,diag%temp_diffusion_coeff_numerical_h, &
                              diag%temp_diffusion_coeff_numerical_v,diag%temperature, &
                              diag%tke,diag%vector_field_placeholder, &
                              diag%viscosity,diag%viscosity_rhombi,diag%viscosity_triangles,grid%volume, &
                              grid%vorticity_indices_triangles, &
                              wind_new,diag%temperature_diffusion_heating,diag%flux_density_div,rho_tend, &
                              diag%phase_trans_rates, &
                              exner_pert_old,grid%exner_bg,diag%condensates_sediment_heat,diag%phase_trans_heating_rate, &
                              diag%radiation_tendency,diag%heating_diss,rk_step)
      endif
      if (rk_step==1) then
        call scalar_tend_expl(rho_new,diag%mass_diff_tendency,diag%scalar_field_placeholder,grid%adjacent_vector_indices_h, &
                              rhotheta_v_tend,grid%adjacent_signs_h,grid%area,diag%flux_density,grid%from_index,grid%to_index, &
                              grid%inner_product_weights,grid%layer_thickness,diag%mass_diffusion_coeff_numerical_h, &
                              diag%mass_diffusion_coeff_numerical_v,diag%molecular_diffusion_coeff,diag%n_squared, &
                              grid%normal_distance,rhotheta_v_new,grid%slope,diag%temp_diffusion_coeff_numerical_h, &
                              diag%temp_diffusion_coeff_numerical_v,diag%temperature, &
                              diag%tke,diag%vector_field_placeholder, &
                              diag%viscosity,diag%viscosity_rhombi,diag%viscosity_triangles, &
                              grid%volume,grid%vorticity_indices_triangles, &
                              wind_new,diag%temperature_diffusion_heating,diag%flux_density_div,rho_tend, &
                              diag%phase_trans_rates, &
                              exner_pert_new,grid%exner_bg,diag%condensates_sediment_heat,diag%phase_trans_heating_rate, &
                              diag%radiation_tendency,diag%heating_diss,rk_step)
      endif

      ! 3.) vertical sound wave solver
      ! ------------------------------
      if (rk_step==0) then
        call three_band_solver_ver_waves(grid%z_vector,grid%area,grid%volume,diag%scalar_flux_resistance,theta_v_pert_new, &
                                         rho_new,rhotheta_v_old,rho_old,rho_old,theta_v_pert_old,grid%theta_v_bg, &
                                         exner_pert_old,grid%exner_bg,theta_v_pert_old,exner_pert_old, &
                                         diag%power_flux_density_sensible,temperature_soil_old,temperature_soil_old, &
                                         rhotheta_v_tend,grid%is_land,rho_tend,rhotheta_v_old,wind_old, &
                                         grid%z_scalar,wind_tend,grid%z_soil_center,grid%t_conduc_soil,grid%sfc_rho_c, &
                                         diag%sfc_lw_out,diag%sfc_sw_in,grid%t_const_soil,grid%z_soil_interface, &
                                         diag%power_flux_density_latent,rhotheta_v_new,exner_pert_new, &
                                         wind_new,temperature_soil_new,rk_step)
      endif
      if (rk_step==1) then
        call three_band_solver_ver_waves(grid%z_vector,grid%area,grid%volume,diag%scalar_flux_resistance,theta_v_pert_new, &
                                         rho_new,rhotheta_v_new,rho_new,rho_old,theta_v_pert_old,grid%theta_v_bg, &
                                         exner_pert_old,grid%exner_bg,theta_v_pert_new,exner_pert_new, &
                                         diag%power_flux_density_sensible,temperature_soil_new,temperature_soil_old, &
                                         rhotheta_v_tend,grid%is_land,rho_tend,rhotheta_v_old,wind_old, &
                                         grid%z_scalar,wind_tend,grid%z_soil_center,grid%t_conduc_soil,grid%sfc_rho_c, &
                                         diag%sfc_lw_out,diag%sfc_sw_in,grid%t_const_soil,grid%z_soil_interface, &
                                         diag%power_flux_density_latent,rhotheta_v_new,exner_pert_new, &
                                         wind_new,temperature_soil_new,rk_step)
      endif
      
      ! 4.) vertical tracer advection
      ! -----------------------------
      if (n_constituents>1) then
        call three_band_solver_gen_densities(wind_old,wind_new,grid%volume,rho_tend,rho_old,rho_new, &
                                             diag%condensates_sediment_heat,grid%area,diag%temperature,rk_step)
      endif
      
    enddo
    
  end subroutine manage_pchevi

end module mo_manage_pchevi





