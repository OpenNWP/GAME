! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_manage_pchevi

  ! This module manages the predictor-corrector HEVI time stepping.
  
  use mo_definitions,            only: wp,t_grid
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
  
  subroutine manage_pchevi(v_squared,wind_old,wind_tend,wind_new, &
                           temperature,wind_div,viscosity_triangles,viscosity,viscosity_rhombi, &
                           condensates_sediment_heat,molecular_diffusion_coeff,v_squared_grad, &
                           time_coordinate,vert_hor_viscosity,vector_field_placeholder,sfc_sw_in, &
                           totally_first_step_bool,curl_of_vorticity,tke, &
                           theta_v_pert_old,theta_v_pert_new,temperature_soil_new,sfc_lw_out, &
                           temperature_soil_old,temperature_diffusion_heating, &
                           temp_diffusion_coeff_numerical_h,temp_diffusion_coeff_numerical_v,friction_acc, &
                           scalar_flux_resistance,scalar_field_placeholder, &
                           roughness_velocity,rhotheta_v_tend,rhotheta_v_old, &
                           rhotheta_v_new,rho_tend,rho_old,rho_new,radiation_tendency, &
                           pot_vort_tend,n_squared, &
                           power_flux_density_latent,power_flux_density_sensible, &
                           rel_vort_on_triangles, &
                           phase_trans_heating_rate,exner_pert_new,exner_pert_old,rel_vort,rad_update, &
                           pressure_gradient_decel_factor,pressure_gradient_acc_neg_l,pressure_gradient_acc_neg_nl, &
                           pot_vort,flux_density,pressure_grad_condensates_v,flux_density_div,dv_hdz, &
                           monin_obukhov_length,heating_diss,pgrad_acc_old,mass_diff_tendency, &
                           mass_diffusion_coeff_numerical_h, &
                           mass_diffusion_coeff_numerical_v,phase_trans_rates,grid)
    
    real(wp), intent(out) :: wind_new(n_vectors),wind_tend(n_vectors),condensates_sediment_heat(n_scalars), &
                             temperature(n_scalars),wind_div(n_scalars),temperature_soil_new(nsoillays*n_scalars_h), &
                             viscosity_rhombi(n_vectors),viscosity(n_scalars),curl_of_vorticity(n_vectors), &
                             molecular_diffusion_coeff(n_scalars),vert_hor_viscosity(n_h_vectors+n_vectors_h), &
                             vector_field_placeholder(n_vectors),v_squared_grad(n_vectors),v_squared(n_scalars), &
                             tke(n_scalars),theta_v_pert_new(n_scalars),temperature_diffusion_heating(n_scalars), &
                             temp_diffusion_coeff_numerical_h(n_scalars),temp_diffusion_coeff_numerical_v(n_scalars), &
                             sfc_sw_in(n_scalars_h),sfc_lw_out(n_scalars_h),scalar_flux_resistance(n_scalars_h), &
                             scalar_field_placeholder(n_scalars),roughness_velocity(n_scalars_h), &
                             rhotheta_v_tend(n_scalars),rhotheta_v_new(n_scalars), &
                             rho_tend(n_scalars*n_constituents),rho_new(n_scalars*n_constituents), &
                             radiation_tendency(n_scalars),pot_vort_tend(n_vectors),n_squared(n_scalars), &
                             power_flux_density_latent(n_scalars_h),power_flux_density_sensible(n_scalars_h), &
                             rel_vort_on_triangles(n_dual_v_vectors),friction_acc(n_vectors), &
                             phase_trans_heating_rate(n_scalars),exner_pert_new(n_scalars), &
                             rel_vort((2*n_layers+1)*n_vectors_h),pressure_gradient_decel_factor(n_scalars), &
                             pressure_gradient_acc_neg_l(n_vectors),pressure_gradient_acc_neg_nl(n_vectors), &
                             pot_vort((2*n_layers+1)*n_vectors_h),flux_density(n_vectors), &
                             pressure_grad_condensates_v(n_vectors),flux_density_div(n_scalars), &
                             heating_diss(n_scalars),pgrad_acc_old(n_vectors), &
                             mass_diff_tendency(n_scalars),mass_diffusion_coeff_numerical_h(n_scalars), &
                             mass_diffusion_coeff_numerical_v(n_scalars),dv_hdz(n_h_vectors+n_vectors_h), &
                             phase_trans_rates((n_condensed_constituents+1)*n_scalars),viscosity_triangles(n_dual_v_vectors), &
                             monin_obukhov_length(n_scalars_h)
    integer,  intent(in)  :: totally_first_step_bool,rad_update
    real(wp), intent(in)  :: time_coordinate,theta_v_pert_old(n_scalars), &
                             wind_old(n_vectors),temperature_soil_old(nsoillays*n_scalars_h), &
                             rhotheta_v_old(n_scalars),rho_old(n_scalars*n_constituents), &
                             exner_pert_old(n_scalars)
    type(t_grid), intent(inout) :: grid
    
    ! local variabels
    integer :: h_index,layer_index,vector_index,rk_step
    
    ! Preparations
    ! ------------
    
    ! diagnosing the temperature
    call  temperature_diagnostics(temperature,grid%theta_v_bg,theta_v_pert_old,grid%exner_bg,exner_pert_old,rho_old)
    
    ! updating surface-related turbulence quantities if it is necessary
    if (lsfc_sensible_heat_flux .or. lsfc_phase_trans .or. pbl_scheme==1) then
      call update_sfc_turb_quantities(grid%is_land,grid%roughness_length,monin_obukhov_length,grid%z_scalar,grid%z_vector, &
                                      grid%theta_v_bg,theta_v_pert_old,v_squared,roughness_velocity,scalar_flux_resistance)
    endif
    
    ! cloud microphysics
    if (lmoist) then
      call calc_h2otracers_source_rates(rho_old,temperature,grid%layer_thickness,temperature_soil_old, &
                                        phase_trans_rates,phase_trans_heating_rate, &
                                        scalar_flux_resistance,grid%is_land,power_flux_density_latent)
    endif
    
    ! Radiation is updated here.
    if (rad_config>0 .and. rad_update==1) then
      call call_radiation(grid%latitude_scalar,grid%longitude_scalar,temperature_soil_old,grid%sfc_albedo,grid%z_scalar, &
                          grid%z_vector,rho_old,temperature,radiation_tendency, &
                          sfc_sw_in,sfc_lw_out,time_coordinate)
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
        call manage_pressure_gradient(pressure_gradient_acc_neg_nl,pressure_gradient_acc_neg_l, &
                                      pressure_gradient_decel_factor,scalar_field_placeholder,exner_pert_old,grid%theta_v_bg, &
                                      rho_old, pgrad_acc_old, &
                                      theta_v_pert_old,grid%from_index,grid%to_index,grid%normal_distance,grid%exner_bg_grad, &
                                      grid%inner_product_weights,grid%slope,totally_first_step_bool)
      endif
      
      if (rk_step==0) then
       call calc_pressure_grad_condensates_v(pressure_gradient_decel_factor,rho_old,grid%gravity_m,pressure_grad_condensates_v)
        ! Only the horizontal momentum is a forward tendency.
       call  vector_tend_expl(wind_old,rel_vort_on_triangles,grid%z_vector,grid%z_vector_dual,rel_vort, &
                              grid%vorticity_indices_triangles,grid%vorticity_signs_triangles,grid%normal_distance, &
                              grid%area_dual,grid%from_index,grid%to_index,grid%from_index_dual, &
                              grid%to_index_dual,grid%inner_product_weights, &
                              grid%slope,temperature,friction_acc,grid%adjacent_signs_h,grid%adjacent_vector_indices_h,grid%area, &
                              molecular_diffusion_coeff,grid%normal_distance_dual, &
                              rho_old,tke,viscosity,viscosity_triangles, &
                              grid%volume,wind_div,viscosity_rhombi,vector_field_placeholder,curl_of_vorticity, &
                              grid%gravity_m,grid%theta_v_bg,theta_v_pert_old,scalar_field_placeholder,n_squared,wind_tend, &
                              grid%density_to_rhombi_indices,grid%density_to_rhombi_weights,dv_hdz,grid%exner_bg, &
                              exner_pert_old,grid%f_vec,flux_density,heating_diss,grid%layer_thickness, &
                              monin_obukhov_length, &
                              pot_vort,grid%roughness_length,grid%trsk_indices,grid%trsk_modified_curl_indices, &
                              grid%trsk_weights, &
                              v_squared,vert_hor_viscosity,grid%z_scalar,pot_vort_tend,v_squared_grad, &
                              pressure_gradient_acc_neg_nl,pressure_gradient_acc_neg_l,pgrad_acc_old, &
                              pressure_grad_condensates_v,rk_step,totally_first_step_bool)
      endif
      if (rk_step==1) then
        call calc_pressure_grad_condensates_v(pressure_gradient_decel_factor,rho_new,grid%gravity_m,pressure_grad_condensates_v)
        ! Only the horizontal momentum is a forward tendency.
        call vector_tend_expl(wind_new,rel_vort_on_triangles,grid%z_vector,grid%z_vector_dual,rel_vort, &
                              grid%vorticity_indices_triangles,grid%vorticity_signs_triangles,grid%normal_distance, &
                              grid%area_dual,grid%from_index,grid%to_index,grid%from_index_dual, &
                              grid%to_index_dual,grid%inner_product_weights, &
                              grid%slope,temperature,friction_acc,grid%adjacent_signs_h,grid%adjacent_vector_indices_h,grid%area, &
                              molecular_diffusion_coeff,grid%normal_distance_dual, &
                              rho_new,tke,viscosity,viscosity_triangles, &
                              grid%volume,wind_div,viscosity_rhombi,vector_field_placeholder,curl_of_vorticity, &
                              grid%gravity_m,grid%theta_v_bg,theta_v_pert_new,scalar_field_placeholder,n_squared,wind_tend, &
                              grid%density_to_rhombi_indices,grid%density_to_rhombi_weights,dv_hdz,grid%exner_bg, &
                              exner_pert_new,grid%f_vec,flux_density,heating_diss,grid%layer_thickness, &
                              monin_obukhov_length, &
                              pot_vort,grid%roughness_length,grid%trsk_indices,grid%trsk_modified_curl_indices, &
                              grid%trsk_weights, &
                              v_squared,vert_hor_viscosity,grid%z_scalar,pot_vort_tend,v_squared_grad, &
                              pressure_gradient_acc_neg_nl,pressure_gradient_acc_neg_l,pgrad_acc_old, &
                              pressure_grad_condensates_v,rk_step,totally_first_step_bool)
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
        call scalar_tend_expl(rho_old,mass_diff_tendency,scalar_field_placeholder,grid%adjacent_vector_indices_h, &
                              rhotheta_v_tend,grid%adjacent_signs_h,grid%area,flux_density,grid%from_index,grid%to_index, &
                              grid%inner_product_weights,grid%layer_thickness,mass_diffusion_coeff_numerical_h, &
                              mass_diffusion_coeff_numerical_v,molecular_diffusion_coeff,n_squared, &
                              grid%normal_distance,rhotheta_v_old,grid%slope,temp_diffusion_coeff_numerical_h, &
                              temp_diffusion_coeff_numerical_v,temperature,tke,vector_field_placeholder, &
                              viscosity,viscosity_rhombi,viscosity_triangles,grid%volume,grid%vorticity_indices_triangles, &
                              wind_new,temperature_diffusion_heating,flux_density_div,rho_tend,phase_trans_rates, &
                              exner_pert_old,grid%exner_bg,condensates_sediment_heat,phase_trans_heating_rate, &
                              radiation_tendency,heating_diss,rk_step)
      endif
      if (rk_step==1) then
        call scalar_tend_expl(rho_new,mass_diff_tendency,scalar_field_placeholder,grid%adjacent_vector_indices_h, &
                              rhotheta_v_tend,grid%adjacent_signs_h,grid%area,flux_density,grid%from_index,grid%to_index, &
                              grid%inner_product_weights,grid%layer_thickness,mass_diffusion_coeff_numerical_h, &
                              mass_diffusion_coeff_numerical_v,molecular_diffusion_coeff,n_squared, &
                              grid%normal_distance,rhotheta_v_new,grid%slope,temp_diffusion_coeff_numerical_h, &
                              temp_diffusion_coeff_numerical_v,temperature,tke,vector_field_placeholder, &
                              viscosity,viscosity_rhombi,viscosity_triangles,grid%volume,grid%vorticity_indices_triangles, &
                              wind_new,temperature_diffusion_heating,flux_density_div,rho_tend,phase_trans_rates, &
                              exner_pert_new,grid%exner_bg,condensates_sediment_heat,phase_trans_heating_rate, &
                              radiation_tendency,heating_diss,rk_step)
      endif

      ! 3.) vertical sound wave solver
      ! ------------------------------
      if (rk_step==0) then
        call three_band_solver_ver_waves(grid%z_vector,grid%area,grid%volume,scalar_flux_resistance,theta_v_pert_new, &
                                         rho_new,rhotheta_v_old,rho_old,rho_old,theta_v_pert_old,grid%theta_v_bg, &
                                         exner_pert_old,grid%exner_bg,theta_v_pert_old,exner_pert_old, &
                                         power_flux_density_sensible,temperature_soil_old,temperature_soil_old, &
                                         rhotheta_v_tend,grid%is_land,rho_tend,rhotheta_v_old,wind_old, &
                                         grid%z_scalar,wind_tend,grid%z_soil_center,grid%t_conduc_soil,grid%sfc_rho_c, &
                                         sfc_lw_out,sfc_sw_in,grid%t_const_soil,grid%z_soil_interface, &
                                         power_flux_density_latent,rhotheta_v_new,exner_pert_new, &
                                         wind_new,temperature_soil_new,rk_step)
      endif
      if (rk_step==1) then
        call three_band_solver_ver_waves(grid%z_vector,grid%area,grid%volume,scalar_flux_resistance,theta_v_pert_new, &
                                         rho_new,rhotheta_v_new,rho_new,rho_old,theta_v_pert_old,grid%theta_v_bg, &
                                         exner_pert_old,grid%exner_bg,theta_v_pert_new,exner_pert_new, &
                                         power_flux_density_sensible,temperature_soil_new,temperature_soil_old, &
                                         rhotheta_v_tend,grid%is_land,rho_tend,rhotheta_v_old,wind_old, &
                                         grid%z_scalar,wind_tend,grid%z_soil_center,grid%t_conduc_soil,grid%sfc_rho_c, &
                                         sfc_lw_out,sfc_sw_in,grid%t_const_soil,grid%z_soil_interface, &
                                         power_flux_density_latent,rhotheta_v_new,exner_pert_new, &
                                         wind_new,temperature_soil_new,rk_step)
      endif
      
      ! 4.) vertical tracer advection
      ! -----------------------------
      if (n_constituents>1) then
        call three_band_solver_gen_densities(wind_old,wind_new,grid%volume,rho_tend,rho_old,rho_new, &
                                             condensates_sediment_heat,grid%area,temperature,rk_step)
      endif
      
    enddo
    
  end subroutine manage_pchevi

end module mo_manage_pchevi





