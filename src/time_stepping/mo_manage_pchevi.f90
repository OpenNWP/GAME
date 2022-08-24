! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_manage_pchevi

  ! This module manages the predictor-corrector HEVI time stepping.
  
  use iso_c_binding
  use mo_definitions,           only: wp
  use constituents_nml,         only: n_constituents,lmoist,n_condensed_constituents
  use grid_nml,                 only: n_layers,n_vectors_per_layer,n_scalars_h,n_vectors,n_dual_vectors,n_scalars, &
                                      n_dual_scalars_h,n_dual_v_vectors,n_h_vectors,n_vectors_h
  use run_nml,                  only: dtime
  use column_solvers,           only: three_band_solver_ver_waves,three_band_solver_gen_densities
  use surface_nml,              only: nsoillays,lsfc_sensible_heat_flux,lsfc_phase_trans,pbl_scheme
  use rad_nml,                  only: rad_config
  use mo_pgrad,                 only: manage_pressure_gradient,calc_pressure_grad_condensates_v
  use derived_quantities,       only: temperature_diagnostics
  use planetary_boundary_layer, only: update_sfc_turb_quantities
  use mo_explicit_scalar_tend,  only: scalar_tendencies_expl
  use mo_explicit_wind_tend,    only: vector_tendencies_expl
  use phase_trans,              only: calc_h2otracers_source_rates
  use manage_radiation_calls,   only: call_radiation
  
  implicit none
  
  contains
  
  subroutine manage_pchevi(adjacent_signs_h,adjacent_vector_indices_h,area,layer_thickness, &
                           z_scalar,z_vector,volume,vorticity_indices_triangles,vorticity_signs_triangles, &
                           z_t_const,z_soil_center,z_soil_interface,v_squared,trsk_weights, &
                           from_index,to_index,from_index_dual,to_index_dual,trsk_modified_curl_indices, &
                           area_dual,z_vector_dual,wind_old,wind_tend,wind_new,trsk_indices, &
                           temperature,wind_div,viscosity_triangles,viscosity,viscosity_rhombi, &
                           condensates_sediment_heat,molecular_diffusion_coeff,v_squared_grad, &
                           time_coordinate,vert_hor_viscosity,vector_field_placeholder,sfc_sw_in, &
                           totally_first_step_bool,gravity_m,curl_of_vorticity,tke,t_const_soil, &
                           theta_v_pert_old,theta_v_pert_new,theta_v_bg,temperature_soil_new,sfc_lw_out, &
                           temperature_soil_old,temperature_diffusion_heating,slope,t_conduc_soil, &
                           temp_diffusion_coeff_numerical_h,temp_diffusion_coeff_numerical_v,friction_acc, &
                           sfc_rho_c,sfc_albedo,scalar_flux_resistance,scalar_field_placeholder, &
                           roughness_velocity,roughness_length,rhotheta_v_tend,rhotheta_v_old, &
                           rhotheta_v_new,rho_tend,rho_old,rho_new,radiation_tendency,exner_bg, &
                           pot_vort_tend,normal_distance,n_squared,inner_product_weights,exner_bg_grad, &
                           normal_distance_dual,power_flux_density_latent,power_flux_density_sensible, &
                           density_to_rhombi_weights,density_to_rhombi_indices,rel_vort_on_triangles, &
                           phase_trans_heating_rate,exner_pert_new,exner_pert_old,rel_vort,f_vec,rad_update, &
                           pressure_gradient_decel_factor,pressure_gradient_acc_neg_l,pressure_gradient_acc_neg_nl, &
                           pot_vort,flux_density,pressure_grad_condensates_v,flux_density_div,dv_hdz, &
                           monin_obukhov_length,heating_diss,is_land,pgrad_acc_old,mass_diff_tendency, &
                           latitude_scalar,longitude_scalar,mass_diffusion_coeff_numerical_h, &
                           mass_diffusion_coeff_numerical_v,phase_trans_rates) &
  bind(c,name = "manage_pchevi")
    
    real(wp), intent(out) :: wind_new(n_vectors),wind_tend(n_vectors),condensates_sediment_heat(n_scalars), &
                             temperature(n_scalars),wind_div(n_scalars),temperature_soil_new(nsoillays*n_scalars_h), &
                             viscosity_rhombi(n_vectors),viscosity(n_scalars),curl_of_vorticity(n_vectors), &
                             molecular_diffusion_coeff(n_scalars),vert_hor_viscosity(n_h_vectors+n_vectors_h), &
                             vector_field_placeholder(n_vectors),v_squared_grad(n_vectors),v_squared(n_scalars), &
                             tke(n_scalars),theta_v_pert_new(n_scalars),temperature_diffusion_heating(n_scalars), &
                             temp_diffusion_coeff_numerical_h(n_scalars),temp_diffusion_coeff_numerical_v(n_scalars), &
                             sfc_sw_in(n_scalars_h),sfc_lw_out(n_scalars_h),scalar_flux_resistance(n_scalars_h), &
                             scalar_field_placeholder(n_scalars),roughness_velocity(n_scalars_h), &
                             roughness_length(n_scalars_h),rhotheta_v_tend(n_scalars),rhotheta_v_new(n_scalars), &
                             rho_tend(n_scalars*n_constituents),rho_new(n_scalars*n_constituents), &
                             radiation_tendency(n_scalars),pot_vort_tend(n_vectors),n_squared(n_scalars), &
                             power_flux_density_latent(n_scalars_h),power_flux_density_sensible(n_scalars_h), &
                             rel_vort_on_triangles(n_dual_v_vectors),friction_acc(n_vectors), &
                             phase_trans_heating_rate(n_scalars),exner_pert_new(n_scalars), &
                             rel_vort((2*n_layers+1)*n_vectors_h),pressure_gradient_decel_factor(n_scalars), &
                             pressure_gradient_acc_neg_l(n_vectors),pressure_gradient_acc_neg_nl(n_vectors), &
                             pot_vort((2*n_layers+1)*n_vectors_h),flux_density(n_vectors), &
                             pressure_grad_condensates_v(n_vectors),flux_density_div(n_scalars), &
                             monin_obukhov_length(n_scalars_h),heating_diss(n_scalars),pgrad_acc_old(n_vectors), &
                             mass_diff_tendency(n_scalars),mass_diffusion_coeff_numerical_h(n_scalars), &
                             mass_diffusion_coeff_numerical_v(n_scalars),dv_hdz(n_h_vectors+n_vectors_h), &
                             phase_trans_rates((n_condensed_constituents+1)*n_scalars),viscosity_triangles(n_dual_v_vectors)
    integer,  intent(in)  :: adjacent_signs_h(6*n_scalars_h),adjacent_vector_indices_h(6*n_scalars_h), &
                             totally_first_step_bool,vorticity_indices_triangles(3*n_dual_scalars_h), &
                             vorticity_signs_triangles(3*n_dual_scalars_h),trsk_indices(10*n_vectors_h), &
                             from_index(n_vectors_h),to_index(n_vectors_h),from_index_dual(n_vectors_h), &
                             to_index_dual(n_vectors_h),trsk_modified_curl_indices(10*n_vectors_h), &
                             density_to_rhombi_indices(4*n_vectors_h),rad_update,is_land(n_scalars_h)
    real(wp), intent(in)  :: area(n_vectors),latitude_scalar(n_scalars_h),longitude_scalar(n_scalars_h), &
                             area_dual(n_dual_vectors),layer_thickness(n_scalars),z_vector(n_vectors), &
                             z_vector_dual(n_dual_vectors),z_t_const,z_soil_center(nsoillays), &
                             z_soil_interface(nsoillays+1),z_scalar(n_scalars),volume(n_scalars), &
                             time_coordinate,gravity_m(n_vectors), &
                             trsk_weights(10*n_vectors_h),theta_v_bg(n_scalars),theta_v_pert_old(n_scalars), &
                             wind_old(n_vectors),temperature_soil_old(nsoillays*n_scalars_h),slope(n_vectors), &
                             t_const_soil(n_scalars_h),t_conduc_soil(n_scalars_h),sfc_rho_c(n_scalars_h), &
                             sfc_albedo(n_scalars_h),rhotheta_v_old(n_scalars),rho_old(n_scalars*n_constituents), &
                             exner_bg(n_scalars),normal_distance(n_vectors),inner_product_weights(8*n_scalars), &
                             exner_bg_grad(n_vectors),normal_distance_dual(n_dual_vectors),f_vec(2*n_vectors_h), &
                             density_to_rhombi_weights(4*n_vectors_h),exner_pert_old(n_scalars)
    
    ! local variabels
    integer :: h_index,layer_index,vector_index,rk_step
    
    ! Preparations
    ! ------------
    
    ! diagnosing the temperature
    call  temperature_diagnostics(temperature,theta_v_bg,theta_v_pert_old,exner_bg,exner_pert_old,rho_old)
    
    ! updating surface-related turbulence quantities if it is necessary
    if (lsfc_sensible_heat_flux .or. lsfc_phase_trans .or. pbl_scheme==1) then
      call update_sfc_turb_quantities(is_land,roughness_length,monin_obukhov_length,z_scalar,z_vector, &
                                      theta_v_bg,theta_v_pert_old,v_squared,roughness_velocity,scalar_flux_resistance)
    endif
    
    ! cloud microphysics
    if (lmoist) then
      call calc_h2otracers_source_rates(rho_old,temperature,layer_thickness,temperature_soil_old, &
                                        phase_trans_rates,phase_trans_heating_rate, &
                                        scalar_flux_resistance,is_land,power_flux_density_latent)
    endif
    
    ! Radiation is updated here.
    if (rad_config>0 .and. rad_update==1) then
      call call_radiation(latitude_scalar,longitude_scalar,temperature_soil_old,sfc_albedo,z_scalar, &
                          z_vector,rho_old,temperature,radiation_tendency, &
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
                                      pressure_gradient_decel_factor,scalar_field_placeholder,exner_pert_old,theta_v_bg, &
                                      rho_old, pgrad_acc_old, &
                                      theta_v_pert_old,from_index,to_index,normal_distance, exner_bg_grad, &
                                      inner_product_weights,slope,totally_first_step_bool)
      endif
      
      if (rk_step==0) then
       call calc_pressure_grad_condensates_v(pressure_gradient_decel_factor,rho_old,gravity_m,pressure_grad_condensates_v)
        ! Only the horizontal momentum is a forward tendency.
       call  vector_tendencies_expl(wind_old,rel_vort_on_triangles,z_vector,z_vector_dual,rel_vort, &
                                    vorticity_indices_triangles,vorticity_signs_triangles,normal_distance, &
                                    area_dual,from_index,to_index,from_index_dual,to_index_dual,inner_product_weights, &
                                    slope,temperature,friction_acc,adjacent_signs_h,adjacent_vector_indices_h,area, &
                                    molecular_diffusion_coeff,normal_distance_dual,rho_old,tke,viscosity,viscosity_triangles, &
                                    volume,wind_div,viscosity_rhombi,vector_field_placeholder,curl_of_vorticity, &
                                    gravity_m,theta_v_bg,theta_v_pert_old,scalar_field_placeholder,n_squared,wind_tend, &
                                    density_to_rhombi_indices,density_to_rhombi_weights,dv_hdz,exner_bg, &
                                    exner_pert_old,f_vec,flux_density,heating_diss,layer_thickness,monin_obukhov_length, &
                                    pot_vort,roughness_length,trsk_indices,trsk_modified_curl_indices,trsk_weights, &
                                    v_squared,vert_hor_viscosity,z_scalar,pot_vort_tend,v_squared_grad, &
                                    pressure_gradient_acc_neg_nl,pressure_gradient_acc_neg_l,pgrad_acc_old, &
                                    pressure_grad_condensates_v,rk_step,totally_first_step_bool)
      endif
      if (rk_step==1) then
        call calc_pressure_grad_condensates_v(pressure_gradient_decel_factor,rho_new,gravity_m,pressure_grad_condensates_v)
        ! Only the horizontal momentum is a forward tendency.
        call vector_tendencies_expl(wind_new,rel_vort_on_triangles,z_vector,z_vector_dual,rel_vort, &
                                    vorticity_indices_triangles,vorticity_signs_triangles,normal_distance, &
                                    area_dual,from_index,to_index,from_index_dual,to_index_dual,inner_product_weights, &
                                    slope,temperature,friction_acc,adjacent_signs_h,adjacent_vector_indices_h,area, &
                                    molecular_diffusion_coeff,normal_distance_dual,rho_new,tke,viscosity,viscosity_triangles, &
                                    volume,wind_div,viscosity_rhombi,vector_field_placeholder,curl_of_vorticity, &
                                    gravity_m,theta_v_bg,theta_v_pert_new,scalar_field_placeholder,n_squared,wind_tend, &
                                    density_to_rhombi_indices,density_to_rhombi_weights,dv_hdz,exner_bg, &
                                    exner_pert_new,f_vec,flux_density,heating_diss,layer_thickness,monin_obukhov_length, &
                                    pot_vort,roughness_length,trsk_indices,trsk_modified_curl_indices,trsk_weights, &
                                    v_squared,vert_hor_viscosity,z_scalar,pot_vort_tend,v_squared_grad, &
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
        call scalar_tendencies_expl(rho_old,mass_diff_tendency,scalar_field_placeholder,adjacent_vector_indices_h, &
                                    rhotheta_v_tend,adjacent_signs_h,area,flux_density,from_index,to_index, &
                                    inner_product_weights,layer_thickness,mass_diffusion_coeff_numerical_h, &
                                    mass_diffusion_coeff_numerical_v,molecular_diffusion_coeff,n_squared, &
                                    normal_distance,rhotheta_v_old,slope,temp_diffusion_coeff_numerical_h, &
                                    temp_diffusion_coeff_numerical_v,temperature,tke,vector_field_placeholder, &
                                    viscosity,viscosity_rhombi,viscosity_triangles,volume,vorticity_indices_triangles, &
                                    wind_new,temperature_diffusion_heating,flux_density_div,rho_tend,phase_trans_rates, &
                                    exner_pert_old,exner_bg,condensates_sediment_heat,phase_trans_heating_rate, &
                                    radiation_tendency,heating_diss,rk_step)
      endif
      if (rk_step==1) then
        call scalar_tendencies_expl(rho_new,mass_diff_tendency,scalar_field_placeholder,adjacent_vector_indices_h, &
                                    rhotheta_v_tend,adjacent_signs_h,area,flux_density,from_index,to_index, &
                                    inner_product_weights,layer_thickness,mass_diffusion_coeff_numerical_h, &
                                    mass_diffusion_coeff_numerical_v,molecular_diffusion_coeff,n_squared, &
                                    normal_distance,rhotheta_v_new,slope,temp_diffusion_coeff_numerical_h, &
                                    temp_diffusion_coeff_numerical_v,temperature,tke,vector_field_placeholder, &
                                    viscosity,viscosity_rhombi,viscosity_triangles,volume,vorticity_indices_triangles, &
                                    wind_new,temperature_diffusion_heating,flux_density_div,rho_tend,phase_trans_rates, &
                                    exner_pert_new,exner_bg,condensates_sediment_heat,phase_trans_heating_rate, &
                                    radiation_tendency,heating_diss,rk_step)
      endif

      ! 3.) vertical sound wave solver
      ! ------------------------------
      if (rk_step==0) then
        call three_band_solver_ver_waves(z_vector,area,volume,scalar_flux_resistance,theta_v_pert_new, &
                                         rho_new,rhotheta_v_old,rho_old,rho_old,theta_v_pert_old,theta_v_bg, &
                                         exner_pert_old,exner_bg,theta_v_pert_old,exner_pert_old, &
                                         power_flux_density_sensible,temperature_soil_old,temperature_soil_old, &
                                         rhotheta_v_tend,is_land,rho_tend,rhotheta_v_old,wind_old, &
                                         z_scalar,wind_tend,z_soil_center,t_conduc_soil,sfc_rho_c, &
                                         sfc_lw_out,sfc_sw_in,t_const_soil,z_soil_interface, &
                                         power_flux_density_latent,rhotheta_v_new,exner_pert_new, &
                                         wind_new,temperature_soil_new,z_t_const,rk_step)
      endif
      if (rk_step==1) then
        call three_band_solver_ver_waves(z_vector,area,volume,scalar_flux_resistance,theta_v_pert_new, &
                                         rho_new,rhotheta_v_new,rho_new,rho_old,theta_v_pert_old,theta_v_bg, &
                                         exner_pert_old,exner_bg,theta_v_pert_new,exner_pert_new, &
                                         power_flux_density_sensible,temperature_soil_new,temperature_soil_old, &
                                         rhotheta_v_tend,is_land,rho_tend,rhotheta_v_old,wind_old, &
                                         z_scalar,wind_tend,z_soil_center,t_conduc_soil,sfc_rho_c, &
                                         sfc_lw_out,sfc_sw_in,t_const_soil,z_soil_interface, &
                                         power_flux_density_latent,rhotheta_v_new,exner_pert_new, &
                                         wind_new,temperature_soil_new,z_t_const,rk_step)
      endif
      
      ! 4.) vertical tracer advection
      ! -----------------------------
      if (n_constituents>1) then
        call three_band_solver_gen_densities(wind_old,wind_new,volume,rho_tend,rho_old,rho_new, &
                                             condensates_sediment_heat,area,temperature,rk_step)
      endif
      
    enddo
    
  end subroutine manage_pchevi

end module mo_manage_pchevi





