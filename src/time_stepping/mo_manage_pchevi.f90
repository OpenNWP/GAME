! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_manage_pchevi

  ! This module manages the predictor-corrector HEVI time stepping.
  
  use iso_c_binding
  use mo_definitions,   only: wp
  use constituents_nml, only: n_constituents
  use grid_nml,         only: n_layers,n_vectors_per_layer,n_scalars_h,n_vectors,n_dual_vectors,n_scalars, &
                              n_dual_scalars_h,n_dual_v_vectors
  use run_nml,          only: dtime
  use column_solvers,   only: three_band_solver_ver_waves,three_band_solver_gen_densities
  use surface_nml,      only: nsoillays
  
  implicit none
  
  contains
  
  subroutine manage_pchevi(adjacent_signs_h,adjacent_vector_indices_h,area,layer_thickness, &
                           z_scalar,z_vector,volume,vorticity_indices_triangles,vorticity_signs_triangles, &
                           z_t_const,z_soil_center,z_soil_interface, &
                           area_dual,z_vector_dual,wind_old,wind_tend,wind_new, &
                           temperature,wind_div,viscosity_triangles,viscosity,viscosity_rhombi, &
                           condensates_sediment_heat, &
                           totally_first_step_bool) &
  bind(c,name = "manage_pchevi")
    
    real(wp), intent(out) :: wind_new(n_vectors),wind_tend(n_vectors),condensates_sediment_heat(n_scalars), &
                             wind_old(n_vectors),temperature(n_scalars),wind_div(n_scalars), &
                             viscosity_rhombi(n_vectors),viscosity(n_scalars)
    integer,  intent(in)  :: adjacent_signs_h(6*n_scalars_h),adjacent_vector_indices_h(6*n_scalars_h), &
                             totally_first_step_bool,vorticity_indices_triangles(3*n_dual_scalars_h), &
                             vorticity_signs_triangles(3*n_dual_scalars_h)
    real(wp), intent(in)  :: area(n_vectors),latitude_scalar(n_scalars_h),longitude_scalar(n_scalars_h), &
                             area_dual(n_dual_vectors),layer_thickness(n_scalars),z_vector(n_vectors), &
                             z_vector_dual(n_dual_vectors),z_t_const,z_soil_center(nsoillays), &
                             z_soil_interface(nsoillays+1),z_scalar(n_scalars),volume(n_scalars), &
                             viscosity_triangles(n_dual_v_vectors)
    
    ! local variabels
    integer :: h_index,layer_index,vector_index,rk_step
    
    ! Preparations
    ! ------------
    
    ! diagnosing the temperature
    call  temperature_diagnostics(temperature,theta_v_bg,theta_v_pert_old,exner_bg,exner_pert_old,rho_old)
    
    ! updating surface-related turbulence quantities if it is necessary
    if (sfc_sensible_heat_flux==1 .or. sfc_phase_trans==1 .or. pbl_scheme==1) then
      call update_sfc_turb_quantities(is_land,roughness_length,monin_obukhov_length,z_scalar,z_vector, &
                                      theta_v_bg,theta_v_pert_old,v_squared,roughness_velocity,scalar_flux_resistance)
    endif
    
    ! cloud microphysics
    if (MOISTURE_ON==1) then
      call calc_h2otracers_source_rates(rho_old,temperature,layer_thickness,temperature_soil_old, &
                                        phase_trans_rates,phase_trans_heating_rate, &
                                        scalar_flux_resistance,is_land,power_flux_density_latent)
    endif
    
    ! Radiation is updated here.
    if (rad_on>0 .and. rad_update==1) then
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
       call alc_pressure_grad_condensates_v(pressure_gradient_decel_factor,rho_old,gravity_m,pressure_grad_condensates_v)
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





