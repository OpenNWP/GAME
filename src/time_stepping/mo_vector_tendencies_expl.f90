! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_explicit_wind_tend

  ! In this module, the calculation of the explicit part of the momentum equation is managed.
  
  use iso_c_binding
  use mo_definitions,           only: wp
  use grid_nml,                 only: n_vectors_per_layer,n_vectors,n_scalars_h,n_scalars,n_dual_vectors,n_vectors_h, &
                                      n_dual_scalars_h,n_dual_v_vectors,n_layers,n_h_vectors
  use gradient_operators,       only: grad
  use constituents_nml,         only: n_condensed_constituents,n_constituents
  use mo_inner_product,         only: inner_product
  use mo_vorticity_flux,        only: vorticity_flux
  use diff_nml,                 only: lmom_diff_h,lmass_diff_h,ltemp_diff_h,lmom_diff_v
  use surface_nml,              only: pbl_scheme
  use mo_tke,                   only: tke_update
  use momentum_diff_diss,       only: hor_momentum_diffusion,vert_momentum_diffusion,simple_dissipation_rate
  use multiplications,          only: scalar_times_vector
  use planetary_boundary_layer, only: pbl_wind_tendency
  use effective_diff_coeffs,    only: update_n_squared
  use vorticities,              only: calc_pot_vort
  
  implicit none
  
  real(wp), parameter :: impl_thermo_weight = 0.75_wp
  
  contains

  subroutine vector_tendencies_expl(wind,rel_vort_on_triangles,z_vector,z_vector_dual,rel_vort, &
                                    vorticity_indices_triangles,vorticity_signs_triangles,normal_distance, &
                                    area_dual,from_index,to_index,from_index_dual,to_index_dual,inner_product_weights, &
                                    slope,temperature,friction_acc,adjacent_signs_h,adjacent_vector_indices_h,area, &
                                    molecular_diffusion_coeff,normal_distance_dual,rho,tke,viscosity,viscosity_triangles, &
                                    volume,wind_div,viscosity_rhombi,vector_field_placeholder,curl_of_vorticity, &
                                    gravity_m,theta_v_bg,theta_v_pert,scalar_field_placeholder,n_squared,wind_tend, &
                                    density_to_rhombi_indices,density_to_rhombi_weights,dv_hdz,exner_bg, &
                                    exner_pert,f_vec,flux_density,heating_diss,layer_thickness,monin_obukhov_length, &
                                    pot_vort,roughness_length,trsk_indices,trsk_modified_curl_indices,trsk_weights, &
                                    v_squared,vert_hor_viscosity,z_scalar,pot_vort_tend,v_squared_grad, &
                                    pressure_gradient_acc_neg_nl,pressure_gradient_acc_neg_l,pgrad_acc_old, &
                                    pressure_grad_condensates_v,rk_step,totally_first_step_bool) &
  bind(c,name = "vector_tendencies_expl")
  
    real(wp), intent(in)    :: wind(n_vectors),z_vector(n_vectors),z_vector_dual(n_dual_vectors),normal_distance(n_vectors), &
                               area_dual(n_dual_vectors),inner_product_weights(8*n_scalars),slope(n_vectors), &
                               temperature(n_scalars),area(n_vectors),volume(n_scalars),gravity_m(n_vectors), &
                               normal_distance_dual(n_dual_vectors),rho(n_constituents*n_scalars), &
                               theta_v_bg(n_scalars),theta_v_pert(n_scalars),density_to_rhombi_weights(4*n_vectors_h), &
                               exner_bg(n_scalars),exner_pert(n_scalars),f_vec(2*n_vectors_h),layer_thickness(n_scalars), &
                               monin_obukhov_length(n_scalars_h),roughness_length(n_scalars_h),trsk_weights(10*n_vectors_h), &
                               z_scalar(n_scalars),pressure_gradient_acc_neg_nl(n_vectors),pressure_gradient_acc_neg_l(n_vectors), &
                               pgrad_acc_old(n_vectors),pressure_grad_condensates_v(n_vectors)
    integer,  intent(in)    :: vorticity_indices_triangles(3*n_dual_scalars_h),density_to_rhombi_indices(4*n_vectors_h), &
                               vorticity_signs_triangles(3*n_dual_scalars_h),rk_step,totally_first_step_bool, &
                               from_index(n_vectors_h),to_index(n_vectors_h),trsk_indices(10*n_vectors_h), &
                               from_index_dual(n_vectors_h),to_index_dual(n_vectors_h), &
                               adjacent_signs_h(6*n_scalars_h),adjacent_vector_indices_h(6*n_scalars_h), &
                               trsk_modified_curl_indices(10*n_vectors_h)
    real(wp), intent(out)   :: friction_acc(n_vectors),molecular_diffusion_coeff(n_scalars),viscosity(n_scalars), &
                               viscosity_triangles(n_dual_v_vectors),wind_div(n_scalars),viscosity_rhombi(n_vectors), &
                               vector_field_placeholder(n_vectors),rel_vort((2*n_layers+1)*n_vectors_h), &
                               rel_vort_on_triangles(n_dual_v_vectors),curl_of_vorticity(n_vectors),v_squared(n_scalars), &
                               scalar_field_placeholder(n_scalars),n_squared(n_scalars),wind_tend(n_vectors), &
                               dv_hdz(n_h_vectors+n_vectors_h),flux_density(n_vectors),pot_vort((2*n_layers+1)*n_vectors_h), &
                               vert_hor_viscosity(n_h_vectors+n_vectors_h),pot_vort_tend(n_vectors),v_squared_grad(n_vectors)
    real(wp), intent(inout) :: heating_diss(n_scalars),tke(n_scalars)
                               
    ! local variables
    integer  :: ji,layer_index,h_index
    real(wp) :: old_weight,new_weight,old_hor_pgrad_weight,current_hor_pgrad_weight,current_ver_pgrad_weight
  
    ! Managing momentum advection
    ! ---------------------------
    
    if (rk_step==1 .or. totally_first_step_bool==1) then
      call scalar_times_vector(rho(n_condensed_constituents*n_scalars+1:(n_condensed_constituents+1)*n_scalars), &
                               wind,flux_density,from_index,to_index)
      ! Now, the "potential vorticity" is evaluated.
      call calc_pot_vort(wind,rel_vort_on_triangles,z_vector,z_vector_dual,rel_vort, &
                         vorticity_indices_triangles,vorticity_signs_triangles,normal_distance, &
                         area_dual,from_index,to_index,from_index_dual,to_index_dual,inner_product_weights, &
                         slope,f_vec,pot_vort,density_to_rhombi_indices,density_to_rhombi_weights, &
                         rho(n_condensed_constituents*n_scalars))
      ! Now, the generalized Coriolis term is evaluated.
      call vorticity_flux(from_index,to_index,pot_vort_tend,trsk_indices,trsk_modified_curl_indices,trsk_weights, &
                          flux_density,pot_vort,inner_product_weights,adjacent_vector_indices_h)
      ! Kinetic energy is prepared for the gradient term of the Lamb transformation.
      call inner_product(wind,wind,v_squared,adjacent_vector_indices_h,inner_product_weights)
      ! Taking the gradient of the kinetic energy
      call grad(v_squared,v_squared_grad,from_index,to_index,normal_distance,inner_product_weights,slope)
    endif
    
    ! Managing momentum diffusion
    ! ---------------------------
    
    if (rk_step==0) then
      ! updating the Brunt-Väisälä frequency and the TKE if any diffusion is switched on because it is required for computing the diffusion coefficients
      if (lmom_diff_h .or. lmass_diff_h .or. ltemp_diff_h) then
        call update_n_squared(theta_v_bg,theta_v_pert,normal_distance,inner_product_weights,gravity_m, &
                              scalar_field_placeholder,vector_field_placeholder,n_squared)
        call tke_update(rho,viscosity,heating_diss, &
                        tke,vector_field_placeholder,wind,scalar_field_placeholder,from_index,to_index, &
                        adjacent_vector_indices_h,normal_distance,inner_product_weights,slope)
      endif
      
      ! momentum diffusion and dissipation (only updated at the first RK step)
      ! horizontal momentum diffusion
      if (lmom_diff_h) then
        call hor_momentum_diffusion(wind,rel_vort_on_triangles,z_vector,z_vector_dual,rel_vort, &
                                    vorticity_indices_triangles,vorticity_signs_triangles,normal_distance, &
                                    area_dual,from_index,to_index,from_index_dual,to_index_dual,inner_product_weights, &
                                    slope,temperature,friction_acc,adjacent_signs_h,adjacent_vector_indices_h,area, &
                                    molecular_diffusion_coeff,normal_distance_dual,rho,tke,viscosity,viscosity_triangles, &
                                    volume,wind_div,viscosity_rhombi,vector_field_placeholder,curl_of_vorticity)
      endif
      ! vertical momentum diffusion
      if (lmom_diff_v) then
        call vert_momentum_diffusion(wind,z_vector,normal_distance,from_index,to_index,inner_product_weights, &
                                     slope,friction_acc,adjacent_signs_h,adjacent_vector_indices_h,area, &
                                     molecular_diffusion_coeff,rho,tke,viscosity,volume,vector_field_placeholder, &
                                     scalar_field_placeholder,layer_thickness,n_squared,dv_hdz,vert_hor_viscosity)
      endif
      ! planetary boundary layer
      if (pbl_scheme>0) then
        call pbl_wind_tendency(wind,z_vector,monin_obukhov_length,exner_bg,exner_pert,v_squared, &
                               from_index,to_index,friction_acc,gravity_m,roughness_length,rho, &
                               temperature,z_scalar)
      endif
      ! calculation of the dissipative heating rate
      if (lmom_diff_h .or. pbl_scheme>0) then
        call simple_dissipation_rate(wind,friction_acc,heating_diss, &
                                     adjacent_vector_indices_h,inner_product_weights,rho)
      endif
    endif
    
    ! Now the explicit forces are added up.
    new_weight = 1._wp
    if (rk_step==1) then
      new_weight = 0.5_wp
    endif
    old_weight = 1._wp-new_weight
    ! the weights for the pressure gradient
    current_hor_pgrad_weight = 0.5_wp + impl_thermo_weight
    old_hor_pgrad_weight = 1._wp - current_hor_pgrad_weight
    current_ver_pgrad_weight = 1._wp - impl_thermo_weight
    !$omp parallel do private(ji,layer_index,h_index)
    do ji=1,n_vectors
      layer_index = (ji-1)/n_vectors_per_layer
      h_index = ji - layer_index*n_vectors_per_layer
      ! upper and lower boundary
      if (ji<=n_scalars_h .or. ji>=n_vectors-n_scalars_h+1) then
        wind_tend(ji) = 0._wp
      ! horizontal case
      elseif (h_index>=n_scalars_h+1) then
        wind_tend(ji) = &
        old_weight*wind_tend(ji) + new_weight*( &
        ! explicit component of pressure gradient acceleration
        ! old time step component
        old_hor_pgrad_weight*pgrad_acc_old(ji) &
        ! current time step component
        - current_hor_pgrad_weight*(pressure_gradient_acc_neg_nl(ji) + pressure_gradient_acc_neg_l(ji)) &
        ! generalized Coriolis term
        + pot_vort_tend(ji) &
        ! kinetic energy term
        - 0.5_wp*v_squared_grad(ji) &
        ! momentum diffusion
        + friction_acc(ji))
      ! vertical case
      elseif (h_index<=n_scalars_h) then
        wind_tend(ji) = &
        old_weight*wind_tend(ji) + new_weight*( &
        ! explicit component of pressure gradient acceleration
        ! current time step component
        -current_ver_pgrad_weight*(pressure_gradient_acc_neg_nl(ji) + pressure_gradient_acc_neg_l(ji)) &
        ! generalized Coriolis term
        + pot_vort_tend(ji) &
        ! kinetic energy term
        - 0.5_wp*v_squared_grad(ji) &
        ! momentum diffusion
        + friction_acc(ji) &
        ! effect of condensates on the pressure gradient acceleration
        + pressure_grad_condensates_v(ji))
      endif
    enddo
    !$omp end parallel do
  
  end subroutine vector_tendencies_expl

end module mo_explicit_wind_tend
    
    
    
      
    
    
    
    
    
    
