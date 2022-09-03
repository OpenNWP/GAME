! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_vector_tend_expl

  ! In this module, the calculation of the explicit part of the momentum equation is managed.
  
  use mo_definitions,           only: wp,t_grid,t_state,t_diag
  use mo_grid_nml,              only: n_vectors_per_layer,n_vectors,n_scalars_h,n_scalars,n_dual_vectors,n_vectors_h, &
                                      n_dual_scalars_h,n_dual_v_vectors,n_layers,n_h_vectors
  use mo_gradient_operators,    only: grad
  use mo_constituents_nml,      only: n_condensed_constituents,n_constituents
  use mo_inner_product,         only: inner_product
  use mo_vorticity_flux,        only: vorticity_flux
  use mo_diff_nml,              only: lmom_diff_h,lmass_diff_h,ltemp_diff_h,lmom_diff_v
  use mo_surface_nml,           only: pbl_scheme
  use mo_tke,                   only: tke_update
  use mo_momentum_diff_diss,    only: hor_momentum_diffusion,vert_momentum_diffusion,simple_dissipation_rate
  use mo_multiplications,       only: scalar_times_vector
  use mo_pbl,                   only: pbl_wind_tendency
  use mo_effective_diff_coeffs, only: update_n_squared
  use mo_vorticities,           only: calc_pot_vort
  
  implicit none
  
  real(wp), parameter :: impl_thermo_weight = 0.75_wp
  
  contains

  subroutine vector_tend_expl(state,state_tend,totally_first_step_bool,diag,grid,rk_step)
  
    integer,  intent(in)    :: rk_step,totally_first_step_bool
    type(t_state), intent(inout) :: state      ! state to use for calculating the tendencies
    type(t_state), intent(inout) :: state_tend ! state containing the tendencies
    type(t_diag),  intent(inout) :: diag
    type(t_grid),  intent(in)    :: grid
                               
    ! local variables
    integer  :: ji,layer_index,h_index
    real(wp) :: old_weight,new_weight,old_hor_pgrad_weight,current_hor_pgrad_weight,current_ver_pgrad_weight
  
    ! Managing momentum advection
    ! ---------------------------
    
    if (rk_step==1 .or. totally_first_step_bool==1) then
      call scalar_times_vector(state%rho(n_condensed_constituents*n_scalars+1:(n_condensed_constituents+1)*n_scalars), &
                               state%wind,diag%flux_density,grid%from_index,grid%to_index)
      ! Now, the "potential vorticity" is evaluated.
      call calc_pot_vort(state%wind,diag%rel_vort_on_triangles,grid%z_vector,grid%z_vector_dual,diag%rel_vort, &
                         grid%vorticity_indices_triangles,grid%vorticity_signs_triangles,grid%normal_distance, &
                         grid%area_dual,grid%from_index,grid%to_index,grid%from_index_dual,grid%to_index_dual, &
                         grid%inner_product_weights, &
                         grid%slope,grid%f_vec,diag%pot_vort,grid%density_to_rhombi_indices,grid%density_to_rhombi_weights, &
                         state%rho(n_condensed_constituents*n_scalars))
      ! Now, the generalized Coriolis term is evaluated.
      call vorticity_flux(grid%from_index,grid%to_index,diag%pot_vort_tend,grid%trsk_indices,grid%trsk_modified_curl_indices, &
                          grid%trsk_weights, &
                          diag%flux_density,diag%pot_vort,grid%inner_product_weights,grid%adjacent_vector_indices_h)
      ! Kinetic energy is prepared for the gradient term of the Lamb transformation.
      call inner_product(state%wind,state%wind,diag%v_squared,grid)
      ! Taking the gradient of the kinetic energy
      call grad(diag%v_squared,diag%v_squared_grad,grid%from_index,grid%to_index, &
                grid%normal_distance,grid%inner_product_weights,grid%slope)
    endif
    
    ! Managing momentum diffusion
    ! ---------------------------
    
    if (rk_step==0) then
      ! updating the Brunt-Väisälä frequency and the TKE if any diffusion is switched on because it is required for computing the diffusion coefficients
      if (lmom_diff_h .or. lmass_diff_h .or. ltemp_diff_h) then
        call update_n_squared(grid%theta_v_bg,state%theta_v_pert,grid%normal_distance,grid%inner_product_weights,grid%gravity_m, &
                              diag%scalar_field_placeholder,diag%vector_field_placeholder,diag%n_squared)
        call tke_update(state%rho,diag%viscosity,diag%heating_diss, &
                        diag%tke,diag%vector_field_placeholder,state%wind,diag%scalar_field_placeholder,grid)
      endif
      
      ! momentum diffusion and dissipation (only updated at the first RK step)
      ! horizontal momentum diffusion
      if (lmom_diff_h) then
        call hor_momentum_diffusion(state%wind,diag%rel_vort_on_triangles,grid%z_vector,grid%z_vector_dual,diag%rel_vort, &
                                    grid%vorticity_indices_triangles,grid%vorticity_signs_triangles,grid%normal_distance, &
                                    grid%area_dual,grid%from_index,grid%to_index,grid%from_index_dual,grid%to_index_dual, &
                                    grid%inner_product_weights, &
                                    grid%slope,diag%temperature,diag%friction_acc,grid%adjacent_signs_h, &
                                    grid%adjacent_vector_indices_h,grid%area, &
                                    diag%molecular_diffusion_coeff,grid%normal_distance_dual,state%rho,diag%tke,diag%viscosity, &
                                    diag%viscosity_triangles, &
                                    grid%volume,diag%wind_div,diag%viscosity_rhombi,diag%vector_field_placeholder, &
                                    diag%curl_of_vorticity,grid)
      endif
      ! vertical momentum diffusion
      if (lmom_diff_v) then
        call vert_momentum_diffusion(state%wind,grid%z_vector,grid%normal_distance,grid%from_index,grid%to_index, &
                                     grid%inner_product_weights, &
                                     grid%slope,diag%friction_acc,grid%adjacent_signs_h,grid%adjacent_vector_indices_h,grid%area, &
                                     diag%molecular_diffusion_coeff,state%rho,diag%tke,diag%viscosity,grid%volume, &
                                     diag%vector_field_placeholder, &
                                     diag%scalar_field_placeholder,grid%layer_thickness, &
                                     diag%n_squared,diag%dv_hdz,diag%vert_hor_viscosity,grid)
      endif
      ! planetary boundary layer
      if (pbl_scheme>0) then
        call pbl_wind_tendency(state%wind,grid%z_vector,diag%monin_obukhov_length,grid%exner_bg,state%exner_pert,diag%v_squared, &
                               grid%from_index,grid%to_index,diag%friction_acc,grid%gravity_m,grid%roughness_length,state%rho, &
                               diag%temperature,grid%z_scalar)
      endif
      ! calculation of the dissipative heating rate
      if (lmom_diff_h .or. pbl_scheme>0) then
        call simple_dissipation_rate(state,diag%friction_acc,diag%heating_diss,grid)
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
        state_tend%wind(ji) = 0._wp
      ! horizontal case
      elseif (h_index>=n_scalars_h+1) then
        state_tend%wind(ji) = &
        old_weight*state_tend%wind(ji) + new_weight*( &
        ! explicit component of pressure gradient acceleration
        ! old time step component
        old_hor_pgrad_weight*diag%pgrad_acc_old(ji) &
        ! current time step component
        - current_hor_pgrad_weight*(diag%pressure_gradient_acc_neg_nl(ji) + diag%pressure_gradient_acc_neg_l(ji)) &
        ! generalized Coriolis term
        + diag%pot_vort_tend(ji) &
        ! kinetic energy term
        - 0.5_wp*diag%v_squared_grad(ji) &
        ! momentum diffusion
        + diag%friction_acc(ji))
      ! vertical case
      elseif (h_index<=n_scalars_h) then
        state_tend%wind(ji) = &
        old_weight*state_tend%wind(ji) + new_weight*( &
        ! explicit component of pressure gradient acceleration
        ! current time step component
        -current_ver_pgrad_weight*(diag%pressure_gradient_acc_neg_nl(ji) + diag%pressure_gradient_acc_neg_l(ji)) &
        ! generalized Coriolis term
        + diag%pot_vort_tend(ji) &
        ! kinetic energy term
        - 0.5_wp*diag%v_squared_grad(ji) &
        ! momentum diffusion
        + diag%friction_acc(ji) &
        ! effect of condensates on the pressure gradient acceleration
        + diag%pressure_grad_condensates_v(ji))
      endif
    enddo
    !$omp end parallel do
  
  end subroutine vector_tend_expl

end module mo_vector_tend_expl
    
    
    
      
    
    
    
    
    
    
