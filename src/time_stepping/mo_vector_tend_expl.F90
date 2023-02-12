! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_vector_tend_expl
  
  ! In this module the calculation of the explicit part of the momentum equation is managed.
  
  use mo_definitions,        only: wp,t_grid,t_state,t_diag
  use mo_constants,          only: impl_thermo_weight
  use mo_gradient_operators, only: grad_hor,grad_vert
  use mo_constituents_nml,   only: n_condensed_constituents
  use mo_inner_product,      only: inner_product
  use mo_vorticity_flux,     only: vorticity_flux
  use mo_run_nml,            only: luse_bg_state
  use mo_diff_nml,           only: lmom_diff_h,lmass_diff_h,ltemp_diff_h,lmom_diff_v,lklemp
  use mo_surface_nml,        only: pbl_scheme
  use mo_tke,                only: tke_update
  use mo_momentum_diff_diss, only: mom_diff_h,mom_diff_v,simple_dissipation_rate
  use mo_multiplications,    only: scalar_times_vector_h,scalar_times_vector_v
  use mo_pbl,                only: pbl_wind_tendency
  use mo_eff_diff_coeffs,    only: update_n_squared
  use mo_vorticities,        only: calc_pot_vort
  
  implicit none
  
  contains
  
  subroutine vector_tend_expl(state,state_tend,diag,grid,ltotally_first_step,rk_step)
    
    ! This subroutine manages the calculation of the explicit part of the wind tendencies.
    
    type(t_state), intent(in)    :: state               ! state to use for calculating the tendencies
    type(t_state), intent(inout) :: state_tend          ! state containing the tendencies
    type(t_diag),  intent(inout) :: diag                ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid                ! grid quantities
    integer,       intent(in)    :: rk_step             ! Runge-Kutta substep
    logical,       intent(in)    :: ltotally_first_step ! switch for the very first step of the whole model run
    
    ! local variables
    real(wp) :: old_weight               ! weight of the old predictor-corrector substep
    real(wp) :: new_weight               ! weight of the new predictor-corrector substep
    integer  :: no_bg_switch             ! set to one if using no hydrostatic background state, zero otherwise
    real(wp) :: old_hor_pgrad_weight     ! old time step horizontal pressure gradient weight
    real(wp) :: current_hor_pgrad_weight ! current time step horizontal pressure gradient weight
    real(wp) :: current_ver_pgrad_weight ! current time step vertical pressure gradient weight
    
    no_bg_switch = 0
    if (.not. luse_bg_state) then
      no_bg_switch = 1
    endif
    
    ! Managing momentum advection
    ! ---------------------------
    
    ! Momentum advection is only computed at the second predictor-corrector substep.
    if (rk_step==2 .or. ltotally_first_step) then
      call scalar_times_vector_h(state%rho(:,:,n_condensed_constituents+1),state%wind_h,diag%flux_density_h,grid)
      call scalar_times_vector_v(state%rho(:,:,n_condensed_constituents+1),state%wind_v,diag%flux_density_v)
      ! Now, the "potential vorticity" is evaluated.
      call calc_pot_vort(state,state%rho(:,:,n_condensed_constituents+1),diag,grid)
      ! Now, the generalized Coriolis term is evaluated.
      call vorticity_flux(diag,grid)
      ! Kinetic energy is prepared for the gradient term of the Lamb transformation.
      call inner_product(state%wind_h,state%wind_v,state%wind_h,state%wind_v,diag%v_squared,grid)
      ! computing the gradient of the kinetic energy
      call grad_vert(diag%v_squared,diag%v_squared_grad_v,grid)
      call grad_hor(diag%v_squared,diag%v_squared_grad_h,diag%v_squared_grad_v,grid)
    endif
    
    ! Managing momentum diffusion
    ! ---------------------------
    
    ! momentum diffusion and dissipation (only updated at the first RK step)
    if (rk_step==1) then
      ! updating the Brunt-Väisälä frequency and the TKE if any diffusion is switched on because it is required for computing the diffusion coefficients
      if (lmom_diff_h .or. lmass_diff_h .or. ltemp_diff_h) then
        call update_n_squared(state,diag,grid)
        call tke_update(state,diag,grid)
      endif
      
      ! momentum diffusion and dissipation (only updated at the first RK step)
      ! horizontal momentum diffusion
      if (lmom_diff_h) then
        call mom_diff_h(state,diag,grid)
      endif
      ! vertical momentum diffusion
      if (lmom_diff_v) then
        call mom_diff_v(state,diag,grid)
      endif
      ! planetary boundary layer
      if (pbl_scheme>0) then
        call pbl_wind_tendency(state,diag,grid)
      endif
      ! calculation of the dissipative heating rate
      if (lmom_diff_h .or. pbl_scheme>0 .or. lklemp) then
        call simple_dissipation_rate(state,diag,grid)
      endif
    endif
    
    ! Now the explicit forces are added up.
    ! the weights for momentum advection
    new_weight = 1._wp
    if (rk_step==2) then
      new_weight = 0.5_wp
    endif
    old_weight = 1._wp-new_weight
    ! the weights for the pressure gradient
    current_hor_pgrad_weight = 0.5_wp + impl_thermo_weight
    old_hor_pgrad_weight = 1._wp - current_hor_pgrad_weight
    current_ver_pgrad_weight = 1._wp - impl_thermo_weight
    
    !$omp parallel workshare
    
    ! horizontal case
    state_tend%wind_h = old_weight*state_tend%wind_h + new_weight*( &
    ! explicit component of pressure gradient acceleration
    ! old time step component
    old_hor_pgrad_weight*diag%pgrad_acc_old_h &
    ! current time step component
    - current_hor_pgrad_weight*(diag%pressure_gradient_acc_neg_nl_h + diag%pressure_gradient_acc_neg_l_h) &
    ! generalized Coriolis term
    + diag%pot_vort_tend_h &
    ! kinetic energy term
    - 0.5_wp*diag%v_squared_grad_h &
    ! momentum diffusion
    + diag%friction_acc_h)
    
    ! vertical case
    state_tend%wind_v = old_weight*state_tend%wind_v + new_weight*( &
    ! explicit component of pressure gradient acceleration
    ! current time step component
    - current_ver_pgrad_weight*(diag%pressure_gradient_acc_neg_nl_v + diag%pressure_gradient_acc_neg_l_v) &
    - no_bg_switch*(1._wp-current_ver_pgrad_weight)*diag%pressure_gradient_acc_neg_l_v &
    ! generalized Coriolis term
    + diag%pot_vort_tend_v &
    ! kinetic energy term
    - 0.5_wp*diag%v_squared_grad_v &
    ! momentum diffusion
    + diag%friction_acc_v &
    ! effect of condensates on the pressure gradient acceleration
    + diag%pressure_grad_condensates_v)
    !$omp end parallel workshare
  
  end subroutine vector_tend_expl

end module mo_vector_tend_expl














