! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_manage_pchevi

  ! This module manages the predictor-corrector HEVI time stepping.
  
  use mo_definitions,            only: wp,t_grid,t_state,t_diag
  use mo_constituents_nml,       only: n_constituents,lmoist
  use mo_grid_setup,             only: dtime
  use mo_column_solvers,         only: three_band_solver_ver_waves,three_band_solver_gen_densities
  use mo_surface_nml,            only: lsfc_sensible_heat_flux,lsfc_phase_trans,pbl_scheme
  use mo_rad_nml,                only: rad_config
  use mo_p_grad,                 only: manage_p_grad,calc_p_grad_condensates_v
  use mo_derived,                only: temperature_diagnostics
  use mo_pbl,                    only: update_sfc_turb_quantities
  use mo_scalar_tend_expl,       only: scalar_tend_expl
  use mo_vector_tend_expl,       only: vector_tend_expl
  use mo_phase_trans,            only: calc_h2otracers_source_rates
  use mo_manage_radiation_calls, only: update_rad_fluxes
  
  implicit none
  
  contains
  
  subroutine manage_pchevi(state_old,state_new,state_tend,diag,grid,lrad_update,ltotally_first_step,time_coordinate)

    ! This subroutine manages the predictor-corrector HEVI time stepping.
    
    type(t_state), intent(in)    :: state_old           ! state of the old time step
    type(t_state), intent(inout) :: state_new           ! state of the new time step (result)
    type(t_state), intent(inout) :: state_tend          ! state containing tendencies
    type(t_diag),  intent(inout) :: diag                ! diagnostic quantities
    type(t_grid),  intent(inout) :: grid                ! grid properties
    logical,       intent(in)    :: lrad_update         ! radiation update switch
    logical,       intent(in)    :: ltotally_first_step ! true only at the very first step of the whole model run
    real(wp),      intent(in)    :: time_coordinate     ! epoch timestamp of the old time step
    
    ! local variabels
    integer :: rk_step ! Runge-Kutta substep index
    
    ! Preparations
    ! ------------
    
    ! diagnosing the temperature
    call temperature_diagnostics(state_old,diag,grid)
    
    ! updating surface-related turbulence quantities if it is necessary
    if (lsfc_sensible_heat_flux .or. lsfc_phase_trans .or. pbl_scheme==1) then
      call update_sfc_turb_quantities(state_old,diag,grid)
    endif
    
    ! cloud microphysics
    if (lmoist) then
      call calc_h2otracers_source_rates(state_old,diag,grid)
    endif
    
    ! updating radiation if necessary
    if (rad_config>0 .and. lrad_update) then
      call update_rad_fluxes(state_old,diag,grid,time_coordinate)
    endif
      
    ! Loop over the RK substeps
    ! -------------------------
    
    do rk_step=1,2
    
      ! state_old remains unchanged the whole time.
      ! At rk_step==1, state_new contains garbage.
    
      ! 1.) explicit component of the momentum equation
      ! -----------------------------------------------
      ! Update of the pressure gradient.
      if (rk_step==1) then
        call manage_p_grad(state_old,diag,grid,ltotally_first_step)
      endif
      
      ! Only the horizontal momentum is a forward tendency.
      if (rk_step==1) then
        call calc_p_grad_condensates_v(state_old,diag,grid)
        call vector_tend_expl(state_old,state_tend,diag,grid,ltotally_first_step,rk_step)
      endif
      if (rk_step==2) then
        call calc_p_grad_condensates_v(state_new,diag,grid)
        call vector_tend_expl(state_new,state_tend,diag,grid,ltotally_first_step,rk_step)
      endif
      
      ! time stepping for the horizontal momentum can be directly executed
      !$omp parallel workshare
      state_new%wind_h = state_old%wind_h + dtime*state_tend%wind_h
      !$omp end parallel workshare
      ! Horizontal velocity can be considered to be updated from now on.

      ! 2.) explicit component of the generalized density equations
      ! -----------------------------------------------------------
      call scalar_tend_expl(state_old,state_new,state_tend,diag,grid,rk_step)

      ! 3.) vertical sound wave solver
      ! ------------------------------
      call three_band_solver_ver_waves(state_old,state_new,state_tend,diag,grid,rk_step)
      
      ! 4.) vertical tracer advection
      ! -----------------------------
      if (n_constituents>1) then
        call three_band_solver_gen_densities(state_old,state_new,state_tend,diag,grid,rk_step)
      endif
      
    enddo
    
  end subroutine manage_pchevi

end module mo_manage_pchevi





