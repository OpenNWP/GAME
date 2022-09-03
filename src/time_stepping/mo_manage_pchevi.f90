! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_manage_pchevi

  ! This module manages the predictor-corrector HEVI time stepping.
  
  use mo_definitions,            only: wp,t_grid,t_state,t_diag
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
  use mo_manage_radiation_calls, only: update_rad_fluxes
  
  implicit none
  
  contains
  
  subroutine manage_pchevi(state_old,state_new,state_tend,ltotally_first_step,diag,grid,lrad_update,time_coordinate)

    ! This subroutine manages the predictor-corrector HEVI time stepping.
    
    type(t_grid),  intent(inout) :: grid
    type(t_state), intent(inout) :: state_old
    type(t_state), intent(inout) :: state_new
    type(t_state), intent(inout) :: state_tend
    type(t_diag),  intent(inout) :: diag
    logical,       intent(in)    :: ltotally_first_step,lrad_update
    real(wp),      intent(in)    :: time_coordinate
    
    ! local variabels
    integer :: h_index,layer_index,vector_index,rk_step
    
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
    
    ! Radiation is updated here.
    if (rad_config>0 .and. lrad_update) then
      call update_rad_fluxes(state_old,diag,grid,time_coordinate)
    endif
      
    ! Loop over the RK substeps
    ! -------------------------
    
    do rk_step=1,2
    
      ! state_old remains unchanged the whole time.
      ! At rk_step == 1, state_new contains garbage.
    
      ! 1.) explicit component of the momentum equation
      ! -----------------------------------------------
      ! Update of the pressure gradient.
      if (rk_step==1) then
        call manage_pressure_gradient(state_old,diag,grid,ltotally_first_step)
      endif
      
      if (rk_step==1) then
       call calc_pressure_grad_condensates_v(diag%pressure_gradient_decel_factor, &
                                             state_old%rho,grid%gravity_m,diag%pressure_grad_condensates_v)
        ! Only the horizontal momentum is a forward tendency.
       call  vector_tend_expl(state_old,state_tend,diag,grid,ltotally_first_step,rk_step)
      endif
      if (rk_step==2) then
        call calc_pressure_grad_condensates_v(diag%pressure_gradient_decel_factor, &
                                              state_new%rho,grid%gravity_m,diag%pressure_grad_condensates_v)
        ! Only the horizontal momentum is a forward tendency.
        call vector_tend_expl(state_new,state_tend,diag,grid,ltotally_first_step,rk_step)
      endif
      
      ! time stepping for the horizontal momentum can be directly executed
      !$omp parallel do private(h_index,layer_index,vector_index)
      do h_index=1,n_vectors_h
        do layer_index=0,n_layers-1
          vector_index = n_scalars_h + layer_index*n_vectors_per_layer + h_index
          state_new%wind(vector_index) = state_old%wind(vector_index) + dtime*state_tend%wind(vector_index)
        enddo
      enddo
      !$omp end parallel do
      ! Horizontal velocity can be considered to be updated from now on.

      ! 2.) explicit component of the generalized density equations
      ! -----------------------------------------------------------
      if (rk_step==1) then
        call scalar_tend_expl(state_old%rho,state_tend%rhotheta_v,state_old%rhotheta_v, &
                              state_new%wind,state_tend%rho,state_old%exner_pert, &
                              diag,grid,rk_step)
      endif
      if (rk_step==2) then
        call scalar_tend_expl(state_new%rho,state_tend%rhotheta_v,state_new%rhotheta_v, &
                              state_new%wind,state_tend%rho,state_new%exner_pert, &
                              diag,grid,rk_step)
      endif

      ! 3.) vertical sound wave solver
      ! ------------------------------
      if (rk_step==1) then
        call three_band_solver_ver_waves(state_new%theta_v_pert, &
                                         state_new%rho,state_old%rhotheta_v,state_old%rho, &
                                         state_old%rho,state_old%theta_v_pert, &
                                         state_old%exner_pert,state_old%theta_v_pert,state_old%exner_pert, &
                                         state_old%temperature_soil,state_old%temperature_soil, &
                                         state_tend%rhotheta_v,state_tend%rho,state_old%rhotheta_v,state_old%wind, &
                                         state_tend%wind, &
                                         state_new%rhotheta_v,state_new%exner_pert, &
                                         state_new%wind,state_new%temperature_soil,diag,grid,rk_step)
      endif
      if (rk_step==2) then
        call three_band_solver_ver_waves(state_new%theta_v_pert, &
                                         state_new%rho,state_new%rhotheta_v,state_new%rho, &
                                         state_old%rho,state_old%theta_v_pert, &
                                         state_old%exner_pert,state_new%theta_v_pert,state_new%exner_pert, &
                                         state_new%temperature_soil, &
                                         state_old%temperature_soil, &
                                         state_tend%rhotheta_v,state_tend%rho,state_old%rhotheta_v,state_old%wind, &
                                         state_tend%wind, &
                                         state_new%rhotheta_v,state_new%exner_pert, &
                                         state_new%wind,state_new%temperature_soil,diag,grid,rk_step)
      endif
      
      ! 4.) vertical tracer advection
      ! -----------------------------
      if (n_constituents>1) then
        call three_band_solver_gen_densities(state_old%wind,state_new%wind,state_tend%rho,state_old%rho,state_new%rho, &
                                             diag%condensates_sediment_heat,diag%temperature,grid,rk_step)
      endif
      
    enddo
    
  end subroutine manage_pchevi

end module mo_manage_pchevi





