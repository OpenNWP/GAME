! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_manage_pchevi

  ! This module manages the predictor-corrector HEVI time stepping.
  
  use mo_definitions,            only: wp,t_grid,t_state,t_diag
  use mo_constituents_nml,       only: n_constituents,lmoist,n_condensed_constituents
  use mo_grid_nml,               only: n_layers,n_vectors_per_layer,n_cells,n_vectors,n_dual_vectors,n_scalars, &
                                       n_dual_scalars_h,n_dual_v_vectors,n_h_vectors,n_edges
  use mo_grid_setup,             only: dtime
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
  
  subroutine manage_pchevi(state_old,state_new,state_tend,diag,grid,lrad_update,ltotally_first_step,time_coordinate)

    ! This subroutine manages the predictor-corrector HEVI time stepping.
    
    type(t_grid),  intent(inout) :: grid                ! grid properties
    type(t_state), intent(inout) :: state_old           ! state of the old time step
    type(t_state), intent(inout) :: state_new           ! state of the new time step (result)
    type(t_state), intent(inout) :: state_tend          ! state containing tendencies
    type(t_diag),  intent(inout) :: diag                ! diagnostic quantities
    logical,       intent(in)    :: lrad_update         ! radiation update switch
    logical,       intent(in)    :: ltotally_first_step ! true only at the very first step of the whole model run
    real(wp),      intent(in)    :: time_coordinate     ! epoch timestamp of the old time step
    
    ! local variabels
    integer :: ji,jl,vector_index,rk_step
    
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
       call calc_pressure_grad_condensates_v(state_old,diag,grid)
        ! Only the horizontal momentum is a forward tendency.
       call  vector_tend_expl(state_old,state_tend,diag,grid,ltotally_first_step,rk_step)
      endif
      if (rk_step==2) then
        call calc_pressure_grad_condensates_v(state_new,diag,grid)
        ! Only the horizontal momentum is a forward tendency.
        call vector_tend_expl(state_new,state_tend,diag,grid,ltotally_first_step,rk_step)
      endif
      
      ! time stepping for the horizontal momentum can be directly executed
      !$omp parallel do private(ji,jl,vector_index)
      do ji=1,n_edges
        do jl=0,n_layers-1
          vector_index = n_cells + jl*n_vectors_per_layer + ji
          state_new%wind(vector_index) = state_old%wind(vector_index) + dtime*state_tend%wind(vector_index)
        enddo
      enddo
      !$omp end parallel do
      ! Horizontal velocity can be considered to be updated from now on.

      ! 2.) explicit component of the generalized density equations
      ! -----------------------------------------------------------
      if (rk_step==1) then
        call scalar_tend_expl(state_old,state_new,state_tend,diag,grid,rk_step)
      endif
      if (rk_step==2) then
        call scalar_tend_expl(state_new,state_new,state_tend,diag,grid,rk_step)
      endif

      ! 3.) vertical sound wave solver
      ! ------------------------------
      if (rk_step==1) then
        call three_band_solver_ver_waves(state_old,state_old,state_new,state_tend,diag,grid,rk_step)
      endif
      if (rk_step==2) then
        call three_band_solver_ver_waves(state_old,state_new,state_new,state_tend,diag,grid,rk_step)
      endif
      
      ! 4.) vertical tracer advection
      ! -----------------------------
      if (n_constituents>1) then
        call three_band_solver_gen_densities(state_old,state_new,state_tend,diag,grid,rk_step)
      endif
      
    enddo
    
  end subroutine manage_pchevi

end module mo_manage_pchevi





