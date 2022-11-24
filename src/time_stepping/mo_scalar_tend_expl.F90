! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME
  
module mo_scalar_tend_expl
  
  ! This module contains the horizontal (explicit) part of the constituent integration.

  use mo_definitions,        only: wp,t_grid,t_state,t_diag
  use mo_constants,          only: c_d_v,c_d_p
  use mo_grid_nml,           only: n_cells,n_layers
  use mo_constituents_nml,   only: n_constituents,n_condensed_constituents
  use mo_derived,            only: c_v_mass_weighted_air
  use mo_diff_nml,           only: lmass_diff_h,lmass_diff_v,ltemp_diff_h,ltemp_diff_v
  use mo_multiplications,    only: scalar_times_vector_h,scalar_times_vector_v,scalar_times_vector_h_upstream
  use mo_divergences,        only: div_h_tracer,div_h,add_vertical_div
  use mo_gradient_operators, only: grad_hor,grad_vert
  use mo_eff_diff_coeffs,    only: scalar_diffusion_coeffs

  implicit none
  
  contains

  subroutine scalar_tend_expl(state_scalar,state_wind,state_tend,diag,grid,rk_step)
    
    ! This subroutine manages the calculation of the explicit part of the scalar tendencies.
    
    type(t_state), intent(in)    :: state_scalar ! state containing the scalar quantities
    type(t_state), intent(in)    :: state_wind   ! state from which to use the wind
    type(t_state), intent(inout) :: state_tend   ! state containing the tendencies
    type(t_diag),  intent(inout) :: diag         ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid         ! grid quantities
    integer,       intent(in)    :: rk_step      ! Runge-Kutta substep
    
    ! local variables
    integer  :: ji                         ! cell index
    integer  :: jl                         ! layer index
    integer  :: jc                         ! constituent index
    real(wp) :: old_weight(n_constituents) ! time stepping weight of the old time step
    real(wp) :: new_weight(n_constituents) ! time stepping weight of the new time step
    
    ! Firstly,some things need to prepared.
    ! -------------------------------------
    
    ! determining the time stepping weights
    do jc=1,n_constituents
      new_weight(jc) = 1._wp
      if (rk_step==2 .and. jc/=n_condensed_constituents+1) then
        new_weight(jc) = 0.5_wp
      endif
      old_weight(jc) = 1._wp-new_weight(jc)
    enddo
    
    ! updating the scalar diffusion coefficient if required
    if (rk_step==1 .and. (lmass_diff_h .or. ltemp_diff_h)) then
      call scalar_diffusion_coeffs(state_scalar,diag,grid)
    endif
    
    ! Temperature diffusion gets updated at the first RK step if required.
    if (ltemp_diff_h .and. rk_step==1) then
      ! The diffusion of the temperature depends on its gradient.
      call grad_vert(diag%temperature,diag%vector_placeholder_v,grid)
      call grad_hor(diag%temperature,diag%vector_placeholder_h,diag%vector_placeholder_v,grid)
      ! Now the diffusive temperature flux density can be obtained.
      call scalar_times_vector_h(diag%temp_diffusion_coeff_numerical_h,diag%vector_placeholder_h,diag%flux_density_h,grid)
      ! The divergence of the diffusive temperature flux density is the diffusive temperature heating.
      call div_h(diag%flux_density_h,diag%temperature_diffusion_heating,grid)
      ! vertical temperature diffusion
      if (ltemp_diff_v) then
        call scalar_times_vector_v(diag%temp_diffusion_coeff_numerical_v,diag%vector_placeholder_v,diag%flux_density_v)
        call add_vertical_div(diag%flux_density_v,diag%temperature_diffusion_heating,grid)
      endif
    endif
    
    ! Mass diffusion gets updated at the first RK step if required.
    if (lmass_diff_h .and. rk_step==1) then
      ! loop over all constituents
      do jc=1,n_constituents
        ! The diffusion of the tracer density depends on its gradient.
        call grad_vert(state_scalar%rho(:,:,jc),diag%vector_placeholder_v,grid)
        call grad_hor(state_scalar%rho(:,:,jc),diag%vector_placeholder_h,diag%vector_placeholder_v,grid)
        ! Now the diffusive mass flux density can be obtained.
        call scalar_times_vector_h(diag%mass_diffusion_coeff_numerical_h,diag%vector_placeholder_h,diag%flux_density_h,grid)
        ! The divergence of the diffusive mass flux density is the diffusive mass source rate.
        call div_h(diag%flux_density_h,diag%mass_diff_tendency(:,:,jc),grid)
        ! vertical mass diffusion
        if (lmass_diff_v) then
          call scalar_times_vector_v(diag%mass_diffusion_coeff_numerical_v,diag%vector_placeholder_v,diag%flux_density_v)
          call add_vertical_div(diag%flux_density_v,diag%mass_diff_tendency(:,:,jc),grid)
        endif
      enddo
    endif
    
    ! Now,the actual scalar tendencies can be computed.
    ! -------------------------------------------------
    
    ! loop over all constituents
    do jc=1,n_constituents
      
      ! This is the mass advection,which needs to be carried out for all constituents.
      ! ------------------------------------------------------------------------------
      ! moist air
      if (jc==n_condensed_constituents+1) then
        call scalar_times_vector_h(state_scalar%rho(:,:,jc),state_wind%wind_h,diag%flux_density_h,grid)
        call div_h(diag%flux_density_h,diag%flux_density_div,grid)
      ! all other constituents
      else
        call scalar_times_vector_h_upstream(state_scalar%rho(:,:,jc),state_wind%wind_h,diag%flux_density_h,grid)
        call div_h_tracer(diag%flux_density_h,state_scalar%rho(:,:,jc),state_wind%wind_h,diag%flux_density_div,grid)
      endif
      
      ! adding the tendencies in all grid boxes
      !$omp parallel workshare
      state_tend%rho(:,:,jc) = old_weight(jc)*state_tend%rho(:,:,jc) + new_weight(jc)*( &
      ! the advection
      - diag%flux_density_div &
      ! the diffusion
      + diag%mass_diff_tendency(:,:,jc) &
      ! phase transitions
      + diag%phase_trans_rates(:,:,min(jc,n_condensed_constituents+1)))
      !$omp end parallel workshare
      
      ! Explicit component of the rho*theta_v integration
      ! -------------------------------------------------
      
      if (jc==n_condensed_constituents+1) then
      
        ! determining the virtual potential temperature
        !$omp parallel workshare
        diag%scalar_placeholder = state_scalar%rhotheta_v/state_scalar%rho(:,:,jc)
        !$omp end parallel workshare
        
        ! calculating the divergence of the horizontal virtual potential temperature flux density
        call scalar_times_vector_h(diag%scalar_placeholder,diag%flux_density_h,diag%vector_placeholder_h,grid)
        call div_h(diag%vector_placeholder_h,diag%flux_density_div,grid)
        
        ! adding the tendencies in all grid boxes
        !$omp parallel do private(ji,jl)
        do jl=1,n_layers
          do ji=1,n_cells
            state_tend%rhotheta_v(ji,jl) = old_weight(jc)*state_tend%rhotheta_v(ji,jl) &
            + new_weight(jc)*( &
            ! the advection (resolved transport)
            -diag%flux_density_div(ji,jl) &
            ! the diabatic forcings
            ! weighting factor accounting for condensates
            + c_d_v*state_scalar%rho(ji,jl,jc)/c_v_mass_weighted_air(state_scalar%rho,diag%temperature,ji,jl)*( &
            ! dissipation through molecular + turbulent momentum diffusion
            diag%heating_diss(ji,jl) &
            ! molecular + turbulent heat transport
            + diag%temperature_diffusion_heating(ji,jl) &
            ! radiation
            + diag%radiation_tendency(ji,jl) &
            ! phase transitions
            + diag%phase_trans_heating_rate(ji,jl) &
            ! heating rate due to falling condensates
            + diag%condensates_sediment_heat(ji,jl)) &
            ! this has to be divided by c_p*exner
            /(c_d_p*(grid%exner_bg(ji,jl) + state_scalar%exner_pert(ji,jl))) &
            ! tendency of due to phase transitions and mass diffusion
            + (diag%phase_trans_rates(ji,jl,jc) + diag%mass_diff_tendency(ji,jl,jc)) &
            *diag%scalar_placeholder(ji,jl))
          enddo
        enddo
        !$omp end parallel do
      endif
    enddo
  
  end subroutine scalar_tend_expl

end module mo_scalar_tend_expl
















