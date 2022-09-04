! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME
  
module mo_scalar_tend_expl
  
  ! This module contains the horizontal (explicit) part of the constituent integration.

  use mo_definitions,           only: wp,t_grid,t_state,t_diag
  use mo_constants,             only: c_d_v,c_d_p
  use mo_grid_nml,              only: n_scalars,n_scalars_h,n_vectors,n_vectors_h,n_dual_v_vectors,n_dual_scalars_h
  use mo_constituents_nml,      only: n_constituents,n_condensed_constituents,lmoist
  use mo_derived,               only: c_v_mass_weighted_air
  use mo_diff_nml,              only: lmass_diff_h,lmass_diff_v,ltemp_diff_h,ltemp_diff_v
  use mo_multiplications,       only: scalar_times_vector_h,scalar_times_vector_v,scalar_times_vector_h_upstream
  use mo_divergences,           only: div_h_tracer,div_h,add_vertical_div
  use mo_gradient_operators,    only: grad
  use mo_effective_diff_coeffs, only: scalar_diffusion_coeffs

  implicit none
  
  contains

  subroutine scalar_tend_expl(state_scalar,state_wind,state_tend,diag,grid,rk_step)
    
    ! This subroutine manages the calculation of the explicit part of the scalar tendencies.
    
    type(t_state), intent(in)    :: state_scalar ! state containing the scalar quantities
    type(t_state), intent(in)    :: state_wind   ! state from which to use the wind
    type(t_state), intent(inout) :: state_tend   ! state containing the tendencies
    type(t_diag),  intent(inout) :: diag         ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid         ! grid quantities
    integer,       intent(in)    :: rk_step      ! Runge-Kutta step
    
    ! local variables
    integer  :: ji,jc,scalar_shift_index,scalar_shift_index_phase_trans,scalar_index
    real(wp) :: old_weight(n_constituents),new_weight(n_constituents)
    
    ! Firstly,some things need to prepared.
    ! --------------------------------------
    
    ! determining the RK weights
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
      call grad(diag%temperature,diag%vector_field_placeholder,grid%from_index,grid%to_index, &
                grid%normal_distance,grid%inner_product_weights,grid%slope)
      ! Now the diffusive temperature flux density can be obtained.
      call scalar_times_vector_h(diag%temp_diffusion_coeff_numerical_h,diag%vector_field_placeholder, &
                                 diag%flux_density,grid)
      ! The divergence of the diffusive temperature flux density is the diffusive temperature heating.
      call div_h(diag%flux_density,diag%temperature_diffusion_heating, &
                 grid%adjacent_signs_h,grid%adjacent_vector_indices_h,grid%inner_product_weights,grid%slope,grid%area,grid%volume)
      ! vertical temperature diffusion
      if (ltemp_diff_v) then
        call scalar_times_vector_v(diag%temp_diffusion_coeff_numerical_v,diag%vector_field_placeholder,diag%flux_density)
        call add_vertical_div(diag%flux_density,diag%temperature_diffusion_heating,grid%area,grid%volume)
      endif
    endif
    
    ! Mass diffusion gets updated at the first RK step if required.
    if (lmass_diff_h .and. rk_step==1) then
      ! loop over all constituents
      do jc=1,n_constituents
        scalar_shift_index = (jc-1)*n_scalars

        ! The diffusion of the tracer density depends on its gradient.
        call grad(state_scalar%rho(scalar_shift_index+1:scalar_shift_index+n_scalars),diag%vector_field_placeholder, &
                  grid%from_index,grid%to_index,grid%normal_distance,grid%inner_product_weights,grid%slope)
        ! Now the diffusive mass flux density can be obtained.
        call scalar_times_vector_h(diag%mass_diffusion_coeff_numerical_h, &
                                   diag%vector_field_placeholder,diag%vector_field_placeholder,grid)
        ! The divergence of the diffusive mass flux density is the diffusive mass source rate.
        call div_h(diag%vector_field_placeholder,diag%mass_diff_tendency(scalar_shift_index+1:scalar_shift_index+n_scalars), &
                   grid%adjacent_signs_h,grid%adjacent_vector_indices_h,grid%inner_product_weights,grid%slope,grid%area,grid%volume)
        ! vertical mass diffusion
        if (lmass_diff_v) then
          call scalar_times_vector_v(diag%mass_diffusion_coeff_numerical_v, &
                                     diag%vector_field_placeholder,diag%vector_field_placeholder)
          call add_vertical_div(diag%vector_field_placeholder, &
                                diag%mass_diff_tendency(scalar_shift_index+1:scalar_shift_index+n_scalars), &
                                grid%area,grid%volume)
        endif
      enddo
    endif
    
    ! Now,the actual scalar tendencies can be computed.
    ! --------------------------------------------------
    
    ! loop over all constituents
    do jc=1,n_constituents
      scalar_shift_index = (jc-1)*n_scalars
      scalar_shift_index_phase_trans = scalar_shift_index
      if (lmoist .and. jc==n_condensed_constituents+2) then
        scalar_shift_index_phase_trans = scalar_shift_index - n_scalars
      endif
      
        ! This is the mass advection,which needs to be carried out for all constituents.
        ! -------------------------------------------------------------------------------
        ! moist air
      if (jc==n_condensed_constituents+1) then
        call scalar_times_vector_h(state_scalar%rho(scalar_shift_index+1:scalar_shift_index+n_scalars), &
                                   state_wind%wind,diag%flux_density,grid)
        call div_h(diag%flux_density,diag%flux_density_div, &
                   grid%adjacent_signs_h,grid%adjacent_vector_indices_h,grid%inner_product_weights,grid%slope,grid%area,grid%volume)
      ! all other constituents
      else
        call scalar_times_vector_h_upstream(state_scalar%rho(scalar_shift_index+1:scalar_shift_index+n_scalars), &
                                            state_wind%wind,diag%flux_density,grid)
        call div_h_tracer(diag%flux_density,state_scalar%rho(scalar_shift_index+1:scalar_shift_index+n_scalars), &
                          state_wind%wind,diag%flux_density_div, &
                          grid%adjacent_signs_h,grid%adjacent_vector_indices_h, &
                          grid%inner_product_weights,grid%slope,grid%area,grid%volume)
      endif
      
      ! adding the tendencies in all grid boxes
      !$omp parallel do private(ji,scalar_index)
      do ji=1,n_scalars
        scalar_index = scalar_shift_index + ji
        state_tend%rho(scalar_index) &
        = old_weight(jc)*state_tend%rho(scalar_index) &
        + new_weight(jc)*( &
        ! the advection
        - diag%flux_density_div(ji) &
        ! the diffusion
        + diag%mass_diff_tendency(scalar_shift_index+ji) &
        ! phase transitions
        + diag%phase_trans_rates(scalar_shift_index_phase_trans+ji))
      enddo
      !$omp end parallel do
      
      ! Explicit component of the rho*theta_v integration
      ! -------------------------------------------------
      
      if (jc==n_condensed_constituents+1) then
        ! determining the virtual potential temperature
        !$omp parallel workshare
        diag%scalar_field_placeholder = state_scalar%rhotheta_v/state_scalar%rho(scalar_shift_index+1:scalar_shift_index+n_scalars)
        !$omp end parallel workshare
        
        call scalar_times_vector_h(diag%scalar_field_placeholder,diag%flux_density,diag%flux_density,grid)
        call div_h(diag%flux_density,diag%flux_density_div,grid%adjacent_signs_h,grid%adjacent_vector_indices_h, &
                   grid%inner_product_weights,grid%slope,grid%area,grid%volume)
        ! adding the tendencies in all grid boxes
        !$omp parallel do private(ji)
        do ji=1,n_scalars
          state_tend%rhotheta_v(ji) &
          = old_weight(jc)*state_tend%rhotheta_v(ji) &
          + new_weight(jc)*( &
          ! the advection (resolved transport)
          -diag%flux_density_div(ji) &
          ! the diabatic forcings
          ! weighting factor accounting for condensates
          + c_d_v*state_scalar%rho(scalar_shift_index+ji)/c_v_mass_weighted_air(state_scalar%rho,diag%temperature,ji-1)*( &
          ! dissipation through molecular + turbulent momentum diffusion
          diag%heating_diss(ji) &
          ! molecular + turbulent heat transport
          + diag%temperature_diffusion_heating(ji) &
          ! radiation
          + diag%radiation_tendency(ji) &
          ! phase transitions
          + diag%phase_trans_heating_rate(ji) &
          ! heating rate due to falling condensates
          + diag%condensates_sediment_heat(ji) &
          ! this has to be divided by c_p*exner
          )/(c_d_p*(grid%exner_bg(ji) + state_scalar%exner_pert(ji))) &
          ! tendency of due to phase transitions and mass diffusion
          + (diag%phase_trans_rates(scalar_shift_index+ji) + diag%mass_diff_tendency(scalar_shift_index+ji)) &
          *diag%scalar_field_placeholder(ji))
        enddo
        !$omp end parallel do
      endif
    enddo
  
  end subroutine scalar_tend_expl

end module mo_scalar_tend_expl




