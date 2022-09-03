! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME
  
module mo_scalar_tend_expl
  
  ! This module contains the horizontal (explicit) part of the constituent integration.

  use mo_definitions,           only: wp
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

  subroutine scalar_tend_expl(rho,mass_diff_tendency,scalar_field_placeholder,adjacent_vector_indices_h, &
                              rhotheta_v_tend,adjacent_signs_h,area,flux_density,from_index,to_index, &
                              inner_product_weights,layer_thickness,mass_diffusion_coeff_numerical_h, &
                              mass_diffusion_coeff_numerical_v,molecular_diffusion_coeff,n_squared, &
                              normal_distance,rhotheta_v,slope,temp_diffusion_coeff_numerical_h, &
                              temp_diffusion_coeff_numerical_v,temperature,tke,vector_field_placeholder, &
                              viscosity,viscosity_rhombi,viscosity_triangles,volume,vorticity_indices_triangles, &
                              wind,temperature_diffusion_heating,flux_density_div,rho_tend,phase_trans_rates, &
                              exner_pert,exner_bg,condensates_sediment_heat,phase_trans_heating_rate, &
                              radiation_tendency,heating_diss,rk_step)
    
    ! This subroutine manages the calculation of the explicit part of the scalar tendencies.
    
    real(wp), intent(in)    :: rho(n_constituents*n_scalars),area(n_vectors),inner_product_weights(8*n_scalars), &
                               layer_thickness(n_scalars),n_squared(n_scalars),normal_distance(n_vectors), &
                               rhotheta_v(n_scalars),slope(n_vectors),temperature(n_scalars),tke(n_scalars), &
                               volume(n_scalars),wind(n_vectors),phase_trans_rates((n_condensed_constituents+1)*n_scalars), &
                               exner_pert(n_scalars),exner_bg(n_scalars),condensates_sediment_heat(n_scalars), &
                               phase_trans_heating_rate(n_scalars),radiation_tendency(n_scalars),heating_diss(n_scalars)
    real(wp), intent(out)   :: mass_diff_tendency(n_constituents*n_scalars),rhotheta_v_tend(n_scalars), &
                               scalar_field_placeholder(n_scalars), &
                               flux_density(n_vectors),mass_diffusion_coeff_numerical_h(n_scalars), &
                               mass_diffusion_coeff_numerical_v(n_scalars),molecular_diffusion_coeff(n_scalars), &
                               temp_diffusion_coeff_numerical_h(n_scalars),temp_diffusion_coeff_numerical_v(n_scalars), &
                               vector_field_placeholder(n_vectors),temperature_diffusion_heating(n_scalars), &
                               flux_density_div(n_scalars),rho_tend(n_constituents*n_scalars)
    integer,  intent(in)    :: adjacent_signs_h(6*n_scalars_h),adjacent_vector_indices_h(6*n_scalars_h), &
                               from_index(n_vectors_h),to_index(n_vectors_h),rk_step, &
                               vorticity_indices_triangles(3*n_dual_scalars_h)
    real(wp), intent(inout) :: viscosity(n_scalars),viscosity_rhombi(n_vectors),viscosity_triangles(n_dual_v_vectors)
    
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
      call scalar_diffusion_coeffs(temperature,tke,rho,from_index,to_index,vorticity_indices_triangles, &
                                   molecular_diffusion_coeff,viscosity_triangles,viscosity,viscosity_rhombi, &
                                   mass_diffusion_coeff_numerical_h,mass_diffusion_coeff_numerical_v, &
                                   temp_diffusion_coeff_numerical_h,temp_diffusion_coeff_numerical_v, &
                                   n_squared,layer_thickness)
    endif
    
    ! Temperature diffusion gets updated at the first RK step if required.
    if (ltemp_diff_h .and. rk_step==1) then
      ! The diffusion of the temperature depends on its gradient.
      call grad(temperature,vector_field_placeholder,from_index,to_index,normal_distance,inner_product_weights,slope)
      ! Now the diffusive temperature flux density can be obtained.
      call scalar_times_vector_h(temp_diffusion_coeff_numerical_h,vector_field_placeholder,flux_density,from_index,to_index)
      ! The divergence of the diffusive temperature flux density is the diffusive temperature heating.
      call div_h(flux_density,temperature_diffusion_heating, &
                 adjacent_signs_h,adjacent_vector_indices_h,inner_product_weights,slope,area,volume)
      ! vertical temperature diffusion
      if (ltemp_diff_v) then
        call scalar_times_vector_v(temp_diffusion_coeff_numerical_v,vector_field_placeholder,flux_density)
        call add_vertical_div(flux_density,temperature_diffusion_heating,area,volume)
      endif
    endif
    
    ! Mass diffusion gets updated at the first RK step if required.
    if (lmass_diff_h .and. rk_step==1) then
      ! loop over all constituents
      do jc=1,n_constituents
        scalar_shift_index = (jc-1)*n_scalars

        ! The diffusion of the tracer density depends on its gradient.
        call grad(rho(scalar_shift_index+1:scalar_shift_index+n_scalars),vector_field_placeholder, &
                  from_index,to_index,normal_distance,inner_product_weights,slope)
        ! Now the diffusive mass flux density can be obtained.
        call scalar_times_vector_h(mass_diffusion_coeff_numerical_h, &
                                   vector_field_placeholder,vector_field_placeholder,from_index,to_index)
        ! The divergence of the diffusive mass flux density is the diffusive mass source rate.
        call div_h(vector_field_placeholder,mass_diff_tendency(scalar_shift_index+1:scalar_shift_index+n_scalars), &
                   adjacent_signs_h,adjacent_vector_indices_h,inner_product_weights,slope,area,volume)
        ! vertical mass diffusion
        if (lmass_diff_v) then
          call scalar_times_vector_v(mass_diffusion_coeff_numerical_v,vector_field_placeholder,vector_field_placeholder)
          call add_vertical_div(vector_field_placeholder,mass_diff_tendency(scalar_shift_index+1:scalar_shift_index+n_scalars), &
                                area,volume)
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
        call scalar_times_vector_h(rho(scalar_shift_index+1:scalar_shift_index+n_scalars),wind,flux_density,from_index,to_index)
        call div_h(flux_density,flux_density_div, &
                   adjacent_signs_h,adjacent_vector_indices_h,inner_product_weights,slope,area,volume)
      ! all other constituents
      else
        call scalar_times_vector_h_upstream(rho(scalar_shift_index+1:scalar_shift_index+n_scalars), &
                                            wind,flux_density,from_index,to_index)
        call div_h_tracer(flux_density,rho(scalar_shift_index+1:scalar_shift_index+n_scalars),wind,flux_density_div, &
                          adjacent_signs_h,adjacent_vector_indices_h,inner_product_weights,slope,area,volume)
      endif
      
      ! adding the tendencies in all grid boxes
      !$omp parallel do private(ji,scalar_index)
      do ji=1,n_scalars
        scalar_index = scalar_shift_index + ji
        rho_tend(scalar_index) &
        = old_weight(jc)*rho_tend(scalar_index) &
        + new_weight(jc)*( &
        ! the advection
        -flux_density_div(ji) &
        ! the diffusion
        + mass_diff_tendency(scalar_shift_index+ji) &
        ! phase transitions
        + phase_trans_rates(scalar_shift_index_phase_trans+ji))
      enddo
      !$omp end parallel do
      
      ! Explicit component of the rho*theta_v integration
      ! -------------------------------------------------
      
      if (jc==n_condensed_constituents+1) then
        ! determining the virtual potential temperature
        !$omp parallel workshare
        scalar_field_placeholder = rhotheta_v/rho(scalar_shift_index+1:scalar_shift_index+n_scalars)
        !$omp end parallel workshare
        
        call scalar_times_vector_h(scalar_field_placeholder,flux_density,flux_density,from_index,to_index)
        call div_h(flux_density,flux_density_div,adjacent_signs_h,adjacent_vector_indices_h, &
                   inner_product_weights,slope,area,volume)
        ! adding the tendencies in all grid boxes
        !$omp parallel do private(ji)
        do ji=1,n_scalars
          rhotheta_v_tend(ji) &
          = old_weight(jc)*rhotheta_v_tend(ji) &
          + new_weight(jc)*( &
          ! the advection (resolved transport)
          -flux_density_div(ji) &
          ! the diabatic forcings
          ! weighting factor accounting for condensates
          + c_d_v*rho(scalar_shift_index+ji)/c_v_mass_weighted_air(rho,temperature,ji-1)*( &
          ! dissipation through molecular + turbulent momentum diffusion
          heating_diss(ji) &
          ! molecular + turbulent heat transport
          + temperature_diffusion_heating(ji) &
          ! radiation
          + radiation_tendency(ji) &
          ! phase transitions
          + phase_trans_heating_rate(ji) &
          ! heating rate due to falling condensates
          + condensates_sediment_heat(ji) &
          ! this has to be divided by c_p*exner
          )/(c_d_p*(exner_bg(ji) + exner_pert(ji))) &
          ! tendency of due to phase transitions and mass diffusion
          + (phase_trans_rates(scalar_shift_index+ji) + mass_diff_tendency(scalar_shift_index+ji)) &
          *scalar_field_placeholder(ji))
        enddo
        !$omp end parallel do
      endif
    enddo
  
  end subroutine scalar_tend_expl

end module mo_scalar_tend_expl




