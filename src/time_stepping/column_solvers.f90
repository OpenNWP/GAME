! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module column_solvers

  ! This module contains the implicit vertical routines (implicit part of the HEVI scheme).

  use iso_c_binding
  use definitions,      only: wp
  use run_nml,          only: dtime
  use grid_nml,         only: n_scalars,n_layers,n_scalars_h,n_vectors_per_layer,n_vectors
  use constituents_nml, only: n_constituents,n_condensed_constituents,cloud_droplets_velocity,rain_velocity, &
                              snow_velocity
  use dictionary,       only: c_p_cond
  
  implicit none
  
  contains

  subroutine three_band_solver_gen_densities(wind_old,wind_new,volume,rho_tend,rho_old,rho_new, &
                                             condensates_sediment_heat,area,temperature,rk_step) &
  bind(c,name = "three_band_solver_gen_densities")
  
    ! This subroutine contains the vertical advection of mass densities (of tracers) with 3-band matrices.
    
    integer,  intent(in)  :: rk_step
    real(wp), intent(in)  :: wind_old(n_vectors),wind_new(n_vectors),volume(n_scalars), &
                             rho_tend(n_constituents*n_scalars),rho_old(n_constituents*n_scalars), &
                             area(n_vectors),temperature(n_scalars)
    real(wp), intent(out) :: condensates_sediment_heat(n_scalars),rho_new(n_constituents*n_scalars)
    
    ! local variables
    integer  :: ji,jl,jc,lower_index,upper_index,base_index
    real(wp) :: impl_weight,expl_weight,density_old_at_interface,temperature_old_at_interface, &
                ! for meanings of these vectors look into the definition of the function thomas_algorithm
                c_vector(n_layers-1),d_vector(n_layers),e_vector(n_layers-1),r_vector(n_layers), &
                vertical_flux_vector_impl(n_layers-1),vertical_flux_vector_rhs(n_layers-1), &
                vertical_enthalpy_flux_vector(n_layers-1),solution_vector(n_layers)
          
    impl_weight = 0.5_wp
    expl_weight = 1._wp - impl_weight
    
    ! loop over all relevant constituents
    do jc=0,n_constituents-1
    
      ! This is done for all tracers apart from the main gaseous constituent.
      if (jc/=n_condensed_constituents) then
        ! loop over all columns
        !$omp parallel do private(ji,jl,lower_index,upper_index,base_index,density_old_at_interface, &
        !$omp temperature_old_at_interface,c_vector,d_vector,e_vector,r_vector,vertical_flux_vector_impl, &
        !$omp vertical_flux_vector_rhs,vertical_enthalpy_flux_vector,solution_vector)
        do ji=1,n_scalars_h
          
          ! diagnozing the vertical fluxes
          do jl=1,n_layers-1
            ! resetting the vertical enthalpy flux density divergence
            if (rk_step==0 .and. jc==0) then
              condensates_sediment_heat((jl-1)*n_scalars_h + ji) = 0._wp
            endif
            base_index = ji + (jl-1)*n_scalars_h
            vertical_flux_vector_impl(jl) = wind_old(ji + jl*n_vectors_per_layer)
            vertical_flux_vector_rhs(jl) = wind_new(ji + jl*n_vectors_per_layer)
            ! preparing the vertical interpolation
            lower_index = ji + jl*n_scalars_h
            upper_index = base_index
            ! For condensed constituents, a sink velocity must be added.
            ! precipitation
            ! snow
            if (jc<n_condensed_constituents/4) then
              vertical_flux_vector_impl(jl) = vertical_flux_vector_impl(jl) - snow_velocity
              vertical_flux_vector_rhs(jl) = vertical_flux_vector_rhs(jl) - snow_velocity
            ! rain
            elseif (jc<n_condensed_constituents/2) then
              vertical_flux_vector_impl(jl) = vertical_flux_vector_impl(jl) - rain_velocity
              vertical_flux_vector_rhs(jl) = vertical_flux_vector_rhs(jl) - rain_velocity
            ! clouds
            elseif (jc<n_condensed_constituents) then
              vertical_flux_vector_impl(jl) = vertical_flux_vector_impl(jl) - cloud_droplets_velocity
              vertical_flux_vector_rhs(jl) = vertical_flux_vector_rhs(jl) - cloud_droplets_velocity
            endif
            ! multiplying the vertical velocity by the
            vertical_flux_vector_impl(jl) = area(ji + jl*n_vectors_per_layer)*vertical_flux_vector_impl(jl)
            vertical_flux_vector_rhs(jl) = area(ji + jl*n_vectors_per_layer)*vertical_flux_vector_rhs(jl)
            ! old density at the interface
            if (vertical_flux_vector_rhs(jl)>=0._wp) then
              density_old_at_interface = rho_old(jc*n_scalars + lower_index)
              temperature_old_at_interface = temperature(lower_index)
            else
              density_old_at_interface = rho_old(jc*n_scalars + upper_index)
              temperature_old_at_interface = temperature(upper_index)
            endif
            vertical_flux_vector_rhs(jl) = density_old_at_interface*vertical_flux_vector_rhs(jl)
            vertical_enthalpy_flux_vector(jl) = c_p_cond(jc,temperature_old_at_interface) &
            *temperature_old_at_interface*vertical_flux_vector_rhs(jl)
          enddo
          if (rk_step==0 .and. jc==0) then
            condensates_sediment_heat((n_layers - 1)*n_scalars_h + ji) = 0._wp
          endif
          
          ! Now we proceed to solving the vertical tridiagonal problems.
          
          ! filling up the original vectors
          do jl=1,n_layers-1
            base_index = ji + (jl-1)*n_scalars_h
            if (vertical_flux_vector_impl(jl)>=0._wp) then
              c_vector(jl) = 0._wp
              e_vector(jl) = -impl_weight*dtime/volume(base_index)*vertical_flux_vector_impl(jl)
            else
              c_vector(jl) = impl_weight*dtime/volume(ji + jl*n_scalars_h)*vertical_flux_vector_impl(jl)
              e_vector(jl) = 0._wp
            endif
          enddo
          do jl=1,n_layers
            base_index = ji + (jl-1)*n_scalars_h
            if (jl==1) then
              if (vertical_flux_vector_impl(1)>=0._wp) then
                d_vector(jl) = 1._wp
              else
                d_vector(jl) = 1._wp - impl_weight*dtime/volume(base_index)*vertical_flux_vector_impl(1)
              endif
            elseif (jl==n_layers) then
              if (vertical_flux_vector_impl(jl-1)>=0._wp) then
                d_vector(jl) = 1._wp + impl_weight*dtime/volume(base_index)*vertical_flux_vector_impl(jl-1)
              else
                d_vector(jl) = 1._wp
              endif
              ! precipitation
              ! snow
              if (jc<n_condensed_constituents/4) then
                d_vector(jl) = d_vector(jl) + impl_weight*snow_velocity*dtime &
                *area(ji + n_vectors - n_scalars_h)/volume(base_index)
              ! rain
              elseif (jc<n_condensed_constituents/2) then
                d_vector(jl) = d_vector(jl) + impl_weight*rain_velocity*dtime &
                *area(ji + n_vectors - n_scalars_h)/volume(base_index)
              ! clouds
              elseif (jc<n_condensed_constituents) then
                d_vector(jl) = d_vector(jl) + impl_weight*cloud_droplets_velocity*dtime &
                *area(ji + n_vectors - n_scalars_h)/volume(base_index)
              endif
            else
              d_vector(jl) = 1._wp
              if (vertical_flux_vector_impl(jl-1)>=0._wp) then
                d_vector(jl) = d_vector(jl) + impl_weight*dtime/volume(base_index)*vertical_flux_vector_impl(jl-1)
              endif
              if (vertical_flux_vector_impl(jl)<0._wp) then
                d_vector(jl) = d_vector(jl) - impl_weight*dtime/volume(base_index)*vertical_flux_vector_impl(jl)
              endif
            endif
            ! the explicit component
            ! mass densities
            r_vector(jl) = rho_old(jc*n_scalars + base_index) + dtime*rho_tend(jc*n_scalars + base_index)
            ! adding the explicit part of the vertical flux divergence
            if (jl==1) then
              r_vector(jl) = r_vector(jl) + expl_weight*dtime*vertical_flux_vector_rhs(jl)/volume(base_index)
              if (rk_step==0 .and. jc<n_condensed_constituents) then
                condensates_sediment_heat(base_index) = condensates_sediment_heat(base_index) &
                + vertical_enthalpy_flux_vector(jl)/volume(base_index)
              endif
            elseif (jl==n_layers) then
              r_vector(jl) = r_vector(jl) - expl_weight*dtime*vertical_flux_vector_rhs(jl-1)/volume(base_index)
              if (rk_step==0 .and. jc<n_condensed_constituents) then
                condensates_sediment_heat(base_index) = condensates_sediment_heat(base_index) &
                - vertical_enthalpy_flux_vector(jl-1)/volume(base_index)
              endif
              ! precipitation
              ! snow
              if (jc<n_condensed_constituents/4) then
                r_vector(jl) = r_vector(jl) - expl_weight*snow_velocity*dtime &
                *rho_old(jc*n_scalars + ji + n_scalars - n_scalars_h) &
                *area(ji + n_vectors - n_scalars_h)/volume(base_index)
                if (rk_step==0) then
                  condensates_sediment_heat(base_index) = condensates_sediment_heat(base_index) &
                  - snow_velocity &
                  *temperature(ji + n_scalars - n_scalars_h)*c_p_cond(jc,temperature(ji + n_scalars - n_scalars_h)) &
                  *rho_old(jc*n_scalars + ji + n_scalars - n_scalars_h) &
                  *area(ji + n_vectors - n_scalars_h)/volume(base_index)
                endif
              ! rain
              elseif (jc<n_condensed_constituents/2) then
                r_vector(jl) = r_vector(jl) - expl_weight*rain_velocity*dtime &
                *rho_old(jc*n_scalars + ji + n_scalars - n_scalars_h) &
                *area(ji + n_vectors - n_scalars_h)/volume(base_index)
                if (rk_step==0) then
                  condensates_sediment_heat(base_index) = condensates_sediment_heat(base_index) &
                  -rain_velocity &
                  *temperature(ji + n_scalars - n_scalars_h)*c_p_cond(jc,temperature(ji + n_scalars - n_scalars_h)) &
                  *rho_old(jc*n_scalars + ji + n_scalars - n_scalars_h) &
                  *area(ji + n_vectors - n_scalars_h)/volume(base_index)
                endif
              ! clouds
              elseif (jc<n_condensed_constituents) then
                r_vector(jl) = r_vector(jl) - expl_weight*cloud_droplets_velocity*dtime &
                *rho_old(jc*n_scalars + ji + n_scalars - n_scalars_h) &
                *area(ji + n_vectors - n_scalars_h)/volume(base_index)
                if (rk_step==0) then
                  condensates_sediment_heat(base_index) = condensates_sediment_heat(base_index) &
                  -cloud_droplets_velocity &
                  *temperature(ji + n_scalars - n_scalars_h)*c_p_cond(jc,temperature(ji + n_scalars - n_scalars_h)) &
                  *rho_old(jc*n_scalars + ji + n_scalars - n_scalars_h) &
                  *area(ji + n_vectors - n_scalars_h)/volume(base_index)
                endif
              endif
            else
              r_vector(jl) = r_vector(jl) + expl_weight*dtime &
              *(-vertical_flux_vector_rhs(jl-1) + vertical_flux_vector_rhs(jl))/volume(base_index)
              if (rk_step==0 .and. jc<n_condensed_constituents) then
                condensates_sediment_heat(base_index) = condensates_sediment_heat(base_index) &
                + (-vertical_enthalpy_flux_vector(jl-1) + vertical_enthalpy_flux_vector(jl))/volume(base_index)
              endif
            endif
          enddo
          
          ! calling the algorithm to solve the system of linear equations
          call thomas_algorithm(c_vector,d_vector,e_vector,r_vector,solution_vector,n_layers)
          
          ! this should account for round-off errors only
          do jl=1,n_layers
            if (solution_vector(jl)<0._wp) then
              solution_vector(jl) = 0._wp
            endif
          enddo
          
          ! writing the result into the new state
          do jl=1,n_layers
            base_index = ji + (jl-1)*n_scalars_h
            rho_new(jc*n_scalars + base_index) = solution_vector(jl)
          enddo
        enddo ! horizontal index
        !$omp end parallel do
      endif
    enddo ! constituent
    
  end subroutine three_band_solver_gen_densities

  subroutine thomas_algorithm(c_vector,d_vector,e_vector,r_vector,solution_vector,solution_length) &
  bind(c,name = "thomas_algorithm")

    ! This subroutine solves a system of linear equations with a three-band matrix.
    
    integer,  intent(in)  :: solution_length                  ! length of the solution vector
    real(wp), intent(in)  :: c_vector(solution_length)        ! lower diagonal vector
    real(wp), intent(in)  :: d_vector(solution_length)        ! main diagonal vector
    real(wp), intent(in)  :: e_vector(solution_length)        ! upper diagonal vector
    real(wp), intent(in)  :: r_vector(solution_length)        ! right hand side vector
    real(wp), intent(out) :: solution_vector(solution_length) ! vector containing the solution
    
    ! local variables
    real(wp) :: e_prime_vector(solution_length-1) ! help vector for solving the matrix equation
    real(wp) :: r_prime_vector(solution_length)   ! help vector for solving the matrix equation
    integer  :: jl                                ! loop index
    
    ! downward sweep (matrix)
    e_prime_vector(1) = e_vector(1)/d_vector(1)
    do jl=2,solution_length-1
      e_prime_vector(jl) = e_vector(jl)/(d_vector(jl) - e_prime_vector(jl-1)*c_vector(jl-1))
    enddo
    ! downward sweep (right-hand side)
    r_prime_vector(1) = r_vector(1)/d_vector(1)
    do jl=2,solution_length
      r_prime_vector(jl) = (r_vector(jl) - r_prime_vector(jl-1)*c_vector(jl-1)) &
      /(d_vector(jl) - e_prime_vector(jl-1)*c_vector(jl-1))
    enddo
    
    ! upward sweep (final solution)
    solution_vector(solution_length) = r_prime_vector(solution_length)
    do jl=solution_length-1,1,-1
      solution_vector(jl) = r_prime_vector(jl) - e_prime_vector(jl)*solution_vector(jl+1)
    enddo
  
  end subroutine thomas_algorithm

end module column_solvers











