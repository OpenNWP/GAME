! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_column_solvers

  ! This module contains the implicit vertical routines (implicit part of the HEVI scheme).

  use iso_c_binding
  use mo_definitions,      only: wp
  use mo_constants,        only: r_d,c_d_v,c_d_p,M_PI
  use mo_run_nml,          only: dtime
  use mo_grid_nml,         only: n_scalars,n_layers,n_scalars_h,n_vectors_per_layer,n_vectors
  use mo_constituents_nml, only: n_constituents,n_condensed_constituents,cloud_droplets_velocity,rain_velocity,&
                                 snow_velocity
  use mo_dictionary,       only: c_p_cond
  use mo_surface_nml,      only: nsoillays,lprog_soil_temp,lsfc_sensible_heat_flux
  use mo_diff_nml,         only: klemp_damp_max,klemp_begin_rel
  use mo_rad_nml,          only: radiation_dtime
  
  implicit none
  
  contains
  
  subroutine three_band_solver_ver_waves(z_vector,area,volume,scalar_flux_resistance,theta_v_pert_target, &
                                         rho_target,rhotheta_v_new,rho_new,rho_old,theta_v_pert_old,theta_v_bg, &
                                         exner_pert_old,exner_bg,theta_v_pert_new,exner_pert_new, &
                                         power_flux_density_sensible,temperature_soil_new,temperature_soil_old, &
                                         rhotheta_v_tend,is_land,rho_tend,rhotheta_v_old,wind_old, &
                                         z_scalar,wind_tend,z_soil_center,t_conduc_soil,sfc_rho_c, &
                                         sfc_lw_out,sfc_sw_in,t_const_soil,z_soil_interface, &
                                         power_flux_density_latent,rhotheta_v_target,exner_pert_target, &
                                         wind_target,temperature_soil_target,z_t_const,rk_step) &
  bind(c,name = "three_band_solver_ver_waves")
    
    ! This subroutine is the implicit vertical solver for the main fluid constituent.
    real(wp), intent(in)    :: z_vector(n_vectors),area(n_vectors),volume(n_scalars), &
                               scalar_flux_resistance(n_scalars_h),rho_old(n_constituents*n_scalars), &
                               rhotheta_v_new(n_scalars),rho_new(n_constituents*n_scalars),z_t_const, &
                               theta_v_pert_old(n_scalars),theta_v_bg(n_scalars),exner_pert_old(n_scalars), &
                               exner_bg(n_scalars),theta_v_pert_new(n_scalars),exner_pert_new(n_scalars), &
                               temperature_soil_new(nsoillays*n_scalars_h),temperature_soil_old(nsoillays*n_scalars_h), &
                               rho_tend(n_constituents*n_scalars),rhotheta_v_old(n_scalars),wind_old(n_vectors), &
                               z_scalar(n_scalars),wind_tend(n_vectors),z_soil_center(nsoillays), &
                               t_conduc_soil(n_scalars_h),sfc_rho_c(n_scalars_h),t_const_soil(n_scalars_h), &
                               sfc_lw_out(n_scalars_h),sfc_sw_in(n_scalars_h),z_soil_interface(nsoillays+1), &
                               power_flux_density_latent(n_scalars_h)
    real(wp), intent(out)   :: theta_v_pert_target(n_scalars),rho_target(n_constituents*n_scalars), &
                               power_flux_density_sensible(n_scalars_h), &
                               rhotheta_v_target(n_scalars),exner_pert_target(n_scalars),wind_target(n_vectors), &
                               temperature_soil_target(nsoillays*n_scalars_h)
    real(wp), intent(inout) :: rhotheta_v_tend(n_scalars)
    integer,  intent(in)    :: rk_step,is_land(n_scalars_h)
    
    ! local variables
    integer  :: i,j,lower_index,base_index,soil_switch,gas_phase_first_index
    real(wp) :: damping_coeff,damping_start_height,z_above_damping,temperature_gas_lowest_layer_old, &
                temperature_gas_lowest_layer_new,radiation_flux_density,resulting_temperature_change, &
                max_rad_temp_change,impl_weight,partial_deriv_new_time_step_weight, &
                c_vector(n_layers-2+nsoillays),d_vector(n_layers-1+nsoillays), &
                e_vector(n_layers-2+nsoillays),r_vector(n_layers-1+nsoillays), &
                rho_expl(n_layers),rhotheta_v_expl(n_layers),theta_v_pert_expl(n_layers),exner_pert_expl(n_layers), &
                theta_v_int_new(n_layers-1),solution_vector(n_layers-1+nsoillays),rho_int_old(n_layers-1), &
                rho_int_expl(n_layers-1),alpha_old(n_layers),beta_old(n_layers),gamma_old(n_layers), &
                alpha_new(n_layers),beta_new(n_layers),gamma_new(n_layers),alpha(n_layers),beta(n_layers), &
                gamma_(n_layers),density_interface_new,heat_flux_density_expl(nsoillays)
                    
    impl_weight = 0.75_wp
    
    ! the maximum temperature change induced by radiation between two radiation time steps in the uppermost soil layer
    max_rad_temp_change = 25._wp
    
    damping_start_height = klemp_begin_rel*z_vector(1)
    
    ! partial derivatives new time step weight
    partial_deriv_new_time_step_weight = 0.5_wp
    
    gas_phase_first_index = n_condensed_constituents*n_scalars
    
    ! calculating the sensible power flux density
    if (lsfc_sensible_heat_flux) then
      !$omp parallel do private(i,base_index,temperature_gas_lowest_layer_old,temperature_gas_lowest_layer_new)
      do i=1,n_scalars_h
        base_index = n_scalars - n_scalars_h + i
        
        ! gas temperature in the lowest layer
        temperature_gas_lowest_layer_old = (exner_bg(base_index) + exner_pert_old(base_index)) &
        *(theta_v_bg(base_index) + theta_v_pert_old(base_index))
        temperature_gas_lowest_layer_new = (exner_bg(base_index) + exner_pert_new(base_index)) &
        *(theta_v_bg(base_index) + theta_v_pert_new(base_index))
        
        ! the sensible power flux density
        power_flux_density_sensible(i) = 0.5_wp*c_d_v*(rho_new(gas_phase_first_index + base_index) &
        *(temperature_gas_lowest_layer_old - temperature_soil_old(i)) &
        + rho_old(gas_phase_first_index + base_index) &
        *(temperature_gas_lowest_layer_new - temperature_soil_new(i)))/scalar_flux_resistance(i)
        
        ! contribution of sensible heat to rhotheta_v
        rhotheta_v_tend(base_index) = rhotheta_v_tend(base_index) &
        - area(n_layers*n_vectors_per_layer + i)*power_flux_density_sensible(i) &
        /((exner_bg(base_index) + exner_pert_new(base_index))*c_d_p)/volume(base_index)
      enddo
      !$omp end parallel do
    endif
    
    ! loop over all columns
    !$omp parallel do private(i,j,base_index,lower_index,damping_coeff,z_above_damping,soil_switch, &
    !$omp radiation_flux_density,resulting_temperature_change,c_vector,d_vector, &
    !$omp e_vector,r_vector,rho_expl,rhotheta_v_expl,theta_v_pert_expl,exner_pert_expl,theta_v_int_new, &
    !$omp solution_vector,rho_int_old,rho_int_expl,alpha_old,beta_old,gamma_old,alpha_new,beta_new, &
    !$omp gamma_new,alpha,beta,gamma_,density_interface_new,heat_flux_density_expl)
    do i=1,n_scalars_h
      
      if (lprog_soil_temp) then
        soil_switch = is_land(i)
      else
        soil_switch = 0
      endif
      
      ! explicit quantities
      do j=1,n_layers
        base_index = i + (j-1)*n_scalars_h
        ! explicit density
        rho_expl(j) = rho_old(gas_phase_first_index + base_index) &
        + dtime*rho_tend(gas_phase_first_index + base_index)
        ! explicit virtual potential temperature density
        rhotheta_v_expl(j) = rhotheta_v_old(base_index) + dtime*rhotheta_v_tend(base_index)
        if (rk_step==0) then
          ! old time step partial derivatives of theta_v and Pi (divided by the volume)
          alpha(j) = -rhotheta_v_old(base_index)/rho_old(gas_phase_first_index + base_index)**2 &
          /volume(base_index)
          beta(j) = 1._wp/rho_old(gas_phase_first_index + base_index)/volume(base_index)
          gamma_(j) = r_d/(c_d_v*rhotheta_v_old(base_index)) &
          *(exner_bg(base_index) + exner_pert_old(base_index))/volume(base_index)
        else
          ! old time step partial derivatives of theta_v and Pi
          alpha_old(j) = -rhotheta_v_old(base_index)/rho_old(gas_phase_first_index + base_index)**2
          beta_old(j) = 1._wp/rho_old(gas_phase_first_index + base_index)
          gamma_old(j) = r_d/(c_d_v*rhotheta_v_old(base_index))*(exner_bg(base_index) + exner_pert_old(base_index))
          ! new time step partial derivatives of theta_v and Pi
          alpha_new(j) = -rhotheta_v_new(base_index)/rho_new(gas_phase_first_index + base_index)**2
          beta_new(j) = 1._wp/rho_new(gas_phase_first_index + base_index)
          gamma_new(j) = r_d/(c_d_v*rhotheta_v_new(base_index))*(exner_bg(base_index) + exner_pert_new(base_index))
          ! interpolation in time and dividing by the volume
          alpha(j) = ((1._wp - partial_deriv_new_time_step_weight)*alpha_old(j) &
          + partial_deriv_new_time_step_weight*alpha_new(j))/volume(base_index)
          beta(j) = ((1._wp - partial_deriv_new_time_step_weight)*beta_old(j) &
          + partial_deriv_new_time_step_weight*beta_new(j))/volume(base_index)
          gamma_(j) = ((1._wp - partial_deriv_new_time_step_weight)*gamma_old(j) &
          + partial_deriv_new_time_step_weight*gamma_new(j))/volume(base_index)
        endif
        ! explicit virtual potential temperature perturbation
        theta_v_pert_expl(j) = theta_v_pert_old(base_index) + dtime*volume(base_index)* &
        (alpha(j)*rho_tend(gas_phase_first_index + base_index) + beta(j)*rhotheta_v_tend(base_index))
        ! explicit Exner pressure perturbation
        exner_pert_expl(j) = exner_pert_old(base_index) + dtime*volume(base_index)*gamma_(j)*rhotheta_v_tend(base_index)
      enddo
      
      ! determining the interface values
      do j=1,n_layers-1
        base_index = i + (j-1)*n_scalars_h
        lower_index = i + j*n_scalars_h
        rho_int_old(j) = 0.5_wp*(rho_old(gas_phase_first_index + base_index) + rho_old(gas_phase_first_index + lower_index))
        rho_int_expl(j) = 0.5_wp*(rho_expl(j) + rho_expl(j+1))
        theta_v_int_new(j) = 0.5_wp*(rhotheta_v_new(base_index)/rho_new(gas_phase_first_index + base_index) &
        + rhotheta_v_new(lower_index)/rho_new(gas_phase_first_index + lower_index))
      enddo
      
      ! filling up the coefficient vectors
      do j=1,n_layers-1
        base_index = i + (j-1)*n_scalars_h
        lower_index = i + j*n_scalars_h
        ! main diagonal
        d_vector(j) = -theta_v_int_new(j)**2*(gamma_(j) + gamma_(j+1)) &
        + 0.5_wp*(exner_bg(base_index) - exner_bg(lower_index)) &
        *(alpha(j+1) - alpha(j) + theta_v_int_new(j)*(beta(j+1) - beta(j))) &
        - (z_scalar(base_index) - z_scalar(lower_index))/(impl_weight*dtime**2*c_d_p*rho_int_old(j)) &
        *(2._wp/area(i + j*n_vectors_per_layer) + dtime*wind_old(i + j*n_vectors_per_layer)*0.5_wp &
        *(-1._wp/volume(base_index) + 1._wp/volume(lower_index)))
        ! right hand side
        r_vector(j) = -(wind_old(i + j*n_vectors_per_layer) + dtime*wind_tend(i + j*n_vectors_per_layer)) &
        *(z_scalar(base_index) - z_scalar(lower_index)) &
        /(impl_weight*dtime**2*c_d_p) &
        + theta_v_int_new(j)*(exner_pert_expl(j) - exner_pert_expl(j+1))/dtime &
        + 0.5_wp/dtime*(theta_v_pert_expl(j) + theta_v_pert_expl(j+1))*(exner_bg(base_index) - exner_bg(lower_index)) &
        - (z_scalar(base_index) - z_scalar(lower_index))/(impl_weight*dtime**2*c_d_p) &
        *wind_old(i + j*n_vectors_per_layer)*rho_int_expl(j)/rho_int_old(j)
      enddo
      do j=1,n_layers-2
        base_index = i + (j-1)*n_scalars_h
        lower_index = i + j*n_scalars_h
        ! lower diagonal
        c_vector(j) = theta_v_int_new(j+1)*gamma_(j+1)*theta_v_int_new(j) &
        + 0.5_wp*(exner_bg(lower_index) - exner_bg((j+1)*n_scalars_h + i)) &
        *(alpha(j+1) + beta(j+1)*theta_v_int_new(j)) &
        - (z_scalar(lower_index) - z_scalar((j+1)*n_scalars_h + i))/(impl_weight*dtime*c_d_p)*0.5_wp &
        *wind_old(i + (j+1)*n_vectors_per_layer)/(volume(lower_index)*rho_int_old(j+1))
        ! upper diagonal
        e_vector(j) = theta_v_int_new(j)*gamma_(j+1)*theta_v_int_new(j+1) &
        - 0.5_wp*(exner_bg(base_index) - exner_bg(lower_index)) &
        *(alpha(j+1) + beta(j+1)*theta_v_int_new(j+1)) &
        + (z_scalar(base_index) - z_scalar(lower_index))/(impl_weight*dtime*c_d_p)*0.5_wp &
        *wind_old(i + j*n_vectors_per_layer)/(volume(lower_index)*rho_int_old(j))
      enddo
      
      ! soil components of the matrix
      if (soil_switch==1) then
        ! calculating the explicit part of the heat flux density
        do j=1,nsoillays-1
          heat_flux_density_expl(j) &
          = -sfc_rho_c(i)*t_conduc_soil(i)*(temperature_soil_old(i + (j-1)*n_scalars_h) &
          - temperature_soil_old(i + j*n_scalars_h))/(z_soil_center(j) - z_soil_center(j+1))
        enddo
        heat_flux_density_expl(nsoillays) &
        = -sfc_rho_c(i)*t_conduc_soil(i)*(temperature_soil_old(i + (nsoillays-1)*n_scalars_h) - t_const_soil(i)) &
        /(2._wp*(z_soil_center(nsoillays) - z_t_const))
        
        radiation_flux_density = sfc_sw_in(i) - sfc_lw_out(i)
        resulting_temperature_change = radiation_flux_density/((z_soil_interface(1) - z_soil_interface(2)) &
        *sfc_rho_c(i))*radiation_dtime
        if (abs(resulting_temperature_change)>max_rad_temp_change) then
          radiation_flux_density = max_rad_temp_change/abs(resulting_temperature_change)*radiation_flux_density
        endif
        
        ! calculating the explicit part of the temperature change
        r_vector(n_layers) &
        ! old temperature
        = temperature_soil_old(i) &
        ! sensible heat flux
        + (power_flux_density_sensible(i) &
        ! latent heat flux
        + power_flux_density_latent(i) &
        ! radiation
        + radiation_flux_density &
        ! heat conduction from below
        + 0.5_wp*heat_flux_density_expl(1)) &
        /((z_soil_interface(1) - z_soil_interface(2))*sfc_rho_c(i))*dtime
        
        ! loop over all soil layers below the first layer
        do j=2,nsoillays
          
          r_vector(j + n_layers - 1) &
          ! old temperature
          = temperature_soil_old(i + (j-1)*n_scalars_h) &
          ! heat conduction from above
          + 0.5_wp*(-heat_flux_density_expl(j-1) &
          ! heat conduction from below
          + heat_flux_density_expl(j)) &
          /((z_soil_interface(j) - z_soil_interface(j+1))*sfc_rho_c(i))*dtime
        enddo
        
        ! the diagonal component
        do j=1,nsoillays
          if (j==1) then
            d_vector(j + n_layers - 1) = 1._wp + 0.5_wp*dtime*sfc_rho_c(i)*t_conduc_soil(i) &
            /((z_soil_interface(j) - z_soil_interface(j+1))*sfc_rho_c(i)) &
            *1._wp/(z_soil_center(j) - z_soil_center(j+1))
          elseif (j==nsoillays) then
            d_vector(j + n_layers - 1) = 1._wp + 0.5_wp*dtime*sfc_rho_c(i)*t_conduc_soil(i) &
            /((z_soil_interface(j) - z_soil_interface(j+1))*sfc_rho_c(i)) &
            *1._wp/(z_soil_center(j-1) - z_soil_center(j))
          else
            d_vector(j + n_layers - 1) = 1._wp + 0.5_wp*dtime*sfc_rho_c(i)*t_conduc_soil(i) &
            /((z_soil_interface(j) - z_soil_interface(j+1))*sfc_rho_c(i)) &
            *(1._wp/(z_soil_center(j-1) - z_soil_center(j)) &
            + 1._wp/(z_soil_center(j) - z_soil_center(j+1)))
          endif
        enddo
        ! the off-diagonal components
        c_vector(n_layers-1) = 0._wp
        e_vector(n_layers-1) = 0._wp
        do j=1,nsoillays-1
          c_vector(j + n_layers - 1) = -0.5_wp*dtime*sfc_rho_c(i)*t_conduc_soil(i) &
          /((z_soil_interface(j+1) - z_soil_interface(j+2))*sfc_rho_c(i)) &
          /(z_soil_center(j) - z_soil_center(j+1))
          e_vector(j + n_layers - 1) = -0.5_wp*dtime*sfc_rho_c(i)*t_conduc_soil(i) &
          /((z_soil_interface(j) - z_soil_interface(j+1))*sfc_rho_c(i)) &
          /(z_soil_center(j) - z_soil_center(j+1))
        enddo
      endif
      
      ! calling the algorithm to solve the system of linear equations
      call thomas_algorithm(c_vector,d_vector,e_vector,r_vector,solution_vector,n_layers-1+soil_switch*nsoillays)
      
      ! Klemp (2008) upper boundary layer
      do j=1,n_layers-1
        base_index = i + (j-1)*n_scalars_h
        z_above_damping = z_vector(i + j*n_vectors_per_layer) - damping_start_height
        if (z_above_damping<0._wp) then
          damping_coeff = 0._wp
        else
          damping_coeff = klemp_damp_max*(sin(0.5_wp*M_PI*z_above_damping/(z_vector(1) - damping_start_height)))**2
        endif
        solution_vector(j) = solution_vector(j)/(1._wp+dtime*damping_coeff)
      enddo
      
      ! Writing the result into the new state.
      ! --------------------------------------
      
      ! mass density
      do j=1,n_layers
        base_index = i + (j-1)*n_scalars_h
        if (j==1) then
          rho_target(gas_phase_first_index + base_index) &
          = rho_expl(j) + dtime*(solution_vector(j))/volume(base_index)
        elseif (j==n_layers) then
          rho_target(gas_phase_first_index + base_index) &
          = rho_expl(j) + dtime*(-solution_vector(j - 1))/volume(base_index)
        else
          rho_target(gas_phase_first_index + base_index) &
          = rho_expl(j) + dtime*(-solution_vector(j - 1) + solution_vector(j))/volume(base_index)
        endif
      enddo
      ! virtual potential temperature density
      do j=1,n_layers
        base_index = i + (j-1)*n_scalars_h
        if (j==1) then
          rhotheta_v_target(base_index) &
          = rhotheta_v_expl(j) + dtime*(theta_v_int_new(j)*solution_vector(j))/volume(base_index)
        elseif (j==n_layers) then
          rhotheta_v_target(base_index) &
          = rhotheta_v_expl(j) + dtime*(-theta_v_int_new(j-1)*solution_vector(j-1))/volume(base_index)
        else
          rhotheta_v_target(base_index) &
          = rhotheta_v_expl(j) + dtime*(-theta_v_int_new(j-1)*solution_vector(j-1) + theta_v_int_new(j)*solution_vector(j)) &
          /volume(base_index)
        endif
      enddo
      ! vertical velocity
      do j=1,n_layers-1
        base_index = i + (j-1)*n_scalars_h
        density_interface_new = 0.5_wp*(rho_target(gas_phase_first_index + base_index) &
        + rho_target(gas_phase_first_index + i + j*n_scalars_h))
        wind_target(i + j*n_vectors_per_layer) &
        = (2._wp*solution_vector(j)/area(i + j*n_vectors_per_layer) &
        - density_interface_new*wind_old(i + j*n_vectors_per_layer))/rho_int_old(j)
      enddo
      ! virtual potential temperature perturbation
      do j=1,n_layers
        base_index = i + (j-1)*n_scalars_h
        theta_v_pert_target(base_index) = rhotheta_v_target(base_index)/rho_target(gas_phase_first_index + base_index) &
        - theta_v_bg(base_index)
      enddo
      ! Exner pressure perturbation
      do j=1,n_layers
        base_index = i + (j-1)*n_scalars_h
        exner_pert_target(base_index) = exner_pert_old(base_index) + volume(base_index) &
        *gamma_(j)*(rhotheta_v_target(base_index) - rhotheta_v_old(base_index))
      enddo
      
      ! soil temperature
      if (soil_switch==1) then
        do j=1,nsoillays
          temperature_soil_target(i + (j-1)*n_scalars_h) = solution_vector(n_layers - 1 + j)
        enddo
      endif
      
    enddo ! end of the column loop
    !$omp end parallel do
    
  end subroutine three_band_solver_ver_waves

  subroutine three_band_solver_gen_densities(wind_old,wind_new,volume,rho_tend,rho_old,rho_new,&
                                             condensates_sediment_heat,area,temperature,rk_step) &
  bind(c,name = "three_band_solver_gen_densities")
  
    ! This subroutine contains the vertical advection of mass densities (of tracers) with 3-band matrices.
    
    integer,  intent(in)  :: rk_step
    real(wp), intent(in)  :: wind_old(n_vectors),wind_new(n_vectors),volume(n_scalars),&
                             rho_tend(n_constituents*n_scalars),rho_old(n_constituents*n_scalars),&
                             area(n_vectors),temperature(n_scalars)
    real(wp), intent(out) :: condensates_sediment_heat(n_scalars),rho_new(n_constituents*n_scalars)
    
    ! local variables
    integer  :: ji,jl,jc,lower_index,upper_index,base_index
    real(wp) :: impl_weight,expl_weight,density_old_at_interface,temperature_old_at_interface,&
                ! for meanings of these vectors look into the definition of the function thomas_algorithm
                c_vector(n_layers-1),d_vector(n_layers),e_vector(n_layers-1),r_vector(n_layers),&
                vertical_flux_vector_impl(n_layers-1),vertical_flux_vector_rhs(n_layers-1),&
                vertical_enthalpy_flux_vector(n_layers-1),solution_vector(n_layers)
          
    impl_weight = 0.5_wp
    expl_weight = 1._wp - impl_weight
    
    ! loop over all relevant constituents
    do jc=0,n_constituents-1
    
      ! This is done for all tracers apart from the main gaseous constituent.
      if (jc/=n_condensed_constituents) then
        ! loop over all columns
        !$omp parallel do private(ji,jl,lower_index,upper_index,base_index,density_old_at_interface,&
        !$omp temperature_old_at_interface,c_vector,d_vector,e_vector,r_vector,vertical_flux_vector_impl,&
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
            condensates_sediment_heat((n_layers-1)*n_scalars_h + ji) = 0._wp
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
    
    integer, intent(in)  :: solution_length                  ! length of the solution vector
    real(wp),intent(in)  :: c_vector(solution_length)        ! lower diagonal vector
    real(wp),intent(in)  :: d_vector(solution_length)        ! main diagonal vector
    real(wp),intent(in)  :: e_vector(solution_length)        ! upper diagonal vector
    real(wp),intent(in)  :: r_vector(solution_length)        ! right hand side vector
    real(wp),intent(out) :: solution_vector(solution_length) ! vector containing the solution
    
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

end module mo_column_solvers











