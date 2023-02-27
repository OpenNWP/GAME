! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_column_solvers
  
  ! This module contains the implicit vertical routines (implicit part of the HEVI scheme).
  
  use mo_definitions,      only: wp,t_grid,t_state,t_diag
  use mo_constants,        only: r_d,c_d_v,c_d_p,M_PI,m_d,m_v,impl_thermo_weight
  use mo_grid_nml,         only: n_layers,n_cells,n_levels
  use mo_grid_setup,       only: z_t_const
  use mo_constituents_nml, only: n_constituents,n_condensed_constituents,lmoist
  use mo_grid_setup,       only: dtime,toa
  use mo_dictionary,       only: c_p_cond,snow_particles_radius,ice_particles_radius,cloud_droplets_radius
  use mo_derived,          only: v_fall_solid,v_fall_liquid
  use mo_surface_nml,      only: nsoillays,lprog_soil_temp,lsfc_sensible_heat_flux
  use mo_diff_nml,         only: klemp_damp_max,klemp_begin_rel
  use mo_rad_nml,          only: radiation_dtime
  
  implicit none
  
  contains
  
  subroutine three_band_solver_ver_waves(state_old,state_new,state_tend,diag,grid,rk_step)
    
    ! This subroutine is the implicit vertical solver for the main fluid constituent.
    
    type(t_state), intent(in),    target :: state_old  ! state variables at the old time step
    type(t_state), intent(inout), target :: state_new  ! state variables at the new time step
    type(t_state), intent(inout)         :: state_tend ! state containing the partial temporal derivatives of the state variables
    type(t_diag),  intent(inout)         :: diag       ! diagnostic quantities
    type(t_grid),  intent(in)            :: grid       ! grid quantities
    integer,       intent(in)            :: rk_step    ! predictor-corrector step index
    
    ! local variables
    integer                :: ji                                    ! cell index
    integer                :: jl                                    ! layer or level index
    integer                :: soil_switch                           ! soil switch: 0 if soil does not have to be calculated, 1 if soil has to be calculated
    real(wp)               :: damping_coeff                         ! damping coefficient of the Klemp layer
    real(wp)               :: damping_start_height                  ! lower boundary height of the Klemp layer
    real(wp)               :: damping_prefactor(n_layers-1)         ! damping coefficient of the Klemp layer
    real(wp)               :: z_above_damping                       ! height of a level above damping_start_height
    real(wp)               :: t_gas_lowest_layer_old                ! temperature of the gas in the lowest layer of the model atmosphere, old time step
    real(wp)               :: t_gas_lowest_layer_new                ! temperature of the gas in the lowest layer of the model atmosphere, new time step
    real(wp)               :: radiation_flux_density                ! radiation flux density into the surface, needed for stabilizing the model at large time steps
    real(wp)               :: resulting_temperature_change          ! temperature change resulting from radiation_flux_density
    real(wp)               :: max_rad_temp_change                   ! the maximum radiation-induced temperature change of the uppermost soil layer
    real(wp)               :: partial_deriv_new_time_step_weight    ! partial derivatives weight of the new time step
    real(wp)               :: c_vector(n_layers-2+nsoillays)        ! needed for the vertical solver
    real(wp)               :: d_vector(n_layers-1+nsoillays)        ! needed for the vertical solver
    real(wp)               :: e_vector(n_layers-2+nsoillays)        ! needed for the vertical solver
    real(wp)               :: r_vector(n_layers-1+nsoillays)        ! needed for the vertical solver
    real(wp)               :: rho_expl(n_layers)                    ! explicit mass density
    real(wp)               :: rhotheta_v_expl(n_layers)             ! explicit virtual potential temperature density
    real(wp)               :: theta_v_pert_expl(n_layers)           ! explicit virtual potential temperature perturbation
    real(wp)               :: exner_pert_expl(n_layers)             ! explicit Exner pressure perturbation
    real(wp)               :: theta_v_int_new(n_layers-1)           ! preliminary new virtual potential temperature interface values
    real(wp)               :: solution_vector(n_layers-1+nsoillays) ! vector containing the solution of the linear problem to solve here
    real(wp)               :: rho_int_old(n_layers-1)               ! old interface mass density
    real(wp)               :: rho_int_expl(n_layers-1)              ! explicit interface mass density
    real(wp)               :: alpha_old(n_layers)                   ! alpha at the old time step
    real(wp)               :: beta_old(n_layers)                    ! beta at the old time step
    real(wp)               :: gamma_old(n_layers)                   ! gamma at the old time step
    real(wp)               :: alpha_new(n_layers)                   ! alpha at the new time step
    real(wp)               :: beta_new(n_layers)                    ! beta at the new time step
    real(wp)               :: gamma_new(n_layers)                   ! gamma at the new time step
    real(wp)               :: alpha(n_layers)                       ! alpha
    real(wp)               :: beta(n_layers)                        ! beta
    real(wp)               :: gammaa(n_layers)                      ! gamma
    real(wp)               :: rho_int_new                           ! new density interface value
    real(wp)               :: heat_flux_density_expl(nsoillays)     ! explicit heat_flux_density in the soil
    type(t_state), pointer :: state_new_used                        ! pointer to the state that is used as the new state in the calculation
    
    if (rk_step==1) then
      state_new_used => state_old
    else
      state_new_used => state_new
    endif
    
    ! the maximum temperature change induced by radiation between two radiation time steps in the uppermost soil layer
    max_rad_temp_change = 25._wp
    
    damping_start_height = klemp_begin_rel*toa
    
    ! partial derivatives new time step weight
    partial_deriv_new_time_step_weight = 0.5_wp
    
    ! calculating the sensible power flux density
    if (lsfc_sensible_heat_flux) then
      !$omp parallel do private(ji,t_gas_lowest_layer_old,t_gas_lowest_layer_new)
      do ji=1,n_cells
        
        ! gas temperature in the lowest layer
        t_gas_lowest_layer_old = (grid%exner_bg(ji,n_layers) + state_old%exner_pert(ji,n_layers)) &
        *(grid%theta_v_bg(ji,n_layers) + state_old%theta_v_pert(ji,n_layers))
        t_gas_lowest_layer_new = (grid%exner_bg(ji,n_layers) + state_new_used%exner_pert(ji,n_layers)) &
        *(grid%theta_v_bg(ji,n_layers) + state_new_used%theta_v_pert(ji,n_layers))
        
        ! converting the virtual temperature to the real temperature
        if (lmoist) then
          t_gas_lowest_layer_old = t_gas_lowest_layer_old &
                                             /(1._wp+state_old%rho(ji,n_layers,n_condensed_constituents+2) &
                                             /state_old%rho(ji,n_layers,n_condensed_constituents+1)*(m_d/m_v-1._wp))
          t_gas_lowest_layer_new = t_gas_lowest_layer_new &
                                             /(1._wp+state_new_used%rho(ji,n_layers,n_condensed_constituents+2) &
                                             /state_new_used%rho(ji,n_layers,n_condensed_constituents+1)*(m_d/m_v-1._wp))
        endif
        
        ! the sensible power flux density
        diag%power_flux_density_sens_soil(ji) = (grid%land_fraction(ji)+grid%lake_fraction(ji)) &
        *0.5_wp*c_d_v*(state_new_used%rho(ji,n_layers,n_condensed_constituents+1) &
        *(t_gas_lowest_layer_old - state_old%temperature_soil(ji,1)) &
        + state_old%rho(ji,n_layers,n_condensed_constituents+1) &
        *(t_gas_lowest_layer_new - state_new_used%temperature_soil(ji,1)))/diag%scalar_flux_resistance(ji)
        diag%power_flux_density_sens_sea(ji) = (1._wp-grid%land_fraction(ji)-grid%lake_fraction(ji)) &
        *0.5_wp*c_d_v*(state_new_used%rho(ji,n_layers,n_condensed_constituents+1) &
        *(t_gas_lowest_layer_old - diag%sst(ji)) &
        + state_old%rho(ji,n_layers,n_condensed_constituents+1) &
        *(t_gas_lowest_layer_new - diag%sst(ji)))/diag%scalar_flux_resistance(ji)
        
        ! contribution of sensible heat to rhotheta_v
        state_tend%rhotheta_v(ji,n_layers) = state_tend%rhotheta_v(ji,n_layers) &
        - grid%area_v(ji,n_levels)*(diag%power_flux_density_sens_soil(ji)+diag%power_flux_density_sens_sea(ji)) &
        /((grid%exner_bg(ji,n_layers) + state_new_used%exner_pert(ji,n_layers))*c_d_p)/grid%volume(ji,n_layers)
      enddo
      !$omp end parallel do
    endif
    
    ! loop over all columns
    !$omp parallel do private(ji,jl,damping_coeff,damping_prefactor,z_above_damping,soil_switch, &
    !$omp radiation_flux_density,resulting_temperature_change,c_vector,d_vector, &
    !$omp e_vector,r_vector,rho_expl,rhotheta_v_expl,theta_v_pert_expl,exner_pert_expl,theta_v_int_new, &
    !$omp solution_vector,rho_int_old,rho_int_expl,alpha_old,beta_old,gamma_old,alpha_new,beta_new, &
    !$omp gamma_new,alpha,beta,gammaa,rho_int_new,heat_flux_density_expl)
    do ji=1,n_cells
      
      ! determining wether soil needs to be calculated
      if (lprog_soil_temp .and. grid%land_fraction(ji)+grid%lake_fraction(ji)>0._wp) then
        soil_switch = 1
      else
        soil_switch = 0
      endif
      
      ! explicit quantities
      do jl=1,n_layers
        ! explicit density
        rho_expl(jl) = state_old%rho(ji,jl,n_condensed_constituents+1) &
        + dtime*state_tend%rho(ji,jl,n_condensed_constituents+1)
        ! explicit virtual potential temperature density
        rhotheta_v_expl(jl) = state_old%rhotheta_v(ji,jl) + dtime*state_tend%rhotheta_v(ji,jl)
        if (rk_step==1) then
          ! old time step partial derivatives of theta_v and Pi (divided by the volume)
          alpha(jl) = -state_old%rhotheta_v(ji,jl)/state_old%rho(ji,jl,n_condensed_constituents+1)**2 &
          /grid%volume(ji,jl)
          beta(jl) = 1._wp/state_old%rho(ji,jl,n_condensed_constituents+1)/grid%volume(ji,jl)
          gammaa(jl) = r_d/(c_d_v*state_old%rhotheta_v(ji,jl)) &
          *(grid%exner_bg(ji,jl) + state_old%exner_pert(ji,jl))/grid%volume(ji,jl)
        else
          ! old time step partial derivatives of theta_v and Pi
          alpha_old(jl) = -state_old%rhotheta_v(ji,jl)/state_old%rho(ji,jl,n_condensed_constituents+1)**2
          beta_old(jl) = 1._wp/state_old%rho(ji,jl,n_condensed_constituents+1)
          gamma_old(jl) = r_d/(c_d_v*state_old%rhotheta_v(ji,jl))*(grid%exner_bg(ji,jl)+state_old%exner_pert(ji,jl))
          ! new time step partial derivatives of theta_v and Pi
          alpha_new(jl) = -state_new_used%rhotheta_v(ji,jl)/state_new_used%rho(ji,jl,n_condensed_constituents+1)**2
          beta_new(jl) = 1._wp/state_new_used%rho(ji,jl,n_condensed_constituents+1)
          gamma_new(jl) = r_d/(c_d_v*state_new_used%rhotheta_v(ji,jl))*(grid%exner_bg(ji,jl)+state_new_used%exner_pert(ji,jl))
          ! interpolation in time and dividing by the volume
          alpha(jl) = ((1._wp - partial_deriv_new_time_step_weight)*alpha_old(jl) &
          + partial_deriv_new_time_step_weight*alpha_new(jl))/grid%volume(ji,jl)
          beta(jl) = ((1._wp - partial_deriv_new_time_step_weight)*beta_old(jl) &
          + partial_deriv_new_time_step_weight*beta_new(jl))/grid%volume(ji,jl)
          gammaa(jl) = ((1._wp - partial_deriv_new_time_step_weight)*gamma_old(jl) &
          + partial_deriv_new_time_step_weight*gamma_new(jl))/grid%volume(ji,jl)
        endif
        ! explicit virtual potential temperature perturbation
        theta_v_pert_expl(jl) = state_old%theta_v_pert(ji,jl) + dtime*grid%volume(ji,jl)* &
        (alpha(jl)*state_tend%rho(ji,jl,n_condensed_constituents+1) + beta(jl)*state_tend%rhotheta_v(ji,jl))
        ! explicit Exner pressure perturbation
        exner_pert_expl(jl) = state_old%exner_pert(ji,jl) &
                              + dtime*grid%volume(ji,jl)*gammaa(jl)*state_tend%rhotheta_v(ji,jl)
      enddo
      
      ! determining the interface values
      do jl=1,n_layers-1
        rho_int_old(jl) = 0.5_wp*(state_old%rho(ji,jl,n_condensed_constituents+1) &
                                + state_old%rho(ji,jl+1,n_condensed_constituents+1))
        rho_int_expl(jl) = 0.5_wp*(rho_expl(jl) + rho_expl(jl+1))
        theta_v_int_new(jl) = 0.5_wp*(state_new_used%rhotheta_v(ji,jl)/state_new_used%rho(ji,jl,n_condensed_constituents+1) &
        + state_new_used%rhotheta_v(ji,jl+1)/state_new_used%rho(ji,jl+1,n_condensed_constituents+1))
      enddo
      
      ! filling up the coefficient vectors
      do jl=1,n_layers-1
        ! main diagonal
        ! Klemp (2008) upper boundary layer
        z_above_damping = grid%z_vector_v(ji,jl+1) - damping_start_height
        if (z_above_damping<0._wp) then
          damping_coeff = 0._wp
        else
          damping_coeff = klemp_damp_max*sin(0.5_wp*M_PI*z_above_damping/(toa - damping_start_height))**2
        endif
        damping_prefactor(jl) = 1._wp + dtime*damping_coeff
        d_vector(jl) = -theta_v_int_new(jl)**2*(gammaa(jl) + gammaa(jl+1)) &
        + 0.5_wp*(grid%exner_bg(ji,jl) - grid%exner_bg(ji,jl+1)) &
        *(alpha(jl+1) - alpha(jl) + theta_v_int_new(jl)*(beta(jl+1) - beta(jl))) &
        - (grid%z_scalar(ji,jl) - grid%z_scalar(ji,jl+1))/(impl_thermo_weight*dtime**2*c_d_p*rho_int_old(jl)) &
        *(2._wp/grid%area_v(ji,jl+1) + dtime*state_old%wind_v(ji,jl+1)*0.5_wp &
        *(-1._wp/grid%volume(ji,jl) + 1._wp/grid%volume(ji,jl+1)))*damping_prefactor(jl)
        ! right hand side
        r_vector(jl) = -(state_old%wind_v(ji,jl+1) + dtime*state_tend%wind_v(ji,jl+1)) &
        *(grid%z_scalar(ji,jl) - grid%z_scalar(ji,jl+1))/(impl_thermo_weight*dtime**2*c_d_p) &
        + theta_v_int_new(jl)*(exner_pert_expl(jl) - exner_pert_expl(jl+1))/dtime &
        + 0.5_wp/dtime*(theta_v_pert_expl(jl) + theta_v_pert_expl(jl+1))*(grid%exner_bg(ji,jl) - grid%exner_bg(ji,jl+1)) &
        - (grid%z_scalar(ji,jl) - grid%z_scalar(ji,jl+1))/(impl_thermo_weight*dtime**2*c_d_p) &
        *state_old%wind_v(ji,jl+1)*rho_int_expl(jl)/rho_int_old(jl)
      enddo
      do jl=1,n_layers-2
        ! lower diagonal
        c_vector(jl) = theta_v_int_new(jl+1)*gammaa(jl+1)*theta_v_int_new(jl) &
        + 0.5_wp*(grid%exner_bg(ji,jl+1) - grid%exner_bg(ji,jl+2)) &
        *(alpha(jl+1) + beta(jl+1)*theta_v_int_new(jl)) &
        - (grid%z_scalar(ji,jl+1) - grid%z_scalar(ji,jl+2))/(impl_thermo_weight*dtime*c_d_p)*0.5_wp &
        *state_old%wind_v(ji,jl+2)/(grid%volume(ji,jl+1)*rho_int_old(jl+1))*damping_prefactor(jl+1)
        ! upper diagonal
        e_vector(jl) = theta_v_int_new(jl)*gammaa(jl+1)*theta_v_int_new(jl+1) &
        - 0.5_wp*(grid%exner_bg(ji,jl) - grid%exner_bg(ji,jl+1)) &
        *(alpha(jl+1) + beta(jl+1)*theta_v_int_new(jl+1)) &
        + (grid%z_scalar(ji,jl) - grid%z_scalar(ji,jl+1))/(impl_thermo_weight*dtime*c_d_p)*0.5_wp &
        *state_old%wind_v(ji,jl+1)/(grid%volume(ji,jl+1)*rho_int_old(jl))*damping_prefactor(jl)
      enddo
      
      ! soil components of the matrix
      if (soil_switch==1) then
        
        ! calculating the explicit part of the heat flux density
        do jl=1,nsoillays-1
          heat_flux_density_expl(jl) &
          = -grid%sfc_rho_c(ji)*grid%t_conduc_soil(ji)*(state_old%temperature_soil(ji,jl) &
          - state_old%temperature_soil(ji,jl+1))/(grid%z_soil_center(jl) - grid%z_soil_center(jl+1))
        enddo
        heat_flux_density_expl(nsoillays) = -grid%sfc_rho_c(ji)*grid%t_conduc_soil(ji)* &
        (state_old%temperature_soil(ji,nsoillays)-grid%t_const_soil(ji))/(2._wp*(grid%z_soil_center(nsoillays) - z_t_const))
        
        radiation_flux_density = diag%sfc_sw_in(ji) - diag%sfc_lw_out(ji)
        resulting_temperature_change = radiation_flux_density/((grid%z_soil_interface(1) - grid%z_soil_interface(2)) &
        *grid%sfc_rho_c(ji))*radiation_dtime
        if (abs(resulting_temperature_change)>max_rad_temp_change) then
          radiation_flux_density = max_rad_temp_change/abs(resulting_temperature_change)*radiation_flux_density
        endif
        
        ! calculating the explicit part of the temperature change
        r_vector(n_layers) &
        ! old temperature
        = state_old%temperature_soil(ji,1) &
        ! sensible heat flux
        + (diag%power_flux_density_sens_soil(ji) &
        ! latent heat flux
        + diag%power_flux_density_lat_lake(ji) &
        ! radiation
        + radiation_flux_density &
        ! heat conduction from below
        + 0.5_wp*heat_flux_density_expl(1)) &
        /((grid%z_soil_interface(1) - grid%z_soil_interface(2))*grid%sfc_rho_c(ji))*dtime
        
        ! loop over all soil layers below the first layer
        do jl=2,nsoillays
          
          r_vector(jl+n_layers-1) &
          ! old temperature
          = state_old%temperature_soil(ji,jl) &
          ! heat conduction from above
          + 0.5_wp*(-heat_flux_density_expl(jl-1) &
          ! heat conduction from below
          + heat_flux_density_expl(jl)) &
          /((grid%z_soil_interface(jl) - grid%z_soil_interface(jl+1))*grid%sfc_rho_c(ji))*dtime
          
        enddo
        
        ! the diagonal component
        do jl=1,nsoillays
          if (jl==1) then
            d_vector(jl+n_layers-1) = 1._wp + 0.5_wp*dtime*grid%sfc_rho_c(ji)*grid%t_conduc_soil(ji) &
            /((grid%z_soil_interface(jl) - grid%z_soil_interface(jl+1))*grid%sfc_rho_c(ji)) &
            *1._wp/(grid%z_soil_center(jl) - grid%z_soil_center(jl+1))
          elseif (jl==nsoillays) then
            d_vector(jl+n_layers-1) = 1._wp + 0.5_wp*dtime*grid%sfc_rho_c(ji)*grid%t_conduc_soil(ji) &
            /((grid%z_soil_interface(jl) - grid%z_soil_interface(jl+1))*grid%sfc_rho_c(ji)) &
            *1._wp/(grid%z_soil_center(jl-1) - grid%z_soil_center(jl))
          else
            d_vector(jl+n_layers-1) = 1._wp + 0.5_wp*dtime*grid%sfc_rho_c(ji)*grid%t_conduc_soil(ji) &
            /((grid%z_soil_interface(jl) - grid%z_soil_interface(jl+1))*grid%sfc_rho_c(ji)) &
            *(1._wp/(grid%z_soil_center(jl-1) - grid%z_soil_center(jl)) &
            + 1._wp/(grid%z_soil_center(jl) - grid%z_soil_center(jl+1)))
          endif
        enddo
        
        ! the off-diagonal components
        c_vector(n_layers-1) = 0._wp
        e_vector(n_layers-1) = 0._wp
        do jl=1,nsoillays-1
          c_vector(jl+n_layers-1) = -0.5_wp*dtime*grid%sfc_rho_c(ji)*grid%t_conduc_soil(ji) &
          /((grid%z_soil_interface(jl+1) - grid%z_soil_interface(jl+2))*grid%sfc_rho_c(ji)) &
          /(grid%z_soil_center(jl) - grid%z_soil_center(jl+1))
          e_vector(jl+n_layers-1) = -0.5_wp*dtime*grid%sfc_rho_c(ji)*grid%t_conduc_soil(ji) &
          /((grid%z_soil_interface(jl) - grid%z_soil_interface(jl+1))*grid%sfc_rho_c(ji)) &
          /(grid%z_soil_center(jl) - grid%z_soil_center(jl+1))
        enddo
        
      endif
      
      ! calling the algorithm to solve the system of linear equations
      call thomas_algorithm(c_vector,d_vector,e_vector,r_vector,solution_vector,n_layers-1+soil_switch*nsoillays)
      
      ! Writing the result into the new state.
      ! --------------------------------------
      
      ! mass density
      do jl=1,n_layers
        if (jl==1) then
          state_new%rho(ji,jl,n_condensed_constituents+1) = rho_expl(jl) + dtime*(solution_vector(jl))/grid%volume(ji,jl)
        elseif (jl==n_layers) then
          state_new%rho(ji,jl,n_condensed_constituents+1) = rho_expl(jl) + dtime*(-solution_vector(jl-1))/grid%volume(ji,jl)
        else
          state_new%rho(ji,jl,n_condensed_constituents+1) &
          = rho_expl(jl) + dtime*(-solution_vector(jl-1) + solution_vector(jl))/grid%volume(ji,jl)
        endif
      enddo
      ! virtual potential temperature density
      do jl=1,n_layers
        if (jl==1) then
          state_new%rhotheta_v(ji,jl) &
          = rhotheta_v_expl(jl) + dtime*(theta_v_int_new(jl)*solution_vector(jl))/grid%volume(ji,jl)
        elseif (jl==n_layers) then
          state_new%rhotheta_v(ji,jl) &
          = rhotheta_v_expl(jl) + dtime*(-theta_v_int_new(jl-1)*solution_vector(jl-1))/grid%volume(ji,jl)
        else
          state_new%rhotheta_v(ji,jl) &
          = rhotheta_v_expl(jl) + dtime*(-theta_v_int_new(jl-1)*solution_vector(jl-1) + theta_v_int_new(jl)*solution_vector(jl)) &
          /grid%volume(ji,jl)
        endif
      enddo
      ! vertical velocity
      do jl=2,n_layers
        rho_int_new = 0.5_wp*(state_new%rho(ji,jl-1,n_condensed_constituents+1) &
                                      + state_new%rho(ji,jl,n_condensed_constituents+1))
        state_new%wind_v(ji,jl)  = (2._wp*solution_vector(jl-1)/grid%area_v(ji,jl) &
                                       - rho_int_new*state_old%wind_v(ji,jl))/rho_int_old(jl-1)
      enddo
      ! virtual potential temperature perturbation
      do jl=1,n_layers
        state_new%theta_v_pert(ji,jl) = state_new%rhotheta_v(ji,jl)/state_new%rho(ji,jl,n_condensed_constituents+1) &
                                           - grid%theta_v_bg(ji,jl)
      enddo
      ! Exner pressure perturbation
      do jl=1,n_layers
        state_new%exner_pert(ji,jl) = state_old%exner_pert(ji,jl) + grid%volume(ji,jl) &
                                         *gammaa(jl)*(state_new%rhotheta_v(ji,jl) - state_old%rhotheta_v(ji,jl))
      enddo
      
      ! soil temperature
      if (soil_switch==1) then
        do jl=1,nsoillays
          state_new%temperature_soil(ji,jl) = solution_vector(n_layers-1+jl)
        enddo
      endif
      
    enddo ! end of the column loop
    !$omp end parallel do
    
  end subroutine three_band_solver_ver_waves
  
  subroutine three_band_solver_gen_densities(state_old,state_new,state_tend,diag,grid,rk_step)
    
    ! This subroutine contains the vertical advection of mass densities (of tracers) with 3-band matrices.
    
    type(t_state), intent(in)    :: state_old  ! state variables at the old time step
    type(t_state), intent(inout) :: state_new  ! state variables at the new time step
    type(t_state), intent(in)    :: state_tend ! state containing the tendencies of the state variables
    type(t_diag),  intent(inout) :: diag       ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid       ! grid quantities
    integer,       intent(in)    :: rk_step    ! predictor-corrector substep
    
    ! local variables
    integer  :: ji                                        ! cell index
    integer  :: jl                                        ! layer index
    integer  :: jc                                        ! constituent index
    real(wp) :: impl_weight                               ! implicit time stepping weight
    real(wp) :: expl_weight                               ! explicit time stepping weight
    real(wp) :: density_old_at_interface                  ! density in a level at the old substep old density in a level
    real(wp) :: temperature_old_at_interface              ! temperature in a level at the old substep
    real(wp) :: c_vector(n_layers-1)                      ! vector for solving the system of linear equations
    real(wp) :: d_vector(n_layers)                        ! vector for solving the system of linear equations
    real(wp) :: e_vector(n_layers-1)                      ! vector for solving the system of linear equations
    real(wp) :: r_vector(n_layers)                        ! vector for solving the system of linear equations
    real(wp) :: vertical_flux_vector_impl(n_layers-1)     ! vertical flux at the new time step
    real(wp) :: vertical_flux_vector_rhs(n_layers-1)      ! vertical flux at the old time step
    real(wp) :: vertical_enthalpy_flux_vector(n_layers-1) ! vertical enthalpy flux density vector
    real(wp) :: solution_vector(n_layers)                 ! solution of the system of linear equations
    real(wp) :: v_fall_upper                              ! fall velocity of a hydrometeor particle in the lower grid box
    real(wp) :: v_fall_lower                              ! fall velocity of a hydrometeor particle in the upper grid box
    real(wp) :: v_fall(n_layers)                          ! fall velocity of a hydrometeor particle in a level
    
    impl_weight = 0.5_wp
    expl_weight = 1._wp - impl_weight
    
    ! loop over all relevant constituents
    do jc=1,n_constituents
      
      ! This is done for all tracers apart from the main gaseous constituent.
      if (jc/=n_condensed_constituents+1) then
        ! loop over all columns
        !$omp parallel do private(ji,jl,density_old_at_interface,v_fall_upper,v_fall_lower,v_fall, &
        !$omp temperature_old_at_interface,c_vector,d_vector,e_vector,r_vector,vertical_flux_vector_impl, &
        !$omp vertical_flux_vector_rhs,vertical_enthalpy_flux_vector,solution_vector)
        do ji=1,n_cells
          
          ! diagnozing the vertical fluxes
          do jl=1,n_layers-1
            ! resetting the vertical enthalpy flux density divergence
            if (rk_step==1 .and. jc==1) then
              diag%condensates_sediment_heat(ji,jl) = 0._wp
            endif
            vertical_flux_vector_impl(jl) = state_old%wind_v(ji,jl+1)
            vertical_flux_vector_rhs(jl) = state_new%wind_v(ji,jl+1)
            ! preparing the vertical interpolation
            ! For condensed constituents, a fall velocity must be added.
            ! precipitation
            ! snow
            if (jc==1) then
              v_fall_upper = v_fall_solid(state_old,diag,snow_particles_radius(),ji,jl)
              v_fall_lower = v_fall_solid(state_old,diag,snow_particles_radius(),ji,jl+1)
              v_fall(jl) = 0.5_wp*(v_fall_upper + v_fall_lower)
            ! rain
            elseif (jc==2) then
              v_fall_upper = v_fall_liquid(state_old,diag,grid,diag%a_rain(ji,jl),ji,jl)
              v_fall_lower = v_fall_liquid(state_old,diag,grid,diag%a_rain(ji,jl+1),ji,jl+1)
              v_fall(jl) = 0.5_wp*(v_fall_upper + v_fall_lower)
            ! ice clouds
            elseif (jc==3) then
              v_fall_upper = v_fall_solid(state_old,diag,ice_particles_radius(),ji,jl)
              v_fall_lower = v_fall_solid(state_old,diag,ice_particles_radius(),ji,jl+1)
              v_fall(jl) = 0.5_wp*(v_fall_upper + v_fall_lower)
            ! water clouds
            elseif (jc==4) then
              v_fall_upper = v_fall_liquid(state_old,diag,grid,cloud_droplets_radius(),ji,jl)
              v_fall_lower = v_fall_liquid(state_old,diag,grid,cloud_droplets_radius(),ji,jl+1)
              v_fall(jl) = 0.5_wp*(v_fall_upper + v_fall_lower)
            else
              v_fall(jl) = 0._wp
            endif
            vertical_flux_vector_impl(jl) = vertical_flux_vector_impl(jl) - v_fall(jl)
            vertical_flux_vector_rhs(jl) = vertical_flux_vector_rhs(jl) - v_fall(jl)
            ! multiplying the vertical velocity by the vertical areas
            vertical_flux_vector_impl(jl) = grid%area_v(ji,jl+1)*vertical_flux_vector_impl(jl)
            vertical_flux_vector_rhs(jl) = grid%area_v(ji,jl+1)*vertical_flux_vector_rhs(jl)
            ! old density and temperature at the interface
            if (vertical_flux_vector_rhs(jl)>=0._wp) then
              density_old_at_interface = state_old%rho(ji,jl+1,jc)
              temperature_old_at_interface = diag%temperature(ji,jl+1)
            else
              density_old_at_interface = state_old%rho(ji,jl,jc)
              temperature_old_at_interface = diag%temperature(ji,jl)
            endif
            vertical_flux_vector_rhs(jl) = density_old_at_interface*vertical_flux_vector_rhs(jl)
            vertical_enthalpy_flux_vector(jl) = c_p_cond(jc,temperature_old_at_interface) &
                                                *temperature_old_at_interface*vertical_flux_vector_rhs(jl)
          enddo
          if (rk_step==1 .and. jc==1) then
            diag%condensates_sediment_heat(ji,n_layers) = 0._wp
          endif
          
          ! sink velocities at the surface
          ! ice
          if (jc==1) then
            v_fall(n_layers) = v_fall_solid(state_old,diag,snow_particles_radius(),ji,n_layers)
          ! rain
          elseif (jc==2) then
            v_fall(n_layers) = v_fall_liquid(state_old,diag,grid,diag%a_rain(ji,jl),ji,n_layers)
          ! ice clouds
          elseif (jc==3) then
            v_fall(n_layers) = v_fall_solid(state_old,diag,ice_particles_radius(),ji,n_layers)
          ! water clouds
          elseif (jc==4) then
            v_fall(n_layers) = v_fall_liquid(state_old,diag,grid,cloud_droplets_radius(),ji,n_layers)
          endif
          
          ! Now we proceed to solving the vertical tridiagonal problem.
          
          ! filling up the original vectors
          do jl=1,n_layers-1
            if (vertical_flux_vector_impl(jl)>=0._wp) then
              c_vector(jl) = 0._wp
              e_vector(jl) = -impl_weight*dtime/grid%volume(ji,jl)*vertical_flux_vector_impl(jl)
            else
              c_vector(jl) = impl_weight*dtime/grid%volume(ji,jl+1)*vertical_flux_vector_impl(jl)
              e_vector(jl) = 0._wp
            endif
          enddo
          do jl=1,n_layers
            if (jl==1) then
              if (vertical_flux_vector_impl(1)>=0._wp) then
                d_vector(jl) = 1._wp
              else
                d_vector(jl) = 1._wp - impl_weight*dtime/grid%volume(ji,jl)*vertical_flux_vector_impl(1)
              endif
            elseif (jl==n_layers) then
              if (vertical_flux_vector_impl(jl-1)>=0._wp) then
                d_vector(jl) = 1._wp + impl_weight*dtime/grid%volume(ji,jl)*vertical_flux_vector_impl(jl-1)
              else
                d_vector(jl) = 1._wp
              endif
              ! precipitation
              ! snow
              if (jc<=n_condensed_constituents/4) then
                d_vector(jl) = d_vector(jl) + impl_weight*v_fall(jl)*dtime*grid%area_v(ji,n_levels)/grid%volume(ji,jl)
              ! rain
              elseif (jc<=n_condensed_constituents/2) then
                d_vector(jl) = d_vector(jl) + impl_weight*v_fall(jl)*dtime*grid%area_v(ji,n_levels)/grid%volume(ji,jl)
              ! clouds
              elseif (jc<=n_condensed_constituents) then
                d_vector(jl) = d_vector(jl) + impl_weight*v_fall(jl)*dtime*grid%area_v(ji,n_levels) &
                               /grid%volume(ji,jl)
              endif
            else
              d_vector(jl) = 1._wp
              if (vertical_flux_vector_impl(jl-1)>=0._wp) then
                d_vector(jl) = d_vector(jl) + impl_weight*dtime/grid%volume(ji,jl)*vertical_flux_vector_impl(jl-1)
              endif
              if (vertical_flux_vector_impl(jl)<0._wp) then
                d_vector(jl) = d_vector(jl) - impl_weight*dtime/grid%volume(ji,jl)*vertical_flux_vector_impl(jl)
              endif
            endif
            ! the explicit component
            ! mass densities
            r_vector(jl) = state_old%rho(ji,jl,jc) + dtime*state_tend%rho(ji,jl,jc)
            ! adding the explicit part of the vertical flux divergence
            if (jl==1) then
              r_vector(jl) = r_vector(jl) + expl_weight*dtime*vertical_flux_vector_rhs(jl)/grid%volume(ji,jl)
              if (rk_step==1 .and. jc<=n_condensed_constituents) then
                diag%condensates_sediment_heat(ji,jl) = diag%condensates_sediment_heat(ji,jl) &
                + vertical_enthalpy_flux_vector(jl)/grid%volume(ji,jl)
              endif
            elseif (jl==n_layers) then
              r_vector(jl) = r_vector(jl) - expl_weight*dtime*vertical_flux_vector_rhs(jl-1)/grid%volume(ji,jl)
              if (rk_step==1 .and. jc<=n_condensed_constituents) then
                diag%condensates_sediment_heat(ji,jl) = diag%condensates_sediment_heat(ji,jl) &
                - vertical_enthalpy_flux_vector(jl-1)/grid%volume(ji,jl)
              endif
              ! precipitation
              ! snow
              if (jc<=n_condensed_constituents/4) then
                r_vector(jl) = r_vector(jl) - expl_weight*v_fall(jl)*dtime &
                *state_old%rho(ji,n_layers,jc)*grid%area_v(ji,n_levels)/grid%volume(ji,jl)
                if (rk_step==1) then
                  diag%condensates_sediment_heat(ji,jl) = diag%condensates_sediment_heat(ji,jl) &
                  - v_fall(jl) &
                  *diag%temperature(ji,n_layers)*c_p_cond(jc,diag%temperature(ji,n_layers)) &
                  *state_old%rho(ji,n_layers,jc)*grid%area_v(ji,n_levels)/grid%volume(ji,jl)
                endif
              ! rain
              elseif (jc<=n_condensed_constituents/2) then
                r_vector(jl) = r_vector(jl) - expl_weight*v_fall(jl)*dtime &
                *state_old%rho(ji,n_layers,jc)*grid%area_v(ji,n_levels)/grid%volume(ji,jl)
                if (rk_step==1) then
                  diag%condensates_sediment_heat(ji,jl) = diag%condensates_sediment_heat(ji,jl) &
                  - v_fall(jl) &
                  *diag%temperature(ji,n_layers)*c_p_cond(jc,diag%temperature(ji,n_layers)) &
                  *state_old%rho(ji,n_layers,jc)*grid%area_v(ji,n_levels)/grid%volume(ji,jl)
                endif
              ! clouds
              elseif (jc<=n_condensed_constituents) then
                r_vector(jl) = r_vector(jl) - expl_weight*v_fall(jl)*dtime &
                *state_old%rho(ji,n_layers,jc)*grid%area_v(ji,n_levels)/grid%volume(ji,jl)
                if (rk_step==1) then
                  diag%condensates_sediment_heat(ji,jl) = diag%condensates_sediment_heat(ji,jl) &
                  - v_fall(jl) &
                  *diag%temperature(ji,n_layers)*c_p_cond(jc,diag%temperature(ji,n_layers)) &
                  *state_old%rho(ji,n_layers,jc)*grid%area_v(ji,n_levels)/grid%volume(ji,jl)
                endif
              endif
            else
              r_vector(jl) = r_vector(jl) + expl_weight*dtime &
              *(-vertical_flux_vector_rhs(jl-1) + vertical_flux_vector_rhs(jl))/grid%volume(ji,jl)
              if (rk_step==1 .and. jc<=n_condensed_constituents) then
                diag%condensates_sediment_heat(ji,jl) = diag%condensates_sediment_heat(ji,jl) &
                + (-vertical_enthalpy_flux_vector(jl-1) + vertical_enthalpy_flux_vector(jl))/grid%volume(ji,jl)
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
            state_new%rho(ji,jl,jc) = solution_vector(jl)
          enddo
        enddo ! horizontal index
        !$omp end parallel do
      endif
    enddo ! constituent
    
  end subroutine three_band_solver_gen_densities
  
  subroutine thomas_algorithm(c_vector,d_vector,e_vector,r_vector,solution_vector,solution_length)
    
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
  
end module mo_column_solvers











