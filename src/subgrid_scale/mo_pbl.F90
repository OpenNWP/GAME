! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_pbl

  ! In this module, quantities referring to the planetary boundary layer are computed.

  use mo_constants,        only: EPSILON_SECURITY,M_PI,gravity,p_0,c_d_p,r_d
  use mo_definitions,      only: wp,t_grid,t_state,t_diag
  use mo_diff_nml,         only: h_prandtl,karman
  use mo_grid_nml,         only: n_cells,n_layers,n_edges,n_levels
  use mo_surface_nml,      only: lprog_soil_temp,pbl_scheme
  use mo_constituents_nml, only: n_constituents
  use mo_grid_setup,       only: dtime
  use mo_derived,          only: gas_constant_diagnostics
  
  implicit none
  
  contains
  
  subroutine pbl_wind_tendency(state,diag,grid)
  
    ! This subroutine computes the interaction of the horizontal wind with the surface.
  
    type(t_state), intent(in)    :: state ! state
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid properties
  
    ! local variables
    integer  :: ji                          ! edge index
    integer  :: jl                          ! layer index
    real(wp) :: flux_resistance             ! momentum flux resistance value
    real(wp) :: wind_speed_lowest_layer     ! 
    real(wp) :: z_agl                       ! height of a gridpoint above ground level
    real(wp) :: layer_thickness             ! 
    real(wp) :: roughness_length_value      ! 
    real(wp) :: monin_obukhov_length_value  ! 
    real(wp) :: wind_rescale_factor         ! 
    real(wp) :: bndr_lr_visc_max            ! 
    real(wp) :: sigma_b                     ! 
    real(wp) :: standard_vert_lapse_rate    ! 
    real(wp) :: exner_from                  ! Exner pressure at the from-gridpoint
    real(wp) :: exner_to                    ! Exner pressure at the to-gridpoint
    real(wp) :: pressure_from               ! pressure at the from-gridpoint
    real(wp) :: pressure_to                 ! pressure at the to-gridpoint
    real(wp) :: pressure                    ! air pressure
    real(wp) :: temp_lowest_layer           ! temperature in the lowest layer
    real(wp) :: pressure_value_lowest_layer ! pressure in the lowest layer (helper variable)
    real(wp) :: temp_surface                ! temperature at the surface
    real(wp) :: surface_p_factor            ! factor used for computing the sea surface pressure
    real(wp) :: pressure_sfc_from           ! surface pressure at the from-gridpoint
    real(wp) :: pressure_sfc_to             ! surface pressure at the to-gridpoint
    real(wp) :: pressure_sfc                ! surface pressure at the edge
    real(wp) :: sigma                       ! pressure/(surface pressure)
    
    ! NWP case
    if (pbl_scheme==1) then
      !$omp parallel do private(ji,flux_resistance,wind_speed_lowest_layer,z_agl, &
      !$omp layer_thickness,monin_obukhov_length_value,wind_rescale_factor,roughness_length_value)
      do ji=1,n_edges
      
        ! averaging some quantities to the vector point
        wind_speed_lowest_layer = 0.5_wp*((diag%v_squared(grid%from_cell(ji),n_layers))**0.5_wp &
                                        + (diag%v_squared(grid%to_cell(ji),n_layers))**0.5_wp)
        z_agl = grid%z_vector_h(ji,n_layers) - 0.5_wp*(grid%z_vector_v(grid%from_cell(ji),n_levels) &
                                                     + grid%z_vector_v(grid%to_cell(ji),n_levels))
        layer_thickness = 0.5_wp*(grid%layer_thickness(grid%from_cell(ji),n_layers) &
                                + grid%layer_thickness(grid%to_cell(ji),n_layers))
        roughness_length_value = 0.5_wp*(grid%roughness_length(grid%from_cell(ji)) + grid%roughness_length(grid%to_cell(ji)))
        monin_obukhov_length_value = 0.5_wp*(diag%monin_obukhov_length(grid%from_cell(ji)) &
                                           + diag%monin_obukhov_length(grid%to_cell(ji)))
      
        ! calculating the flux resistance at the vector point
        flux_resistance = momentum_flux_resistance(wind_speed_lowest_layer, &
                                                   z_agl,roughness_length_value,monin_obukhov_length_value)
      
        ! rescaling the wind if the lowest wind vector is above the height of the Prandtl layer
        wind_rescale_factor = 1._wp
        if (z_agl>h_prandtl) then
          wind_rescale_factor = log(h_prandtl/roughness_length_value)/log(z_agl/roughness_length_value)
        endif
      
        ! adding the momentum flux into the surface as an acceleration
        diag%friction_acc_h(ji,n_layers) = diag%friction_acc_h(ji,n_layers) &
                                           - wind_rescale_factor*state%wind_h(ji,n_layers)/flux_resistance/layer_thickness
      enddo
      !$omp end parallel do
    endif
  
    ! This is the explicit friction ansatz in the boundary layer from the Held-Suarez (1994) test case.
    if (pbl_scheme==2) then
      ! some parameters
      bndr_lr_visc_max = 1._wp/86400._wp ! maximum friction coefficient in the boundary layer
      sigma_b = 0.7_wp ! boundary layer height in sigma-p coordinates
      standard_vert_lapse_rate = 0.0065_wp
      !$omp parallel do private(ji,jl,exner_from,exner_to,pressure_from,pressure_to,pressure, &
      !$omp temp_lowest_layer,pressure_value_lowest_layer,temp_surface,surface_p_factor, &
      !$omp pressure_sfc_from,pressure_sfc_to,pressure_sfc,sigma)
      do jl=1,n_layers
        do ji=1,n_edges
          ! calculating the pressure at the horizontal vector point
          exner_from = grid%exner_bg(grid%from_cell(ji),jl) + state%exner_pert(grid%from_cell(ji),jl)
          exner_to = grid%exner_bg(grid%to_cell(ji),jl) + state%exner_pert(grid%to_cell(ji),jl)
          pressure_from = p_0*exner_from**(c_d_p/r_d)
          pressure_to = p_0*exner_to**(c_d_p/r_d)
          pressure = 0.5_wp*(pressure_from+pressure_to)
        
          ! calculating the surface pressure at the horizontal vector point
          ! calculating the surface pressure at the from scalar point
          temp_lowest_layer = diag%temperature(grid%from_cell(ji),n_layers)
          exner_from = grid%exner_bg(grid%from_cell(ji),n_layers) + state%exner_pert(grid%from_cell(ji),n_layers)
          pressure_value_lowest_layer = p_0*exner_from**(c_d_p/r_d)
          temp_surface = temp_lowest_layer  + standard_vert_lapse_rate*(grid%z_scalar(grid%from_cell(ji),n_layers) &
                                                                      - grid%z_vector_v(grid%from_cell(ji),n_levels))
          surface_p_factor = (1._wp - (temp_surface-temp_lowest_layer)/temp_surface) &
                             **(grid%gravity_m_v(grid%from_cell(ji),n_layers)/ &
                             (gas_constant_diagnostics(state%rho,grid%from_cell(ji),n_layers)*standard_vert_lapse_rate))
          pressure_sfc_from = pressure_value_lowest_layer/surface_p_factor
          ! calculating the surface pressure at the to scalar point
          temp_lowest_layer = diag%temperature(grid%to_cell(ji),n_layers)
          exner_to = grid%exner_bg(grid%to_cell(ji),n_layers) + state%exner_pert(grid%to_cell(ji),n_layers)
          pressure_value_lowest_layer = p_0*exner_to**(c_d_p/r_d)
          temp_surface = temp_lowest_layer + standard_vert_lapse_rate*(grid%z_scalar(grid%to_cell(ji),n_layers) &
                                                                     - grid%z_vector_v(grid%to_cell(ji),n_levels))
          surface_p_factor = (1._wp - (temp_surface-temp_lowest_layer)/temp_surface) &
                             **(grid%gravity_m_v(grid%to_cell(ji),n_layers)/ &
                             (gas_constant_diagnostics(state%rho,grid%to_cell(ji),n_layers)*standard_vert_lapse_rate))
          pressure_sfc_to = pressure_value_lowest_layer/surface_p_factor
          ! averaging the surface pressure to the vector point
          pressure_sfc = 0.5_wp*(pressure_sfc_from+pressure_sfc_to)
        
          ! calculating sigma
          sigma = pressure/pressure_sfc
          ! finally calculating the friction acceleration
          diag%friction_acc_h(ji,jl) = diag%friction_acc_h(ji,jl) &
                                       - bndr_lr_visc_max*max(0._wp,(sigma-sigma_b)/(1._wp-sigma_b))*state%wind_h(ji,jl)
        enddo
      enddo
      !$omp end parallel do
    endif
    
  end subroutine
  
  subroutine update_sfc_turb_quantities(state,diag,grid)
    
    ! This subroutine updates surface-related turbulence quantities.
    
    type(t_state), intent(in)    :: state ! state
    type(t_diag),  intent(inout) :: diag  ! type containing diagnostic quantities
    type(t_grid),  intent(inout) :: grid  ! grid properties
    
    ! local variables
    integer  :: ji                       ! cell index
    real(wp) :: u_lowest_layer           ! 
    real(wp) :: u10                      ! 
    real(wp) :: z_agl                    ! 
    real(wp) :: theta_v_lowest_layer     ! 
    real(wp) :: theta_v_second_layer     ! 
    real(wp) :: dz                       ! 
    real(wp) :: dtheta_v_dz              ! 
    real(wp) :: w_pert                   ! 
    real(wp) :: theta_v_pert_value       ! 
    real(wp) :: w_pert_theta_v_pert_avg  ! 
    real(wp) :: w_theta_v_corr           ! 
    
    ! semi-empirical coefficient
    w_theta_v_corr = 0.2_wp
    
    ! loop over all cells in the lowest layer
    !$omp parallel do private(ji,u_lowest_layer,u10,z_agl,theta_v_lowest_layer,theta_v_second_layer, &
    !$omp dz,dtheta_v_dz,w_pert,theta_v_pert_value,w_pert_theta_v_pert_avg)
    do ji=1,n_cells
    
      ! computing the height of the gridpoint above the surface
      z_agl = grid%z_scalar(ji,n_layers) - grid%z_vector_v(ji,n_levels)
      
      ! wind speed in the lowest layer
      u_lowest_layer = diag%v_squared(ji,n_layers)**0.5_wp
      
      ! calculating the 10 m wind velocity from the logarithmic wind profile
      u10 = u_lowest_layer*log(10._wp/grid%roughness_length(ji))/log(z_agl/grid%roughness_length(ji))
      
      ! only over the sea the roughness length is time-dependant (because of the waves)
      if (grid%is_land(ji)==0) then
        ! calculating the roughness length fom the wind velocity
        grid%roughness_length(ji) = roughness_length_from_u10_sea(u10)
      endif
      
      ! updating the roughness velocity
      diag%roughness_velocity(ji) = calc_roughness_velocity(u_lowest_layer,z_agl,grid%roughness_length(ji))
      
      ! theta_v in the lowest layer
      theta_v_lowest_layer = grid%theta_v_bg(ji,n_layers) + state%theta_v_pert(ji,n_layers)
      ! theta_v in the second-lowest layer
      theta_v_second_layer = grid%theta_v_bg(ji,n_layers-1) + state%theta_v_pert(ji,n_layers-1)
      
      ! delta z
      dz = grid%z_scalar(ji,n_layers-1)-grid%z_scalar(ji,n_layers)
      
      ! vertical gradient of theta_v
      dtheta_v_dz = (theta_v_second_layer-theta_v_lowest_layer)/dz
      
      ! the perturbation of the vertical velocity is assumed to be proportional to the 10 m wind speed
      ! times a stability-dependant factor
      w_pert = u10*max(0.001_wp,0.02_wp*(1._wp-dtheta_v_dz/0.01_wp))
      theta_v_pert_value = -0.2_wp*dtime*w_pert*dtheta_v_dz
      w_pert_theta_v_pert_avg = w_theta_v_corr*w_pert*theta_v_pert_value
      
      ! security
      if (abs(w_pert_theta_v_pert_avg)<EPSILON_SECURITY) then
        w_pert_theta_v_pert_avg = EPSILON_SECURITY
      endif
      
      ! finally computing the Monin-Obukhov length
      diag%monin_obukhov_length(ji) = -theta_v_lowest_layer*diag%roughness_velocity(ji)**3/(karman*gravity*w_pert_theta_v_pert_avg)
    enddo
    !$omp end parallel do
    
    ! updating the surface flux resistance acting on scalar quantities (moisture and sensible heat)
    if (lprog_soil_temp) then
      !$omp parallel do private(ji)
      do ji=1,n_cells
        diag%scalar_flux_resistance(ji) = calc_scalar_flux_resistance(diag%roughness_velocity(ji), &
                                          grid%z_scalar(ji,n_layers)-grid%z_vector_v(ji,n_levels), &
                                          grid%roughness_length(ji),diag%monin_obukhov_length(ji))
      enddo
      !$omp end parallel do
    endif
    
  end subroutine update_sfc_turb_quantities

  function roughness_length_from_u10_sea(u10)
    
    ! This function returns the roughness length as a function of the mean wind speed at 10 m above a fully developed sea.
    ! refer to Stensrud,Parameterization schemes (2007), p.130

    real(wp), intent(in) :: u10                           ! wind speed 10 m above the surface (m/s)e
    real(wp)             :: roughness_length_from_u10_sea ! result

    ! local variables
    real(wp) :: swh        ! significant wave height
    real(wp) :: period     ! mean sea surface wave period
    real(wp) :: wavelength ! mean sea surface wave length

    ! empirically determined formula for the SWH
    swh = 0.0248_wp*u10**2

    ! empirically determined period of the waves
    period = 0.729_wp*u10

    ! deep-water gravity waves
    wavelength = gravity*period**2/(2._wp*M_PI)

    ! final result
    roughness_length_from_u10_sea = 1200._wp*swh*(swh/max(wavelength,EPSILON_SECURITY))**4.5_wp

    ! avoid too small values for stability
    roughness_length_from_u10_sea = max(0.0001_wp,roughness_length_from_u10_sea)
    
  end function roughness_length_from_u10_sea

  function calc_scalar_flux_resistance(roughness_velocity_value,z_agl,roughness_length_value,monin_obukhov_length_value)

    ! This function returns the surface flux resistance for scalar quantities.
    
    real(wp), intent(in) :: roughness_velocity_value    ! roughness length
    real(wp), intent(in) :: z_agl                       ! height above ground level of the lowest layer
    real(wp), intent(in) :: roughness_length_value      ! roughness length
    real(wp), intent(in) :: monin_obukhov_length_value  ! Monin-Obukhov length
    real(wp)             :: calc_scalar_flux_resistance ! result (m)
    
    ! local variables
    real(wp) :: used_vertical_height ! vertical height (above ground level) actually used in the calculation
    
    ! height of the prandtl layer
    used_vertical_height = min(z_agl,h_prandtl)
    
    calc_scalar_flux_resistance = 1._wp/(karman*roughness_velocity_value) &
    ! neutral conditions
    *(log(used_vertical_height/roughness_length_value) &
    ! non-neutral conditions
    - psi_h(used_vertical_height,monin_obukhov_length_value) &
    ! interfacial sublayer
    + log(7._wp))
    
    ! limitting the result for security
    if (calc_scalar_flux_resistance<dtime/z_agl) then
      calc_scalar_flux_resistance = dtime/z_agl
    endif
    
  end function calc_scalar_flux_resistance

  function momentum_flux_resistance(wind_h_lowest_layer,z_agl,roughness_length_value,monin_obukhov_length_value)

    ! This function returns the surface flux resistance for momentum.

    real(wp), intent(in) :: wind_h_lowest_layer        ! horizontal wind speed in the lowest layer
    real(wp), intent(in) :: z_agl                      ! height above ground level
    real(wp), intent(in) :: roughness_length_value     ! roughness length
    real(wp), intent(in) :: monin_obukhov_length_value ! Monin-Obukhov length
    real(wp)             :: momentum_flux_resistance   ! result

    ! local variables
    real(wp) :: used_vertical_height ! vertical height (above ground level) actually used in the calculation

    ! height of the prandtl layer
    used_vertical_height = min(z_agl,h_prandtl)

    momentum_flux_resistance = 1._wp/(karman*calc_roughness_velocity(wind_h_lowest_layer,z_agl,roughness_length_value)) &
    ! neutral conditions
    *(log(used_vertical_height/roughness_length_value) &
    ! non-neutral conditions
    - psi_m(used_vertical_height,monin_obukhov_length_value))

    ! limitting the result for security
    if (momentum_flux_resistance<dtime/z_agl) then
      momentum_flux_resistance = dtime/z_agl
    endif

  end function momentum_flux_resistance

  function calc_roughness_velocity(wind_speed,z_agl,roughness_length_value)

    ! This function returns the roughness velocity.
    
    real(wp), intent(in) :: wind_speed              ! wind speed at a certain height
    real(wp), intent(in) :: z_agl                   ! height at which the wind speed is valid
    real(wp), intent(in) :: roughness_length_value  ! roughness length at this point
    real(wp)             :: calc_roughness_velocity ! the result

    ! local variables
    real(wp) :: denominator ! helper variable

    denominator = log(z_agl/roughness_length_value)

    ! security
    if (abs(denominator)<EPSILON_SECURITY) then
      denominator = EPSILON_SECURITY
    endif

    calc_roughness_velocity = wind_speed*karman/denominator

    calc_roughness_velocity = max(EPSILON_SECURITY,calc_roughness_velocity)

  end function calc_roughness_velocity

  function psi_h(z_eff,l)

    ! This is a helper function for the correction to the surface scalar flux resistance for non-neutral conditions.
    
    real(wp), intent(in) :: z_eff ! effective height above the surface
    real(wp), intent(in) :: l     ! Monin-Obukhov length
    real(wp)             :: psi_h ! the value of the helper function
    
    ! local variables
    real(wp) :: x       ! helper variable
    real(wp) :: l_local ! used to avoid l==0
    
    ! avoiding l==0
    l_local = l
    if (abs(l_local)<EPSILON_SECURITY) then
      l_local = EPSILON_SECURITY
    endif
    
    ! unstable conditions
    if (l_local<0._wp) then
      x = (1._wp-15._wp*z_eff/l_local)**0.25_wp
      psi_h = 2._wp*log((1._wp + x**2)/2._wp)
    ! neutral and stable conditions
    else
      psi_h = -4.7_wp*z_eff/l_local
    endif
    
  end function psi_h

  function psi_m(z_eff,l)

    ! This is a helper function for the correction to the surface momentum flux resistance for non-neutral conditions.
    
    real(wp), intent(in) :: z_eff ! effective height above the surface (m)
    real(wp), intent(in) :: l     ! Monin-Obukhov length (m)
    real(wp)             :: psi_m ! the value of the helper function
    
    ! local variables
    real(wp) :: x       ! helper variable
    real(wp) :: l_local ! used to avoid l==0
    
    ! avoiding l==0
    l_local = l
    if (abs(l_local)<EPSILON_SECURITY) then
      l_local = EPSILON_SECURITY
    endif
    
    ! unstable conditions
    if (l_local<0._wp) then
      x = (1._wp-15._wp*z_eff/l_local)**0.25_wp
    
      psi_m = 2.0_wp*log((1._wp + x)/2._wp) + log((1._wp + x**2)/2._wp) - 2._wp*atan(x) + M_PI/2._wp
    ! neutral and stable conditions
    else
      psi_m = -4.7_wp*z_eff/l_local
    endif
    
  end function psi_m

end module mo_pbl








