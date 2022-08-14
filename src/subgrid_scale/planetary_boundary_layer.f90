! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module planetary_boundary_layer

  ! In this module, quantities referring to the planetary boundary layer are computed.

  use iso_c_binding
  use constants,   only: EPSILON_SECURITY,M_PI,gravity
  use definitions, only: wp
  use run_nml,     only: dtime
  use diff_nml,    only: h_prandtl
  use grid_nml,    only: n_scalars,n_scalars_h,n_vectors_per_layer,n_layers,n_vectors
  use surface_nml, only: lprog_soil_temp
  
  implicit none
  
  real(wp), parameter :: KARMAN = 0.4_wp ! von Karman's constant
  
  contains
  
  subroutine update_sfc_turb_quantities(is_land,roughness_length,monin_obukhov_length,z_scalar,z_vector, &
                                        theta_v_bg,theta_v_pert,v_squared,roughness_velocity,scalar_flux_resistance) &
  bind(c,name = "update_sfc_turb_quantities")
  
    ! This subroutine updates surface-related turbulence quantities.
    
    integer,  intent(in)    :: is_land(n_scalars_h)
    real(wp), intent(in)    :: z_scalar(n_scalars),z_vector(n_vectors),theta_v_bg(n_scalars), &
                               theta_v_pert(n_scalars),v_squared(n_scalars)
    real(wp), intent(out)   :: monin_obukhov_length(n_scalars_h),roughness_velocity(n_scalars_h), &
                               scalar_flux_resistance(n_scalars_h)
    real(wp), intent(inout) :: roughness_length(n_scalars_h)
    
    ! local variables
    integer  :: ji
    real(wp) :: u_lowest_layer,u10,z_agl,theta_v_lowest_layer,theta_v_second_layer, &
                dz,dtheta_v_dz,w_pert,theta_v_pert_value,w_pert_theta_v_pert_avg,w_theta_v_corr
    
    ! semi-empirical coefficient
    w_theta_v_corr = 0.2_wp
    
    !$omp parallel do private(ji,u_lowest_layer,u10,z_agl,theta_v_lowest_layer,theta_v_second_layer, &
    !$omp dz,dtheta_v_dz,w_pert,theta_v_pert_value,w_pert_theta_v_pert_avg)
    do ji=1,n_scalars_h
      z_agl = z_scalar(n_scalars-n_scalars_h+ji)-z_vector(n_vectors-n_scalars_h+ji)
      
      ! wind speed in the lowest layer
      u_lowest_layer = v_squared(n_scalars-n_scalars_h+ji)**0.5_wp
        
      ! calculating the 10 m wind velocity from the logarithmic wind profile
      u10 = u_lowest_layer*log(10._wp/roughness_length(ji))/log(z_agl/roughness_length(ji))
      
      ! only over the sea the roughness length is time-dependant (because of the waves)
      if (is_land(ji)==0) then
        ! calculating the roughness length fom the wind velocity
        roughness_length(ji) = roughness_length_from_u10_sea(u10)
      endif
      
      ! updating the roughness velocity
      roughness_velocity(ji) = calc_roughness_velocity(u_lowest_layer,z_agl,roughness_length(ji))
      
      ! theta_v in the lowest layer
      theta_v_lowest_layer = theta_v_bg(n_scalars-n_scalars_h+ji) + theta_v_pert(n_scalars-n_scalars_h+ji)
      ! theta_v in the second-lowest layer
      theta_v_second_layer = theta_v_bg(n_scalars-2*n_scalars_h+ji) + theta_v_pert(n_scalars-2*n_scalars_h+ji)
      
      ! delta z
      dz = z_scalar(n_scalars-2*n_scalars_h+ji)-z_scalar(n_scalars-n_scalars_h+ji)
      
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
      monin_obukhov_length(ji) = -theta_v_lowest_layer*roughness_velocity(ji)**3/(KARMAN*gravity*w_pert_theta_v_pert_avg)
    enddo
    !$omp end parallel do
    
    ! updating the surface flux resistance acting on scalar quantities (moisture and sensible heat)
    if (lprog_soil_temp) then
      !$omp parallel do private(ji)
      do ji=1,n_scalars_h
        scalar_flux_resistance(ji) = calc_scalar_flux_resistance(roughness_velocity(ji), &
                                            z_scalar(n_scalars-n_scalars_h+ji)-z_vector(n_layers*n_vectors_per_layer+ji), &
                                            roughness_length(ji),monin_obukhov_length(ji))
      enddo
      !$omp end parallel do
    endif
    
  end subroutine update_sfc_turb_quantities

  function roughness_length_from_u10_sea(u10) &
  bind(c,name = "roughness_length_from_u10_sea")
  
    ! This function returns the roughness length as a function of the mean wind speed at 10 m above a fully developed sea.

    ! input variable
    real(wp), intent(in) :: u10
    ! output variable
    real(wp)             :: roughness_length_from_u10_sea

    ! local variables
    real(wp) :: swh,period,wavelength ! properties of the wave field

    ! refer to Stensrud,Parameterization schemes (2007), p.130

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

  function calc_scalar_flux_resistance(roughness_velocity_value,z_agl,roughness_length_value,monin_obukhov_length_value) &
  bind(c,name = "calc_scalar_flux_resistance")

    ! This function returns the surface flux resistance for scalar quantities.

    ! input variables
    real(wp), intent(in) :: roughness_velocity_value,z_agl,roughness_length_value,monin_obukhov_length_value
    ! output variable
    real(wp)             :: calc_scalar_flux_resistance

    ! local variables
    real(wp) :: used_vertical_height
    
    ! height of the prandtl layer
    used_vertical_height = min(z_agl,h_prandtl)
    
    calc_scalar_flux_resistance = 1._wp/(KARMAN*roughness_velocity_value) &
    ! neutral conditions
    *(log(used_vertical_height/roughness_length_value) &
    ! non-neutral conditions
   -psi_h(used_vertical_height,monin_obukhov_length_value) &
    ! interfacial sublayer
    + log(7._wp))

    ! limitting the result for security
    if (calc_scalar_flux_resistance<dtime/z_agl) then
      calc_scalar_flux_resistance = dtime/z_agl
    endif 
    
  end function calc_scalar_flux_resistance

  function momentum_flux_resistance(wind_h_lowest_layer,z_agl,roughness_length_value,monin_obukhov_length_value) &
  bind(c,name = "momentum_flux_resistance")

    ! This function returns the surface flux resistance for momentum.

    ! input variables
    real(wp), intent(in) :: wind_h_lowest_layer,z_agl,roughness_length_value,monin_obukhov_length_value
    ! output variable
    real(wp)             :: momentum_flux_resistance

    ! local variables
    real(wp) :: used_vertical_height

    ! height of the prandtl layer
    used_vertical_height = min(z_agl,h_prandtl)

    momentum_flux_resistance = 1._wp/(KARMAN*calc_roughness_velocity(wind_h_lowest_layer,z_agl,roughness_length_value)) &
    ! neutral conditions
    *(log(used_vertical_height/roughness_length_value) &
    ! non-neutral conditions
   -psi_m(used_vertical_height,monin_obukhov_length_value))

    ! limitting the result for security
    if (momentum_flux_resistance<dtime/z_agl) then
      momentum_flux_resistance = dtime/z_agl
    endif

  end function momentum_flux_resistance

  function calc_roughness_velocity(wind_speed,z_agl,roughness_length_value) &
  bind(c,name = "calc_roughness_velocity")

    ! This function returns the roughness velocity.

    ! input variables
    real(wp), intent(in) :: wind_speed              ! wind speed at a certain height
    real(wp), intent(in) :: z_agl                   ! height at which the wind speed is valid
    real(wp), intent(in) :: roughness_length_value  ! roughness length at this point
    ! output variable
    real(wp)             :: calc_roughness_velocity ! the result

    ! local variables
    real(wp) :: denominator ! helper variable

    denominator = log(z_agl/roughness_length_value)

    ! security
    if (abs(denominator)<EPSILON_SECURITY) then
      denominator = EPSILON_SECURITY
    endif

    calc_roughness_velocity = wind_speed*KARMAN/denominator

    calc_roughness_velocity = max(EPSILON_SECURITY,calc_roughness_velocity)

  end function calc_roughness_velocity

  function psi_h(eff,l)

    ! This is a helper function for the correction to the surface scalar flux resistance for non-neutral conditions.

    ! input variables
    real(wp), intent(in) :: eff   ! effective height above the surface
    real(wp), intent(in) :: l     ! Monin-Obukhov length
    ! output variable
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
      x = (1._wp-15._wp*eff/l_local)**0.25_wp
      psi_h = 2._wp*log((1._wp + x**2)/2._wp)     
    ! neutral and stable conditions
    else
      psi_h = -4._wp*eff/l_local
    endif
    
  end function psi_h

  function psi_m(eff,l)

    ! This is a helper function for the correction to the surface momentum flux resistance for non-neutral conditions.

    ! input variables
    real(wp), intent(in) :: eff   ! effective height above the surface
    real(wp), intent(in) :: l     ! Monin-Obukhov length
    ! output variable
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
      x = (1._wp-15._wp*eff/l_local)**0.25_wp

      psi_m = 2.0_wp*log((1._wp + x)/2._wp) + log((1._wp + x**2)/2._wp)-2._wp*atan(x) + M_PI/2._wp
    ! neutral and stable conditions
    else
      psi_m = -4.7_wp*eff/l_local
    endif
    
 end function psi_m
    
end module planetary_boundary_layer
    
    
    
    
    
    
    
    
