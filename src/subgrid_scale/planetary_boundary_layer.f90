! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module planetary_boundary_layer

  ! In this module, quantities referring to the planetary boundary layer are computed.

  use iso_c_binding
  use constants,   only: EPSILON_SECURITY,M_PI,gravity
  use definitions, only: wp
  use run_nml,     only: dtime
  use diff_nml,    only: h_prandtl
  
  implicit none
  
  real(wp), parameter :: KARMAN = 0.4_wp ! von Karman's constant
  
  contains

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

  function scalar_flux_resistance(roughness_velocity_value,z_agl,roughness_length_value,monin_obukhov_length_value) &
  bind(c,name = "scalar_flux_resistance")

    ! This function returns the surface flux resistance for scalar quantities.

    ! input variables
    real(wp), intent(in) :: roughness_velocity_value,z_agl,roughness_length_value,monin_obukhov_length_value
    ! output variable
    real(wp)             :: scalar_flux_resistance

    ! local variables
    real(wp) :: used_vertical_height
    
    ! height of the prandtl layer
    used_vertical_height = min(z_agl,h_prandtl)
    
    scalar_flux_resistance = 1._wp/(KARMAN*roughness_velocity_value) &
    ! neutral conditions
    *(log(used_vertical_height/roughness_length_value) &
    ! non-neutral conditions
    - psi_h(used_vertical_height,monin_obukhov_length_value) &
    ! interfacial sublayer
    + log(7._wp))

    ! limitting the result for security
    if (scalar_flux_resistance<dtime/z_agl) then
      scalar_flux_resistance = dtime/z_agl
    endif 
    
  end function scalar_flux_resistance

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

    momentum_flux_resistance = 1._wp/(KARMAN*roughness_velocity(wind_h_lowest_layer,z_agl,roughness_length_value)) &
    ! neutral conditions
    *(log(used_vertical_height/roughness_length_value) &
    ! non-neutral conditions
    - psi_m(used_vertical_height,monin_obukhov_length_value))

    ! limitting the result for security
    if (momentum_flux_resistance<dtime/z_agl) then
      momentum_flux_resistance = dtime/z_agl
    endif

  end function momentum_flux_resistance

  function roughness_velocity(wind_speed,z_agl,roughness_length_value) &
  bind(c,name = "roughness_velocity")

    ! This function returns the roughness velocity.

    ! input variables
    real(wp), intent(in) :: wind_speed             ! wind speed at a certain height
    real(wp), intent(in) :: z_agl                  ! height at which the wind speed is valid
    real(wp), intent(in) :: roughness_length_value ! roughness length at this point
    ! output variable
    real(wp)             :: roughness_velocity     ! the result

    ! local variables
    real(wp) :: denominator ! helper variable

    denominator = log(z_agl/roughness_length_value)

    ! security
    if (abs(denominator)<EPSILON_SECURITY) then
      denominator = EPSILON_SECURITY
    endif

    roughness_velocity = wind_speed*KARMAN/denominator

    roughness_velocity = max(EPSILON_SECURITY,roughness_velocity)

  end function roughness_velocity

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
      x = (1._wp - 15._wp*eff/l_local)**0.25_wp
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
      x = (1._wp - 15._wp*eff/l_local)**0.25_wp

      psi_m = 2.0_wp*log((1._wp + x)/2._wp) + log((1._wp + x**2)/2._wp) - 2._wp*atan(x) + M_PI/2._wp
    ! neutral and stable conditions
    else
      psi_m = -4.7_wp*eff/l_local
    endif
    
 end function psi_m
    
end module planetary_boundary_layer
    
    
    
    
    
    
    
    
