! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_dictionary

  ! This module contains look-up functions for properties of the atmosphere.

  use mo_definitions, only: wp
  use mo_constants,   only: t_0,m_v,rho_h2o,M_PI
  
  implicit none
  
  contains
  
  function molar_fraction_in_dry_air(gas_number)
    
    ! This function returns the molar fraction of certain gases in dry air.
    
    integer, intent(in) :: gas_number                ! defines the gas (definition see below)
    real(wp)            :: molar_fraction_in_dry_air ! result
    
    ! N2
    if (gas_number==2) then
      molar_fraction_in_dry_air = 0.7809_wp
    endif
    ! O2
    if (gas_number==3) then
      molar_fraction_in_dry_air = 0.2095_wp
    endif
    ! Ar
    if (gas_number==4) then
      molar_fraction_in_dry_air = 0.0093_wp
    endif
    ! CO2
    if (gas_number==5) then
      molar_fraction_in_dry_air = 0.0003_wp
    endif
    ! Ne
    if (gas_number==6) then
      molar_fraction_in_dry_air = 1.8e-5_wp
    endif
    ! He
    if (gas_number==7) then
      molar_fraction_in_dry_air = 5.2e-6_wp
    endif
    ! CH4
    if (gas_number==8) then
      molar_fraction_in_dry_air = 1.5e-6_wp
    endif
    ! CO
    if (gas_number==9) then
      molar_fraction_in_dry_air = 1.0e-7_wp
    endif
    ! O3
    if (gas_number==10) then
      molar_fraction_in_dry_air = 1e-6_wp
    endif
    ! N2O
    if (gas_number==11) then
      ! source: https://www.epa.gov/climate-indicators/climate-change-indicators-atmospheric-concentrations-greenhouse-gases
      molar_fraction_in_dry_air = 0.3e-6_wp
    endif
    
  end function molar_fraction_in_dry_air
  
  function calc_o3_vmr(height)
    
    ! This function calculates the ozone VMR as a function of height.
    ! assumes a Gaussian distribution
    
    real(wp), intent(in) :: height      ! height above MSL
    real(wp)             :: calc_o3_vmr ! result
    
    ! local variables
    real(wp) :: fwhm = 20e3_wp       ! full width at half maximum
    real(wp) :: max_height = 34e3_wp ! height of the maximum of the distribution
    real(wp) :: max_vmr = 8.5e-6_wp  ! maximum volume mixing ratio
    real(wp) :: sigma                ! standard deviation
    real(wp) :: distance             ! distance from the maximum
    
    ! calculation of the result
    sigma = fwhm/(8._wp*log(2._wp))**0.5_wp
    distance = height-max_height
    calc_o3_vmr = max_vmr*exp(-distance**2/(2._wp*sigma**2))
    
  end function calc_o3_vmr

  function phase_trans_heat(direction,temperature)
    
    ! This function calculates the phase transition heat.
    
    integer  :: direction        ! defines the kind of phase transition (definition see below)
    real(wp) :: temperature      ! temperature
    real(wp) :: phase_trans_heat ! result

    phase_trans_heat = 0._wp
    ! 1: gas to liquid
    if (direction==1) then
      phase_trans_heat = enthalpy_evaporation(temperature)
    endif
    ! 2: gas to solid
    if (direction==2) then
      phase_trans_heat = enthalpy_sublimation(temperature)
    endif
    ! 3: liquid to solid
    if (direction==3) then
      phase_trans_heat = enthalpy_sublimation(temperature) - enthalpy_evaporation(temperature)
    endif
    
  end function phase_trans_heat
  
  function c_p_water(temperature)
    
    ! This function returns c_p of water.
  
    real(wp), intent(in) :: temperature ! temperature
    real(wp)             :: c_p_water   ! result
  
    ! local variables
    real(wp) :: temp_c ! temperature in degrees Celsius
  
    ! calculating the temperature in degrees Celsius
    temp_c = temperature - t_0
  
    ! For "positive" temperatures we use the formula cited in Pruppacher and Klett (2010), p. 93, Eq. (3-15).
  
    if (temp_c>=0._wp) then
      ! clipping values that are too extreme for this approximation
      if (temp_c>35._wp) then
        temp_c = 35._wp
      endif
      c_p_water = 0.9979_wp + 3.1e-6_wp*(temp_c-35._wp)**2 + 3.8e-9_wp*(temp_c-35._wp)**4
    ! This is the case of supercooled water. We use the formula cited in Pruppacher and Klett (2010), p. 93, Eq. (3-16).
    else
      ! clipping values that are too extreme for this approximation
      if (temp_c<-37._wp) then
        temp_c = -37._wp
      endif
      c_p_water = 1.000938_wp - 2.7052e-3_wp*temp_c - 2.3235e-5_wp*temp_c**2 + 4.3778e-6_wp*temp_c**3 + 2.7136e-7_wp*temp_c**4
    endif
    
    ! unit conversion from IT cal/(g*K) to J/(kg*K)
    c_p_water = 4186.8_wp*c_p_water
  
  end function c_p_water
  
  function c_p_ice(temperature)
    
    ! This function returns c_p of ice.
    ! It follows Eq. (4) in Murphy DM, Koop T. Review of the vapour pressures of ice and supercooled water for atmospheric applications.
    ! QUARTERLY JOURNAL OF THE ROYAL METEOROLOGICAL SOCIETY. 2005;131(608):1539-1565.
  
    real(wp), intent(in) :: temperature ! temperature
    real(wp)             :: c_p_ice     ! result
    
    ! local variables
    real(wp) :: temperature_local ! local copy of the temperature
    
    temperature_local = temperature
  
    ! ice cannot exist in equilibrium at temperatures > T_0
    if (temperature_local>t_0) then
      temperature_local = t_0
    endif
    ! clipping values that are too extreme for this approximation
    if (temperature_local<20._wp) then
      temperature_local = 20._wp
    endif
    c_p_ice = -2.0572_wp + 0.14644_wp*temperature_local + 0.06163_wp*temperature_local*exp(-(temperature_local/125.1_wp)**2)
    ! unit conversion from J/(mol*K) to J/(kg*K)
    c_p_ice = c_p_ice/m_v
    
  end function c_p_ice

  function c_p_cond(const_id,temperature)
  
    ! This function resturns c_p of a specific condensed constituent.
    
    integer,  intent(in) :: const_id    ! index of the constituent
    real(wp), intent(in) :: temperature ! temperature
    real(wp)             :: c_p_cond    ! result
  
    if (mod(const_id-1,2)==0) then
      c_p_cond = c_p_ice(temperature)
    else
      c_p_cond = c_p_water(temperature)
    endif
  
  end function c_p_cond
  
  function enthalpy_evaporation(temperature)
    
    ! This function returns the enthalpy of evaporation depending on the temperature.
    
    real(wp), intent(in) :: temperature          ! temperature
    real(wp)             :: enthalpy_evaporation ! result
    
    ! local variables
    real(wp) :: temperature_local ! local copy of the temperature
    
    temperature_local = temperature
    
    if (temperature_local<t_0) then
      ! This follows Eq. (9) in Murphy DM, Koop T. Review of the vapour pressures of ice and supercooled water for atmospheric applications.
      ! QUARTERLY JOURNAL OF THE ROYAL METEOROLOGICAL SOCIETY. 2005;131(608):1539-1565.
      ! clipping values that are too extreme for these approximations
      if (temperature_local<30._wp) then
        temperature_local = 30._wp
      endif
      enthalpy_evaporation = 56579._wp - 42.212_wp*temperature_local + exp(0.1149_wp*(281.6_wp - temperature_local))
      ! unit conversion from J/mol to J/kg
      enthalpy_evaporation = enthalpy_evaporation/m_v
    else
      ! This follows the formula (Eq. (8)) cited by Huang:
      ! A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice, 2018, DOI: 10.1175/JAMC-D-17-0334.1.
      ! clipping values that are too extreme for these approximations
      if (temperature_local>t_0+100._wp) then
        temperature_local = t_0 + 100._wp
      endif
      enthalpy_evaporation = 3151378._wp - 2386._wp*temperature_local
    endif

  end function enthalpy_evaporation
  
  function enthalpy_sublimation(temperature)
  
    ! This function returns the enthalpy of sublimation depending on the temperature.
    ! It follows Eq. (5) in Murphy DM, Koop T. Review of the vapour pressures of ice and supercooled water for atmospheric applications.
    ! QUARTERLY JOURNAL OF THE ROYAL METEOROLOGICAL SOCIETY. 2005;131(608):1539-1565.
    
    real(wp), intent(in) :: temperature          ! temperature
    real(wp)             :: enthalpy_sublimation ! result

    ! local variables
    real(wp) :: temperature_local ! local copy of the temperature
    
    temperature_local = temperature
     
    ! clipping values that are too extreme for this approximation
    if (temperature_local<30._wp) then
      temperature_local = 30._wp
    endif
    ! sublimation is not happening in thermodynamic equilibrium at temperatures > t_0
    if (temperature_local>t_0) then
      temperature_local = t_0
    endif

    enthalpy_sublimation = 46782.5_wp + 35.8925_wp*temperature_local - 0.07414_wp*temperature_local**2 &
    + 541.5_wp*exp(-(temperature_local/123.75_wp)**2)

    ! unit conversion from J/mol to J/kg
    enthalpy_sublimation = enthalpy_sublimation/m_v
    
  end function enthalpy_sublimation
  
  function saturation_pressure_over_water(temperature)

    ! This function returns the saturation pressure in Pa over liquid water as a function of the temperature in K.
    ! It uses the formula by Huang: A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice, 2018, DOI: 10.1175/JAMC-D-17-0334.1.
    
    real(wp), intent(in) :: temperature                    ! temperature
    real(wp)             :: saturation_pressure_over_water ! result
    
    ! local variables
    real(wp)  :: temp_c ! temperature in degrees Celsius

    temp_c = temperature - t_0
    ! clipping too extreme values for this approximation
    if (temp_c>100._wp) then
      temp_c = 100._wp
    endif
    
    if (temp_c>0._wp) then
      saturation_pressure_over_water = exp(34.494_wp-4924.99_wp/(temp_c+237.1_wp))/(temp_c+105._wp)**1.57_wp
    ! For super-cooled water we use the formula cited in Pruppacher and Klett (2010), p. 854, Eq. (A.4-1).
    else
      ! Clipping values that are too extreme for this approximation.
      if (temp_c<-50._wp) then
        temp_c = -50._wp
      endif
      saturation_pressure_over_water &
      = 6.107799961_wp &
      + 4.436518521e-1_wp*temp_c &
      + 1.428945805e-2_wp*temp_c**2 &
      + 2.650648471e-4_wp*temp_c**3 &
      + 3.031240396e-6_wp*temp_c**4 &
      + 2.034080948e-8_wp*temp_c**5 &
      + 6.136820929e-11_wp*temp_c**6
    endif

  end function saturation_pressure_over_water
  
  function dsaturation_pressure_over_water_dT(temperature)
    
    ! This function computes the derivative of the function saturation_pressure_over_water.
    
    real(wp), intent(in) :: temperature                        ! temperature
    real(wp)             :: dsaturation_pressure_over_water_dT ! result
    
    ! local variables
    real(wp) :: temp_c ! temperature in degrees Celsius
    
    ! calculating the temperature in degrees Celsius
    temp_c = temperature - t_0
    
    ! these are the limits of this approximation
    if (temp_c>100._wp) then
      dsaturation_pressure_over_water_dT = 0._wp
    elseif (temp_c<0._wp) then
      dsaturation_pressure_over_water_dT = 0._wp
    else
      dsaturation_pressure_over_water_dT = saturation_pressure_over_water(temperature) &
      *(4924.99_wp/(temp_c + 237.1_wp)**2 - 1.57_wp/(temp_c + 105._wp))
    endif
  
  end function dsaturation_pressure_over_water_dT

  function saturation_pressure_over_ice(temperature)
    
    ! This function returns the saturation pressure in Pa over ice as a function of the temperature in K.
    ! It blends the two formulas of Huang and Murphy.
    
    real(wp), intent(in) :: temperature                  ! temperature
    real(wp)             :: saturation_pressure_over_ice ! result
    
    ! local variables
    real(wp) :: t_local      ! local copy of temperature
    real(wp) :: temp_c       ! temperature in degrees Celsius
    real(wp) :: huang_weight ! weight of the Huang formula

    t_local = temperature

    temp_c = t_local - t_0
    
    if (temp_c>=-80._wp) then
      ! at temperatures > 0 degrees Celsius ice cannot exist in equilibrium which is why this is clipped
      if (t_local>t_0) then
        t_local = t_0
      endif
      saturation_pressure_over_ice = saturation_pressure_ice_huang(t_local)
    elseif (temp_c>=-100._wp) then
      huang_weight = (temp_c+100._wp)/20._wp
      saturation_pressure_over_ice = huang_weight*saturation_pressure_ice_huang(t_local) &
                                     + (1._wp-huang_weight)+saturation_pressure_ice_murphy(t_local)
    else
      ! clipping too extreme values for this approximation
      if (t_local<110._wp) then
        t_local = 110._wp
      endif
      saturation_pressure_over_ice = saturation_pressure_ice_murphy(t_local)
    endif
    
  end function saturation_pressure_over_ice
  
  function dsaturation_pressure_over_ice_dT(temperature)

    ! This function computes the derivative of the function saturation_pressure_over_ice.
    
    real(wp), intent(in) :: temperature                      ! temperature
    real(wp)             :: dsaturation_pressure_over_ice_dT ! result
    
    ! local variables
    real(wp) :: t_local      ! local copy of temperature
    real(wp) :: temp_c       ! temperature in degrees Celsius
    real(wp) :: huang_weight ! weight of the Huang formula
    
    t_local = temperature
    
    ! calculating the temperature in degrees Celsius
    temp_c = t_local - t_0
    
    ! at temperatures > 0 degrees Celsius ice cannot exist in equilibrium which is why this is clipped
    if (temp_c>0._wp) then
       dsaturation_pressure_over_ice_dT = 0._wp
    elseif (temp_c>=-80._wp) then
      dsaturation_pressure_over_ice_dT = dsaturation_pressure_ice_huang_dT(t_local)
    elseif (temp_c>=-100._wp) then
      huang_weight = (temp_c+100._wp)/20._wp
      dsaturation_pressure_over_ice_dT = huang_weight*dsaturation_pressure_ice_huang_dT(t_local) &
                                         + (1._wp-huang_weight)+dsaturation_pressure_ice_murphy_dT(t_local)
    elseif (t_local>=110._wp) then
      dsaturation_pressure_over_ice_dT = dsaturation_pressure_ice_murphy_dT(t_local)
    ! clipping too extreme values for this approximation
    else
      dsaturation_pressure_over_ice_dT = 0._wp
    endif

  end function dsaturation_pressure_over_ice_dT
  
  function enhancement_factor_over_water(air_pressure)

    ! This function calculates the enhancement factor over water, which accounts for the fact the the saturation vapour pressure is different in moist air compared to pure water vapour.
    ! It uses the formula by Huang: A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice, 2018, DOI: 10.1175/JAMC-D-17-0334.1.

    real(wp), intent(in) :: air_pressure                  ! air pressure
    real(wp)             :: enhancement_factor_over_water ! result
    
    enhancement_factor_over_water = 1.00071_wp*exp(0.000000045_wp*air_pressure)

    end function enhancement_factor_over_water
  
  function enhancement_factor_over_ice(air_pressure)

    ! This function calculates the enhancement factor over ice, which accounts for the fact the the saturation vapour pressure is different in moist air compared to pure water vapour.
    ! It uses the formula by Huang: A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice, 2018, DOI: 10.1175/JAMC-D-17-0334.1.

    real(wp), intent(in) :: air_pressure                ! air pressure
    real(wp)             :: enhancement_factor_over_ice ! result
    
    enhancement_factor_over_ice = 0.99882_wp*exp(0.00000008_wp*air_pressure)

  end function enhancement_factor_over_ice
  
  function saturation_pressure_ice_huang(temperature)
  
    ! This function computes the saturation pressure over ice.
    ! It follows the formula by Huang: A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice, 2018, DOI: 10.1175/JAMC-D-17-0334.1.
    
    real(wp), intent(in) :: temperature                   ! temperature
    real(wp)             :: saturation_pressure_ice_huang ! result
    
    ! local variables
    real(wp) :: temp_c  ! temperature in degrees Celsius
    
    temp_c = temperature - t_0
    
    saturation_pressure_ice_huang = exp(43.494_wp-6545.8_wp/(temp_c+278._wp))/(temp_c+868._wp)**2
  
  end function saturation_pressure_ice_huang
  
  function saturation_pressure_ice_murphy(temperature)
  
    ! This function computes the saturation pressure over ice.
    ! It follows Eq. (7) in Murphy DM, Koop T. Review of the vapour pressures of ice and supercooled water for atmospheric applications.
    ! QUARTERLY JOURNAL OF THE ROYAL METEOROLOGICAL SOCIETY. 2005;131(608):1539-1565.

    real(wp), intent(in) :: temperature                    ! temperature
    real(wp)             :: saturation_pressure_ice_murphy ! result
    
    ! computing the result
    saturation_pressure_ice_murphy = exp(9.550426_wp-5723.265_wp/temperature+3.53068_wp*log(temperature)-0.00728332_wp*temperature)
  
  end function saturation_pressure_ice_murphy
  
  function dsaturation_pressure_ice_huang_dT(temperature)
  
    ! This function computes the derivative of the function saturation_pressure_ice_huang.
  
    real(wp), intent(in) :: temperature                       ! temperature
    real(wp)             :: dsaturation_pressure_ice_huang_dT ! result
    
    ! local variables
    real(wp) :: temp_c  ! temperature in degrees Celsius
    
    temp_c = temperature - t_0
    
    dsaturation_pressure_ice_huang_dT = saturation_pressure_ice_huang(temperature) &
    *(6545.8_wp/(temp_c + 278._wp)**2 - 2._wp/(temp_c + 868._wp))
  
  end function dsaturation_pressure_ice_huang_dT
  
  function dsaturation_pressure_ice_murphy_dT(temperature)
  
    ! This function computes the derivative of the function saturation_pressure_ice_murphy.
  
    real(wp), intent(in) :: temperature                        ! temperature
    real(wp)             :: dsaturation_pressure_ice_murphy_dT ! result
  
    ! computing the result
    dsaturation_pressure_ice_murphy_dT = saturation_pressure_ice_murphy(temperature) &
    *(5723.265_wp/temperature**2+3.53068_wp/temperature-0.00728332_wp)
      
  end function dsaturation_pressure_ice_murphy_dT
  
  function cloud_droplets_radius()
  
    ! This function returns the mean radius of cloud droplets.
  
    real(wp) :: cloud_droplets_radius ! result
    
    cloud_droplets_radius = 12.e-6_wp
  
  end function cloud_droplets_radius
  
  function rain_drops_radius(liquid_density)
    
    ! This function returns the mean radius of rain drops assuming a Marshall-Palmer (1948) distribution.
    
    real(wp) :: liquid_density    ! the density of liquid water (clouds + rain)
    real(wp) :: rain_drops_radius ! result
    
    ! local variables
    real(wp) :: n_0 = 8.e6_wp ! prefactor of the Marshall-Palmer distribution
    
    rain_drops_radius = 0.5_wp*(liquid_density/(rho_h2o*M_PI*n_0))**0.25_wp
    
  end function rain_drops_radius
  
  function ice_particles_radius()
  
    ! This function returns the mean radius of ice particles.
  
    real(wp) :: ice_particles_radius ! result
    
    ice_particles_radius = 95.e-6_wp
  
  end function ice_particles_radius
  
  function snow_particles_radius()
  
    ! This function returns the mean radius of snowflakes.
  
    real(wp) :: snow_particles_radius ! result
    
    snow_particles_radius = 500.e-6_wp
  
  end function snow_particles_radius

end module mo_dictionary











