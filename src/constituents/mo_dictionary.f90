! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_dictionary

  ! This module contains look-up functions for properties of the atmosphere.

  use mo_definitions, only: wp
  use mo_constants,   only: t_0,n_a,r_v,m_v
  
  implicit none
  
  contains
  
  function molar_fraction_in_dry_air(gas_number)
    
    ! This function returns the molar fraction of certain gases in dry air.
    
    ! gaseous constituents IDs:
    ! 0: dry air
    ! 1: H2O
    ! 2: N2
    ! 3: O2
    ! 4: Ar
    ! 5: CO2
    ! 6: Ne
    ! 7: He
    ! 8: CH4
    ! 9: CO
    ! 10: O3
    ! 11: N2O
    
    integer, intent(in) :: gas_number
    real(wp)            :: molar_fraction_in_dry_air
  
    if (gas_number==2) then
      molar_fraction_in_dry_air = 0.7809_wp
    endif
    if (gas_number==3) then
      molar_fraction_in_dry_air = 0.2095_wp
    endif
    if (gas_number==4) then
      molar_fraction_in_dry_air = 0.0093_wp
    endif
    if (gas_number==5) then
      molar_fraction_in_dry_air = 0.0003_wp
    endif
    if (gas_number==6) then
      molar_fraction_in_dry_air = 1.8e-5_wp
    endif
    if (gas_number==7) then
      molar_fraction_in_dry_air = 5.2e-6_wp
    endif
    if (gas_number==8) then
      molar_fraction_in_dry_air = 1.5e-6_wp
    endif
    if (gas_number==9) then
      molar_fraction_in_dry_air = 1.0e-7_wp
    endif
    if (gas_number==10) then
      molar_fraction_in_dry_air = 1e-6_wp
    endif
    if (gas_number==11) then
    ! https://www.epa.gov/climate-indicators/climate-change-indicators-atmospheric-concentrations-greenhouse-gases
      molar_fraction_in_dry_air = 0.3e-6_wp
    endif
  
  end function molar_fraction_in_dry_air
  
  function calc_o3_vmr(height)

    ! This function calculates the ozone VMR as a function of height.
    ! assumes a Gaussian distribution
    
    real(wp), intent(in) :: height ! height above MSL
    real(wp)             :: calc_o3_vmr
    
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

    ! input arguments
    integer  :: direction
    real(wp) :: temperature

    real(wp) :: phase_trans_heat
    
    ! directions:
    ! 0: gas to liquid
    ! 1: gas to solid
    ! 2: liquid to solid

    phase_trans_heat = 0._wp
    if (direction==0) then
      phase_trans_heat = enthalpy_evaporation(temperature)
    endif
    if (direction==1) then
      phase_trans_heat = enthalpy_sublimation(temperature)
    endif
    if (direction==2) then
      phase_trans_heat = enthalpy_sublimation(temperature) - enthalpy_evaporation(temperature)
    endif
    
  end function phase_trans_heat
  
  function c_p_water(temperature)

    ! This function returns c_p of water.
  
    real(wp), intent(in) :: temperature
    real(wp)             :: c_p_water
  
    ! local variables
    real(wp) :: temp_c
  
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
  
    real(wp), intent(in) :: temperature
    real(wp)             :: c_p_ice
    
    ! local variables
    real(wp) :: temperature_local
    
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
    
    integer, intent(in)  :: const_id
    real(wp), intent(in) :: temperature
    real(wp)             :: c_p_cond
  
    if (mod(const_id-1,2)==0) then
      c_p_cond = c_p_ice(temperature)
    else
      c_p_cond = c_p_water(temperature)
    endif
  
  end function c_p_cond
  
  function enthalpy_evaporation(temperature)
    
    ! This function returns the enthalpy of evaporation depending on the temperature.
    
    real(wp), intent(in) :: temperature
    real(wp)             :: enthalpy_evaporation
    
    ! local variables
    real(wp) :: temperature_local
    
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
    
    real(wp), intent(in) :: temperature
    real(wp)             :: enthalpy_sublimation

    ! local variables
    real(wp) :: temperature_local
    
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
    
    real(wp), intent(in) :: temperature
    real(wp)             :: saturation_pressure_over_water
    
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
  
    ! This function returns the derivative of the saturation pressure in Pa of pure water vapour over plane liquid water
    ! as a function of the temperature in K.
    
    real(wp), intent(in) :: temperature
    real(wp)             :: dsaturation_pressure_over_water_dT
    
    ! local variables
    real(wp) :: temp_c
    
    ! calculating the temperature in degrees Celsius
    temp_c = temperature - t_0
    
    ! these are the limits of this approximation
    if (temp_c>100._wp) then
      temp_c = 100._wp
    endif
    if (temp_c<0._wp) then
      temp_c = 0._wp
    endif
    
    dsaturation_pressure_over_water_dT = saturation_pressure_over_water(temperature) &
    *(4924.99_wp/(temp_c + 237.1_wp)**2 - 1.57_wp/(temp_c + 105._wp))
  
  end function dsaturation_pressure_over_water_dT

  function saturation_pressure_over_ice(temperature)
    
    ! This function returns the saturation pressure in Pa over ice as a function of the temperature in K.
    ! It uses the formula by Huang: A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice, 2018, DOI: 10.1175/JAMC-D-17-0334.1.
    
    real(wp), intent(in) :: temperature
    real(wp)             :: saturation_pressure_over_ice
    
    ! local variables
    real(wp)             :: temp_c

    temp_c = temperature - t_0

    ! clipping too extreme values for this approximation
    if (temp_c<-80._wp) then
      temp_c = -80._wp
    endif
    if (temp_c>0._wp) then
      temp_c = 0._wp
    endif
    
    saturation_pressure_over_ice = exp(43.494_wp-6545.8_wp/(temp_c+278._wp))/(temp_c+868._wp)**2
    
  end function saturation_pressure_over_ice
  
  function dsaturation_pressure_over_ice_dT(temperature)

    ! This function returns derivative of the the saturation pressure in Pa of pure water vapour over plane ice
    ! as a function of the temperature in K.
    
    real(wp), intent(in) :: temperature
    real(wp)             :: dsaturation_pressure_over_ice_dT
    
    ! local variables
    real(wp) :: temp_c
    
    ! calculating the temperature in degrees Celsius
    temp_c = temperature - t_0
    
    ! this is the stability limit
    if (temp_c<-80._wp) then
      temp_c = -80._wp
    endif
    ! at temperatures > 0 degrees Celsius ice cannot exist in equilibrium which is why this is clipped
    if (temp_c>0._wp) then
      temp_c = 0._wp
    endif
    
    dsaturation_pressure_over_ice_dT = saturation_pressure_over_ice(temperature) &
    *(6545.8_wp/(temp_c + 278._wp)**2 - 2._wp/(temp_c + 868._wp))

  end function dsaturation_pressure_over_ice_dT
  
  function enhancement_factor_over_water(air_pressure)

    ! This function calculates the enhancement factor over water, which accounts for the fact the the saturation vapour pressure is different in moist air compared to pure water vapour.
    ! It uses the formula by Huang: A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice, 2018, DOI: 10.1175/JAMC-D-17-0334.1.

    real(wp), intent(in) :: air_pressure
    real(wp)             :: enhancement_factor_over_water
    
    enhancement_factor_over_water = 0.99882_wp*exp(0.00000008_wp*air_pressure)

    end function enhancement_factor_over_water
  
  function enhancement_factor_over_ice(air_pressure)

    ! This function calculates the enhancement factor over ice, which accounts for the fact the the saturation vapour pressure is different in moist air compared to pure water vapour.
    ! It uses the formula by Huang: A Simple Accurate Formula for Calculating Saturation Vapor Pressure of Water and Ice, 2018, DOI: 10.1175/JAMC-D-17-0334.1.

    real(wp), intent(in) :: air_pressure
    real(wp)             :: enhancement_factor_over_ice
    
    enhancement_factor_over_ice = 0.99882_wp*exp(0.00000008_wp*air_pressure)

  end function enhancement_factor_over_ice

end module mo_dictionary











