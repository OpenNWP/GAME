! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module phase_trans

  ! This file contains functions calculating everything related to phase transition rates.

  use iso_c_binding
  use mo_definitions,     only: wp
  use run_nml,            only: dtime
  use grid_nml,           only: n_scalars,n_scalars_h,n_layers
  use mo_constants,       only: r_v,t_0,r_d
  use constituents_nml,   only: n_condensed_constituents,n_constituents
  use dictionary,         only: saturation_pressure_over_water,saturation_pressure_over_ice, &
                                dsaturation_pressure_over_water_dT,dsaturation_pressure_over_ice_dT, &
                                phase_trans_heat,enhancement_factor_over_water,enhancement_factor_over_ice
  use derived_quantities, only: c_v_mass_weighted_air
  use surface_nml,        only: nsoillays,lsfc_phase_trans
  use manage_radiation_calls, only: call_radiation

  implicit none
  
  contains

  subroutine calc_h2otracers_source_rates(rho,temperature,layer_thickness,temperature_soil, &
                                          phase_trans_rates,phase_trans_heating_rate, &
                                          scalar_flux_resistance,is_land,power_flux_density_latent) &
  bind(c,name = "calc_h2otracers_source_rates")
    
    ! This subroutine calculates phase transition rates and associated heat source rates.
    ! It assumes the following order for the constituents:
    ! precipitating ice - precipitating liquid water - cloud ice - liquid cloud water - moist air - water vapour
    
    real(wp), intent(in)  :: rho(n_constituents*n_scalars),temperature(n_scalars),layer_thickness(n_scalars), &
                             temperature_soil(nsoillays*n_scalars_h),scalar_flux_resistance(n_scalars_h)
    real(wp), intent(out) :: phase_trans_rates((n_condensed_constituents+1)*n_scalars), &
                             phase_trans_heating_rate(n_scalars),power_flux_density_latent(n_scalars_h)
    integer,  intent(in)  :: is_land(n_scalars_h)
    
    ! local variables
    integer  :: ji,layer_index,h_index
    real(wp) :: diff_density,phase_trans_density,saturation_pressure,water_vapour_pressure, &
                diff_density_sfc,saturation_pressure_sfc,dry_pressure,air_pressure, &
                a,b,c,p,q,enhancement_factor,maximum_cloud_water_content
    
    ! maximum cloud water content in (kg cloud)/(kg dry air).
    maximum_cloud_water_content = 0.2e-3_wp
    
    ! loop over all grid boxes
    !$omp parallel do private(ji,diff_density,phase_trans_density,saturation_pressure,water_vapour_pressure, &
    !$omp layer_index,h_index,diff_density_sfc,saturation_pressure_sfc,dry_pressure,air_pressure, &
    !$omp a,b,c,p,q,enhancement_factor)
    do ji=1,n_scalars
      ! Preparation
      ! -----------
      layer_index = (ji-1)/n_scalars_h
      
      ! determining the saturation pressure
      ! "positive" temperatures (the saturation pressure is different over water compared to over ice)
      if (temperature(ji)>=t_0) then
        saturation_pressure = saturation_pressure_over_water(temperature(ji))
      ! "negative" temperatures
      else
        saturation_pressure = saturation_pressure_over_ice(temperature(ji))
      endif
      
      ! determining the water vapour pressure (using the EOS)
      water_vapour_pressure = rho((n_condensed_constituents+1)*n_scalars+ji)*r_v*temperature(ji)
      
      ! determining the water vapour pressure (using the EOS)
      dry_pressure = (rho(n_condensed_constituents*n_scalars+ji) - rho((n_condensed_constituents+1)*n_scalars+ji)) &
      *r_d*temperature(ji)
        
      ! calculating the total air pressure
      air_pressure = dry_pressure + water_vapour_pressure
        
      ! multiplying the saturation pressure by the enhancement factor
      if (temperature(ji)>=t_0) then
        enhancement_factor = enhancement_factor_over_water(air_pressure)
      ! "negative" temperatures
      else
        enhancement_factor = enhancement_factor_over_ice(air_pressure)
      endif
      
      saturation_pressure = enhancement_factor*saturation_pressure
        
      ! Clouds
      ! ------
      ! the case where the air is not over-saturated
      if (saturation_pressure>=water_vapour_pressure) then
        ! temperature>=0° C
        if (temperature(ji)>=t_0) then
          ! It is assumed that the still present ice vanishes within one time step.
          phase_trans_rates(2*n_scalars+ji) = -rho(2*n_scalars+ji)/dtime
                
          ! The amount of liquid water per volume that will evaporate.
          ! In case the air cannot take all the water, not everything will evaporate.
          a = -r_v*phase_trans_heat(0,temperature(ji))/c_v_mass_weighted_air(rho,temperature,ji-1)
          b = r_v*temperature(ji) - r_v*rho((n_condensed_constituents+1)*n_scalars+ji) &
          *phase_trans_heat(0,temperature(ji))/c_v_mass_weighted_air(rho,temperature,ji-1) &
          + enhancement_factor*dsaturation_pressure_over_water_dT(temperature(ji)) &
          *phase_trans_heat(0,temperature(ji))/c_v_mass_weighted_air(rho,temperature,ji-1)
          c = water_vapour_pressure - saturation_pressure
          p = b/a
          q = c/a
          diff_density = -0.5_wp*p - (0.25_wp*p**2 - q)**0.5_wp
          phase_trans_density = min(rho(3*n_scalars+ji), diff_density)
                
          ! the tendency for the water vapour
          phase_trans_rates(4*n_scalars+ji) = phase_trans_density/dtime
                
          ! The source rate for the liquid water consists of two terms:
          ! 1.) the melting
          ! 2.) the evaporation
                
          phase_trans_rates(3*n_scalars+ji) = rho(2*n_scalars+ji)/dtime - phase_trans_density/dtime
                
          ! the heat source rates
          phase_trans_heating_rate(ji) &
          ! melting
          = phase_trans_rates(2*n_scalars+ji)*phase_trans_heat(2,temperature(ji)) &
          ! evaporation
          - phase_trans_density*phase_trans_heat(0,temperature(ji))/dtime
        ! temperature<0° C
        else
          ! It is assumed that the still present liquid water vanishes within one time step.
          phase_trans_rates(3*n_scalars+ji) = -rho(3*n_scalars+ji)/dtime
                
          ! The amount of ice per volume that will sublimate.
          ! In case the air cannot take all the water, not everything will sublimate.
                
          a = -r_v*phase_trans_heat(1,temperature(ji))/c_v_mass_weighted_air(rho,temperature,ji-1)
          b = r_v*temperature(ji) - r_v*rho((n_condensed_constituents+1)*n_scalars+ji) &
          *phase_trans_heat(1,temperature(ji))/c_v_mass_weighted_air(rho,temperature,ji-1) &
          + enhancement_factor*dsaturation_pressure_over_ice_dT(temperature(ji)) &
          *phase_trans_heat(1,temperature(ji))/c_v_mass_weighted_air(rho,temperature,ji-1)
          c = water_vapour_pressure - saturation_pressure
          p = b/a
          q = c/a
          diff_density = -0.5_wp*p - (0.25_wp*p**2 - q)**0.5_wp
          phase_trans_density = min(rho(2*n_scalars+ji), diff_density)
                
          ! the tendency for the water vapour
          phase_trans_rates(4*n_scalars+ji) = phase_trans_density/dtime
          
          ! the tendency for the ice contains two terms:
          ! 1.) the freezing
          ! 2.) the phase transition through sublimation
          phase_trans_rates(2*n_scalars+ji) = rho(3*n_scalars+ji)/dtime - phase_trans_density/dtime
                
          ! the heat source rates
          phase_trans_heating_rate(ji) &
          ! the freezing
          = -phase_trans_rates(3*n_scalars+ji)*phase_trans_heat(2,temperature(ji)) &
          ! the sublimation
          - phase_trans_density*phase_trans_heat(1,temperature(ji))/dtime
        endif
      ! the case where the air is over-saturated
      else
        ! temperature>=0° C
        if (temperature(ji)>=t_0) then
          ! It is assumed that the still present ice vanishes within one time step.
          phase_trans_rates(2*n_scalars+ji) = -rho(2*n_scalars+ji)/dtime
                
          ! the vanishing of water vapour through the phase transition
          a = -r_v*phase_trans_heat(0,temperature(ji))/c_v_mass_weighted_air(rho,temperature,ji-1)
          b = r_v*temperature(ji) - r_v*rho((n_condensed_constituents+1)*n_scalars+ji) &
          *phase_trans_heat(0,temperature(ji))/c_v_mass_weighted_air(rho,temperature,ji-1) &
          + enhancement_factor*dsaturation_pressure_over_water_dT(temperature(ji)) &
          *phase_trans_heat(0,temperature(ji))/c_v_mass_weighted_air(rho,temperature,ji-1)
          c = water_vapour_pressure - saturation_pressure
          p = b/a
          q = c/a
          diff_density = -0.5_wp*p - (0.25_wp*p**2 - q)**0.5_wp
                
          ! the tendency for the water vapour
          phase_trans_rates(4*n_scalars+ji) = diff_density/dtime
                
          ! The source rate for the liquid water consists of two terms:
          ! 1.) the melting
          ! 2.) the condensation
          phase_trans_rates(3*n_scalars+ji) = rho(2*n_scalars+ji)/dtime - diff_density/dtime
                
          ! the heat source rates
          phase_trans_heating_rate(ji) &
          ! melting
          = phase_trans_rates(2*n_scalars+ji)*phase_trans_heat(2,temperature(ji)) &
          ! condensation
          - diff_density*phase_trans_heat(0,temperature(ji))/dtime
        ! temperature<0° C
        else   
          ! It is assumed that the liquid water disappears within one time step.
          phase_trans_rates(3*n_scalars+ji) = -rho(3*n_scalars+ji)/dtime
                
          ! the vanishing of water vapour through the phase transition
          a = -r_v*phase_trans_heat(1,temperature(ji))/c_v_mass_weighted_air(rho,temperature,ji-1)
          b = r_v*temperature(ji) - r_v*rho((n_condensed_constituents+1)*n_scalars+ji) &
          *phase_trans_heat(1,temperature(ji))/c_v_mass_weighted_air(rho,temperature,ji-1) &
          + enhancement_factor*dsaturation_pressure_over_ice_dT(temperature(ji)) &
          *phase_trans_heat(1,temperature(ji))/c_v_mass_weighted_air(rho,temperature,ji-1)
          c = water_vapour_pressure - saturation_pressure
          p = b/a
          q = c/a
          diff_density = -0.5_wp*p - (0.25_wp*p**2 - q)**0.5_wp
               
          ! the tendency for the water vapour
          phase_trans_rates(4*n_scalars+ji) = diff_density/dtime
            
          ! The source rate for the cloud ice consists of two terms:
          ! 1.) the freezing
          ! 2.) the resublimation
                
          phase_trans_rates(2*n_scalars+ji) = rho(3*n_scalars+ji)/dtime - diff_density/dtime
                
          ! the heat source rates
          phase_trans_heating_rate(ji) &
          ! freezing
          = -phase_trans_rates(3*n_scalars+ji)*phase_trans_heat(2,temperature(ji)) &
          ! resublimation
          - diff_density*phase_trans_heat(1,temperature(ji))/dtime
        endif
      endif
      
      ! Precipitation
      ! -------------
      phase_trans_rates(ji) = 0._wp
      phase_trans_rates(n_scalars+ji) = 0._wp
      ! snow
      if (temperature(ji)<t_0) then
        phase_trans_rates(ji) = max(rho(2*n_scalars+ji) - maximum_cloud_water_content*rho(4*n_scalars+ji),0._wp)/1000._wp
        ! the snow creation comes at the cost of cloud ice particles
        phase_trans_rates(2*n_scalars+ji) = phase_trans_rates(2*n_scalars+ji) - phase_trans_rates(ji)
      ! rain
      elseif (temperature(ji)>=t_0) then
        phase_trans_rates(n_scalars+ji) = max(rho(3*n_scalars+ji) &
                                          - maximum_cloud_water_content*rho(4*n_scalars+ji),0._wp)/1000._wp
        ! the rain creation comes at the cost of cloud water particles
        phase_trans_rates(3*n_scalars+ji) = phase_trans_rates(3*n_scalars+ji) - phase_trans_rates(n_scalars+ji)
      endif
        
      ! turning of snow to rain
      if (temperature(ji)>=t_0 .and. rho(ji)>0._wp) then
        phase_trans_rates(ji) = -rho(ji)/dtime
        phase_trans_rates(n_scalars+ji) = phase_trans_rates(n_scalars+ji) - phase_trans_rates(ji)
        phase_trans_heating_rate(ji) = phase_trans_heating_rate(ji)+phase_trans_rates(ji)*phase_trans_heat(2,temperature(ji))
      endif
      ! turning of rain to snow
      if (temperature(ji)<t_0 .and. rho(n_scalars+ji)>0._wp) then
        phase_trans_rates(n_scalars+ji) = -rho(n_scalars+ji)/dtime
        phase_trans_rates(ji) = phase_trans_rates(ji) - phase_trans_rates(n_scalars+ji)
        phase_trans_heating_rate(ji) = phase_trans_heating_rate(ji) - &
                                       phase_trans_rates(n_scalars+ji)*phase_trans_heat(2,temperature(ji))
      endif
        
      ! Surface effects
      ! ---------------
      if (layer_index==n_layers - 1 .and. lsfc_phase_trans) then
        h_index = ji - layer_index*n_scalars_h
        
        ! evaporation and latent heat rates
        if (is_land(h_index)==0) then
          ! saturation pressure at surface temperature
          if (temperature_soil(h_index)>=t_0) then
            saturation_pressure_sfc = saturation_pressure_over_water(temperature_soil(h_index))
            saturation_pressure_sfc = enhancement_factor_over_water(air_pressure)*saturation_pressure_sfc
          else
            saturation_pressure_sfc = saturation_pressure_over_ice(temperature_soil(h_index))
            saturation_pressure_sfc = enhancement_factor_over_ice(air_pressure)*saturation_pressure_sfc
          endif
            
          ! difference water vapour density between saturation at ground temperature and actual absolute humidity in the lowest model layer
          diff_density_sfc = saturation_pressure_sfc/(r_v*temperature_soil(h_index)) &
          - rho((n_condensed_constituents+1)*n_scalars+ji)
            
          ! evporation, sublimation
          phase_trans_rates(n_condensed_constituents*n_scalars+ji) = &
          phase_trans_rates(n_condensed_constituents*n_scalars+ji) + &
          max(0._wp,diff_density_sfc/scalar_flux_resistance(h_index))/layer_thickness(ji)
          
          ! calculating the latent heat flux density affecting the surface
          if (temperature_soil(h_index)>=t_0) then
            power_flux_density_latent(h_index) = -phase_trans_heat(0,temperature_soil(h_index)) &
            *max(0._wp, diff_density_sfc/scalar_flux_resistance(h_index))
          else
            power_flux_density_latent(h_index) = -phase_trans_heat(1,temperature_soil(h_index)) &
            *max(0._wp, diff_density_sfc/scalar_flux_resistance(h_index))
          endif
        endif
      endif
    enddo
  
  end subroutine calc_h2otracers_source_rates

end module phase_trans












