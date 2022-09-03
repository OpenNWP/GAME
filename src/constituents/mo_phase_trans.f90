! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_phase_trans

  ! This file contains functions calculating everything related to phase transition rates.

  use mo_definitions,            only: wp,t_grid,t_state,t_diag
  use mo_run_nml,                only: dtime
  use mo_grid_nml,               only: n_scalars,n_scalars_h,n_layers
  use mo_constants,              only: r_v,t_0,r_d
  use mo_constituents_nml,       only: n_condensed_constituents,n_constituents
  use mo_dictionary,             only: saturation_pressure_over_water,saturation_pressure_over_ice, &
                                       dsaturation_pressure_over_water_dT,dsaturation_pressure_over_ice_dT, &
                                       phase_trans_heat,enhancement_factor_over_water,enhancement_factor_over_ice
  use mo_derived,                only: c_v_mass_weighted_air
  use mo_surface_nml,            only: nsoillays,lsfc_phase_trans

  implicit none
  
  contains

  subroutine calc_h2otracers_source_rates(state,diag,grid)
    
    ! This subroutine calculates phase transition rates and associated heat source rates.
    ! It assumes the following order for the constituents:
    ! precipitating ice - precipitating liquid water - cloud ice - liquid cloud water - moist air - water vapour
    
    type(t_state), intent(in)    :: state
    type(t_diag),  intent(inout) :: diag
    type(t_grid),  intent(in)    :: grid
    
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
      if (diag%temperature(ji)>=t_0) then
        saturation_pressure = saturation_pressure_over_water(diag%temperature(ji))
      ! "negative" temperatures
      else
        saturation_pressure = saturation_pressure_over_ice(diag%temperature(ji))
      endif
      
      ! determining the water vapour pressure (using the EOS)
      water_vapour_pressure = state%rho((n_condensed_constituents+1)*n_scalars+ji)*r_v*diag%temperature(ji)
      
      ! determining the water vapour pressure (using the EOS)
      dry_pressure = (state%rho(n_condensed_constituents*n_scalars+ji) - state%rho((n_condensed_constituents+1)*n_scalars+ji)) &
      *r_d*diag%temperature(ji)
        
      ! calculating the total air pressure
      air_pressure = dry_pressure + water_vapour_pressure
        
      ! multiplying the saturation pressure by the enhancement factor
      if (diag%temperature(ji)>=t_0) then
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
        if (diag%temperature(ji)>=t_0) then
          ! It is assumed that the still present ice vanishes within one time step.
          diag%phase_trans_rates(2*n_scalars+ji) = -state%rho(2*n_scalars+ji)/dtime
                
          ! The amount of liquid water per volume that will evaporate.
          ! In case the air cannot take all the water, not everything will evaporate.
          a = -r_v*phase_trans_heat(0,diag%temperature(ji))/c_v_mass_weighted_air(state%rho,diag%temperature,ji-1)
          b = r_v*diag%temperature(ji) - r_v*state%rho((n_condensed_constituents+1)*n_scalars+ji) &
          *phase_trans_heat(0,diag%temperature(ji))/c_v_mass_weighted_air(state%rho,diag%temperature,ji-1) &
          + enhancement_factor*dsaturation_pressure_over_water_dT(diag%temperature(ji)) &
          *phase_trans_heat(0,diag%temperature(ji))/c_v_mass_weighted_air(state%rho,diag%temperature,ji-1)
          c = water_vapour_pressure - saturation_pressure
          p = b/a
          q = c/a
          diff_density = -0.5_wp*p - (0.25_wp*p**2 - q)**0.5_wp
          phase_trans_density = min(state%rho(3*n_scalars+ji), diff_density)
                
          ! the tendency for the water vapour
          diag%phase_trans_rates(4*n_scalars+ji) = phase_trans_density/dtime
                
          ! The source rate for the liquid water consists of two terms:
          ! 1.) the melting
          ! 2.) the evaporation
                
          diag%phase_trans_rates(3*n_scalars+ji) = state%rho(2*n_scalars+ji)/dtime - phase_trans_density/dtime
                
          ! the heat source rates
          diag%phase_trans_heating_rate(ji) &
          ! melting
          = diag%phase_trans_rates(2*n_scalars+ji)*phase_trans_heat(2,diag%temperature(ji)) &
          ! evaporation
          - phase_trans_density*phase_trans_heat(0,diag%temperature(ji))/dtime
        ! temperature<0° C
        else
          ! It is assumed that the still present liquid water vanishes within one time step.
          diag%phase_trans_rates(3*n_scalars+ji) = -state%rho(3*n_scalars+ji)/dtime
                
          ! The amount of ice per volume that will sublimate.
          ! In case the air cannot take all the water, not everything will sublimate.
                
          a = -r_v*phase_trans_heat(1,diag%temperature(ji))/c_v_mass_weighted_air(state%rho,diag%temperature,ji-1)
          b = r_v*diag%temperature(ji) - r_v*state%rho((n_condensed_constituents+1)*n_scalars+ji) &
          *phase_trans_heat(1,diag%temperature(ji))/c_v_mass_weighted_air(state%rho,diag%temperature,ji-1) &
          + enhancement_factor*dsaturation_pressure_over_ice_dT(diag%temperature(ji)) &
          *phase_trans_heat(1,diag%temperature(ji))/c_v_mass_weighted_air(state%rho,diag%temperature,ji-1)
          c = water_vapour_pressure - saturation_pressure
          p = b/a
          q = c/a
          diff_density = -0.5_wp*p - (0.25_wp*p**2 - q)**0.5_wp
          phase_trans_density = min(state%rho(2*n_scalars+ji), diff_density)
                
          ! the tendency for the water vapour
          diag%phase_trans_rates(4*n_scalars+ji) = phase_trans_density/dtime
          
          ! the tendency for the ice contains two terms:
          ! 1.) the freezing
          ! 2.) the phase transition through sublimation
          diag%phase_trans_rates(2*n_scalars+ji) = state%rho(3*n_scalars+ji)/dtime - phase_trans_density/dtime
                
          ! the heat source rates
          diag%phase_trans_heating_rate(ji) &
          ! the freezing
          = -diag%phase_trans_rates(3*n_scalars+ji)*phase_trans_heat(2,diag%temperature(ji)) &
          ! the sublimation
          - phase_trans_density*phase_trans_heat(1,diag%temperature(ji))/dtime
        endif
      ! the case where the air is over-saturated
      else
        ! temperature>=0° C
        if (diag%temperature(ji)>=t_0) then
          ! It is assumed that the still present ice vanishes within one time step.
          diag%phase_trans_rates(2*n_scalars+ji) = -state%rho(2*n_scalars+ji)/dtime
                
          ! the vanishing of water vapour through the phase transition
          a = -r_v*phase_trans_heat(0,diag%temperature(ji))/c_v_mass_weighted_air(state%rho,diag%temperature,ji-1)
          b = r_v*diag%temperature(ji) - r_v*state%rho((n_condensed_constituents+1)*n_scalars+ji) &
          *phase_trans_heat(0,diag%temperature(ji))/c_v_mass_weighted_air(state%rho,diag%temperature,ji-1) &
          + enhancement_factor*dsaturation_pressure_over_water_dT(diag%temperature(ji)) &
          *phase_trans_heat(0,diag%temperature(ji))/c_v_mass_weighted_air(state%rho,diag%temperature,ji-1)
          c = water_vapour_pressure - saturation_pressure
          p = b/a
          q = c/a
          diff_density = -0.5_wp*p - (0.25_wp*p**2 - q)**0.5_wp
                
          ! the tendency for the water vapour
          diag%phase_trans_rates(4*n_scalars+ji) = diff_density/dtime
                
          ! The source rate for the liquid water consists of two terms:
          ! 1.) the melting
          ! 2.) the condensation
          diag%phase_trans_rates(3*n_scalars+ji) = state%rho(2*n_scalars+ji)/dtime - diff_density/dtime
                
          ! the heat source rates
          diag%phase_trans_heating_rate(ji) &
          ! melting
          = diag%phase_trans_rates(2*n_scalars+ji)*phase_trans_heat(2,diag%temperature(ji)) &
          ! condensation
          - diff_density*phase_trans_heat(0,diag%temperature(ji))/dtime
        ! temperature<0° C
        else   
          ! It is assumed that the liquid water disappears within one time step.
          diag%phase_trans_rates(3*n_scalars+ji) = -state%rho(3*n_scalars+ji)/dtime
                
          ! the vanishing of water vapour through the phase transition
          a = -r_v*phase_trans_heat(1,diag%temperature(ji))/c_v_mass_weighted_air(state%rho,diag%temperature,ji-1)
          b = r_v*diag%temperature(ji) - r_v*state%rho((n_condensed_constituents+1)*n_scalars+ji) &
          *phase_trans_heat(1,diag%temperature(ji))/c_v_mass_weighted_air(state%rho,diag%temperature,ji-1) &
          + enhancement_factor*dsaturation_pressure_over_ice_dT(diag%temperature(ji)) &
          *phase_trans_heat(1,diag%temperature(ji))/c_v_mass_weighted_air(state%rho,diag%temperature,ji-1)
          c = water_vapour_pressure - saturation_pressure
          p = b/a
          q = c/a
          diff_density = -0.5_wp*p - (0.25_wp*p**2 - q)**0.5_wp
               
          ! the tendency for the water vapour
          diag%phase_trans_rates(4*n_scalars+ji) = diff_density/dtime
            
          ! The source rate for the cloud ice consists of two terms:
          ! 1.) the freezing
          ! 2.) the resublimation
                
          diag%phase_trans_rates(2*n_scalars+ji) = state%rho(3*n_scalars+ji)/dtime - diff_density/dtime
                
          ! the heat source rates
          diag%phase_trans_heating_rate(ji) &
          ! freezing
          = -diag%phase_trans_rates(3*n_scalars+ji)*phase_trans_heat(2,diag%temperature(ji)) &
          ! resublimation
          - diff_density*phase_trans_heat(1,diag%temperature(ji))/dtime
        endif
      endif
      
      ! Precipitation
      ! -------------
      diag%phase_trans_rates(ji) = 0._wp
      diag%phase_trans_rates(n_scalars+ji) = 0._wp
      ! snow
      if (diag%temperature(ji)<t_0) then
        diag%phase_trans_rates(ji) = max(state%rho(2*n_scalars+ji) &
                                         - maximum_cloud_water_content*state%rho(4*n_scalars+ji),0._wp)/1000._wp
        ! the snow creation comes at the cost of cloud ice particles
        diag%phase_trans_rates(2*n_scalars+ji) = diag%phase_trans_rates(2*n_scalars+ji) - diag%phase_trans_rates(ji)
      ! rain
      elseif (diag%temperature(ji)>=t_0) then
        diag%phase_trans_rates(n_scalars+ji) = max(state%rho(3*n_scalars+ji) &
                                          - maximum_cloud_water_content*state%rho(4*n_scalars+ji),0._wp)/1000._wp
        ! the rain creation comes at the cost of cloud water particles
        diag%phase_trans_rates(3*n_scalars+ji) = diag%phase_trans_rates(3*n_scalars+ji) - diag%phase_trans_rates(n_scalars+ji)
      endif
        
      ! turning of snow to rain
      if (diag%temperature(ji)>=t_0 .and. state%rho(ji)>0._wp) then
        diag%phase_trans_rates(ji) = -state%rho(ji)/dtime
        diag%phase_trans_rates(n_scalars+ji) = diag%phase_trans_rates(n_scalars+ji) - diag%phase_trans_rates(ji)
        diag%phase_trans_heating_rate(ji) = diag%phase_trans_heating_rate(ji) &
                                            + diag%phase_trans_rates(ji)*phase_trans_heat(2,diag%temperature(ji))
      endif
      ! turning of rain to snow
      if (diag%temperature(ji)<t_0 .and. state%rho(n_scalars+ji)>0._wp) then
        diag%phase_trans_rates(n_scalars+ji) = -state%rho(n_scalars+ji)/dtime
        diag%phase_trans_rates(ji) = diag%phase_trans_rates(ji) - diag%phase_trans_rates(n_scalars+ji)
        diag%phase_trans_heating_rate(ji) = diag%phase_trans_heating_rate(ji) - &
                                       diag%phase_trans_rates(n_scalars+ji)*phase_trans_heat(2,diag%temperature(ji))
      endif
        
      ! Surface effects
      ! ---------------
      if (layer_index==n_layers - 1 .and. lsfc_phase_trans) then
        h_index = ji - layer_index*n_scalars_h
        
        ! evaporation and latent heat rates
        if (grid%is_land(h_index)==0) then
          ! saturation pressure at surface temperature
          if (state%temperature_soil(h_index)>=t_0) then
            saturation_pressure_sfc = saturation_pressure_over_water(state%temperature_soil(h_index))
            saturation_pressure_sfc = enhancement_factor_over_water(air_pressure)*saturation_pressure_sfc
          else
            saturation_pressure_sfc = saturation_pressure_over_ice(state%temperature_soil(h_index))
            saturation_pressure_sfc = enhancement_factor_over_ice(air_pressure)*saturation_pressure_sfc
          endif
            
          ! difference water vapour density between saturation at ground temperature and actual absolute humidity in the lowest model layer
          diff_density_sfc = saturation_pressure_sfc/(r_v*state%temperature_soil(h_index)) &
          - state%rho((n_condensed_constituents+1)*n_scalars+ji)
            
          ! evporation, sublimation
          diag%phase_trans_rates(n_condensed_constituents*n_scalars+ji) = &
          diag%phase_trans_rates(n_condensed_constituents*n_scalars+ji) + &
          max(0._wp,diff_density_sfc/diag%scalar_flux_resistance(h_index))/grid%layer_thickness(ji)
          
          ! calculating the latent heat flux density affecting the surface
          if (state%temperature_soil(h_index)>=t_0) then
            diag%power_flux_density_latent(h_index) = -phase_trans_heat(0,state%temperature_soil(h_index)) &
            *max(0._wp, diff_density_sfc/diag%scalar_flux_resistance(h_index))
          else
            diag%power_flux_density_latent(h_index) = -phase_trans_heat(1,state%temperature_soil(h_index)) &
            *max(0._wp, diff_density_sfc/diag%scalar_flux_resistance(h_index))
          endif
        endif
      endif
    enddo
  
  end subroutine calc_h2otracers_source_rates

end module mo_phase_trans












