! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_phase_trans

  ! This file contains functions calculating everything related to phase transition rates.

  use mo_definitions,      only: wp,t_grid,t_state,t_diag
  use mo_grid_nml,         only: n_cells,n_layers
  use mo_constants,        only: r_v,t_0,r_d
  use mo_constituents_nml, only: n_condensed_constituents,n_constituents
  use mo_grid_setup,       only: dtime
  use mo_dictionary,       only: saturation_pressure_over_water,saturation_pressure_over_ice,rain_drops_radius, &
                                 dsaturation_pressure_over_water_dT,dsaturation_pressure_over_ice_dT, &
                                 phase_trans_heat,enhancement_factor_over_water,enhancement_factor_over_ice
  use mo_derived,          only: c_v_mass_weighted_air
  use mo_surface_nml,      only: nsoillays,lsfc_phase_trans

  implicit none
  
  contains

  subroutine calc_h2otracers_source_rates(state,diag,grid)
    
    ! This subroutine calculates phase transition rates and associated heat source rates.
    ! It assumes the following order for the constituents:
    ! precipitating ice - precipitating liquid water - cloud ice - liquid cloud water - moist air - water vapour
    
    type(t_state), intent(in)    :: state ! the state which to use for computing the phase transition rates
    type(t_diag),  intent(inout) :: diag  ! the diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities
    
    ! local variables
    integer  :: ji                          ! cell index
    integer  :: jl                          ! layer index
    real(wp) :: diff_density                ! difference between saturation water vapour density and actual water vapour density and 
    real(wp) :: phase_trans_density         ! actual phase transition density
    real(wp) :: saturation_pressure         ! saturation water vapour pressure
    real(wp) :: water_vapour_pressure       ! actual water vapour pressure
    real(wp) :: diff_density_sfc            ! diff_density at the surface
    real(wp) :: saturation_pressure_sfc     ! saturation water vapour pressure at the surface
    real(wp) :: dry_pressure                ! dry air pressure
    real(wp) :: air_pressure                ! complete air pressure
    real(wp) :: a                           ! helper variable for computing the second-order phase transition rates
    real(wp) :: b                           ! helper variable for computing the second-order phase transition rates
    real(wp) :: c                           ! helper variable for computing the second-order phase transition rates
    real(wp) :: p                           ! helper variable for computing the second-order phase transition rates
    real(wp) :: q                           ! helper variable for computing the second-order phase transition rates
    real(wp) :: enhancement_factor          ! factor taking into account non-ideal effects of air
    real(wp) :: maximum_cloud_water_content ! maximum cloud water content in (kg cloud)/(kg dry air)
    
    maximum_cloud_water_content = 0.2e-3_wp
    
    ! loop over all grid boxes
    !$omp parallel do private(ji,jl,diff_density,phase_trans_density,saturation_pressure,water_vapour_pressure, &
    !$omp diff_density_sfc,saturation_pressure_sfc,dry_pressure,air_pressure, &
    !$omp a,b,c,p,q,enhancement_factor)
    do jl=1,n_layers
      do ji=1,n_cells
        
        ! Preparation
        ! -----------
        
        ! determining the radius of rain drops
        diag%a_rain(ji,jl) = rain_drops_radius(state%rho(ji,jl,2)+state%rho(ji,jl,4))
        
        ! determining the saturation pressure
        ! "positive" temperatures (the saturation pressure is different over water compared to over ice)
        if (diag%temperature(ji,jl)>=t_0) then
          saturation_pressure = saturation_pressure_over_water(diag%temperature(ji,jl))
        ! "negative" temperatures
        else
          saturation_pressure = saturation_pressure_over_ice(diag%temperature(ji,jl))
        endif
        
        ! determining the water vapour pressure (using the EOS)
        water_vapour_pressure = state%rho(ji,jl,n_condensed_constituents+2)*r_v*diag%temperature(ji,jl)
        
        ! determining the water vapour pressure (using the EOS)
        dry_pressure = (state%rho(ji,jl,n_condensed_constituents+1) - state%rho(ji,jl,n_condensed_constituents+2)) &
        *r_d*diag%temperature(ji,jl)
          
        ! calculating the total air pressure
        air_pressure = dry_pressure + water_vapour_pressure
          
        ! multiplying the saturation pressure by the enhancement factor
        if (diag%temperature(ji,jl)>=t_0) then
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
          if (diag%temperature(ji,jl)>=t_0) then
            ! It is assumed that the still present ice vanishes within one time step.
            diag%phase_trans_rates(ji,jl,3) = -state%rho(ji,jl,3)/dtime
                  
            ! The amount of liquid water per volume that will evaporate.
            ! In case the air cannot take all the water, not everything will evaporate.
            a = -r_v*phase_trans_heat(1,diag%temperature(ji,jl))/c_v_mass_weighted_air(state%rho,diag%temperature,ji,jl)
            b = r_v*diag%temperature(ji,jl) - r_v*state%rho(ji,jl,n_condensed_constituents+2) &
            *phase_trans_heat(1,diag%temperature(ji,jl))/c_v_mass_weighted_air(state%rho,diag%temperature,ji,jl) &
            + enhancement_factor*dsaturation_pressure_over_water_dT(diag%temperature(ji,jl)) &
            *phase_trans_heat(1,diag%temperature(ji,jl))/c_v_mass_weighted_air(state%rho,diag%temperature,ji,jl)
            c = water_vapour_pressure - saturation_pressure
            p = b/a
            q = c/a
            diff_density = -0.5_wp*p - (0.25_wp*p**2 - q)**0.5_wp
            phase_trans_density = min(state%rho(ji,jl,4),diff_density)
            
            ! the tendency for the water vapour
            diag%phase_trans_rates(ji,jl,5) = phase_trans_density/dtime
            
            ! The source rate for the liquid water consists of two terms:
            ! 1.) the melting
            ! 2.) the evaporation
            
            diag%phase_trans_rates(ji,jl,4) = state%rho(ji,jl,3)/dtime - phase_trans_density/dtime
            
            ! the heat source rates
            diag%phase_trans_heating_rate(ji,jl) &
            ! melting
            = diag%phase_trans_rates(ji,jl,3)*phase_trans_heat(3,diag%temperature(ji,jl)) &
            ! evaporation
            - phase_trans_density*phase_trans_heat(1,diag%temperature(ji,jl))/dtime
          ! temperature<0° C
          else
            ! It is assumed that the still present liquid water vanishes within one time step.
            diag%phase_trans_rates(ji,jl,4) = -state%rho(ji,jl,4)/dtime
            
            ! The amount of ice per volume that will sublimate.
            ! In case the air cannot take all the water, not everything will sublimate.
            
            a = -r_v*phase_trans_heat(2,diag%temperature(ji,jl))/c_v_mass_weighted_air(state%rho,diag%temperature,ji,jl)
            b = r_v*diag%temperature(ji,jl) - r_v*state%rho(ji,jl,n_condensed_constituents+2) &
            *phase_trans_heat(2,diag%temperature(ji,jl))/c_v_mass_weighted_air(state%rho,diag%temperature,ji,jl) &
            + enhancement_factor*dsaturation_pressure_over_ice_dT(diag%temperature(ji,jl)) &
            *phase_trans_heat(2,diag%temperature(ji,jl))/c_v_mass_weighted_air(state%rho,diag%temperature,ji,jl)
            c = water_vapour_pressure - saturation_pressure
            p = b/a
            q = c/a
            diff_density = -0.5_wp*p - (0.25_wp*p**2 - q)**0.5_wp
            phase_trans_density = min(state%rho(ji,jl,3), diff_density)
            
            ! the tendency for the water vapour
            diag%phase_trans_rates(ji,jl,5) = phase_trans_density/dtime
            
            ! the tendency for the ice contains two terms:
            ! 1.) the freezing
            ! 2.) the phase transition through sublimation
            diag%phase_trans_rates(ji,jl,3) = state%rho(ji,jl,4)/dtime - phase_trans_density/dtime
            
            ! the heat source rates
            diag%phase_trans_heating_rate(ji,jl) &
            ! the freezing
            = -diag%phase_trans_rates(ji,jl,4)*phase_trans_heat(3,diag%temperature(ji,jl)) &
            ! the sublimation
            - phase_trans_density*phase_trans_heat(2,diag%temperature(ji,jl))/dtime
          endif
        ! the case where the air is over-saturated
        else
          ! temperature>=0° C
          if (diag%temperature(ji,jl)>=t_0) then
            ! It is assumed that the still present ice vanishes within one time step.
            diag%phase_trans_rates(ji,jl,3) = -state%rho(ji,jl,3)/dtime
            
            ! the vanishing of water vapour through the phase transition
            a = -r_v*phase_trans_heat(1,diag%temperature(ji,jl))/c_v_mass_weighted_air(state%rho,diag%temperature,ji,jl)
            b = r_v*diag%temperature(ji,jl) - r_v*state%rho(ji,jl,n_condensed_constituents+2) &
            *phase_trans_heat(1,diag%temperature(ji,jl))/c_v_mass_weighted_air(state%rho,diag%temperature,ji,jl) &
            + enhancement_factor*dsaturation_pressure_over_water_dT(diag%temperature(ji,jl)) &
            *phase_trans_heat(1,diag%temperature(ji,jl))/c_v_mass_weighted_air(state%rho,diag%temperature,ji,jl)
            c = water_vapour_pressure - saturation_pressure
            p = b/a
            q = c/a
            diff_density = -0.5_wp*p - (0.25_wp*p**2 - q)**0.5_wp
            
            ! the tendency for the water vapour
            diag%phase_trans_rates(ji,jl,5) = diff_density/dtime
            
            ! The source rate for the liquid water consists of two terms:
            ! 1.) the melting
            ! 2.) the condensation
            diag%phase_trans_rates(ji,jl,4) = state%rho(ji,jl,3)/dtime - diff_density/dtime
            
            ! the heat source rates
            diag%phase_trans_heating_rate(ji,jl) &
            ! melting
            = diag%phase_trans_rates(ji,jl,3)*phase_trans_heat(3,diag%temperature(ji,jl)) &
            ! condensation
            - diff_density*phase_trans_heat(1,diag%temperature(ji,jl))/dtime
          ! temperature<0° C
          else
            ! It is assumed that the liquid water disappears within one time step.
            diag%phase_trans_rates(ji,jl,4) = -state%rho(ji,jl,4)/dtime
            
            ! the vanishing of water vapour through the phase transition
            a = -r_v*phase_trans_heat(2,diag%temperature(ji,jl))/c_v_mass_weighted_air(state%rho,diag%temperature,ji,jl)
            b = r_v*diag%temperature(ji,jl) - r_v*state%rho(ji,jl,n_condensed_constituents+2) &
            *phase_trans_heat(2,diag%temperature(ji,jl))/c_v_mass_weighted_air(state%rho,diag%temperature,ji,jl) &
            + enhancement_factor*dsaturation_pressure_over_ice_dT(diag%temperature(ji,jl)) &
            *phase_trans_heat(2,diag%temperature(ji,jl))/c_v_mass_weighted_air(state%rho,diag%temperature,ji,jl)
            c = water_vapour_pressure - saturation_pressure
            p = b/a
            q = c/a
            diff_density = -0.5_wp*p - (0.25_wp*p**2 - q)**0.5_wp
            
            ! the tendency for the water vapour
            diag%phase_trans_rates(ji,jl,5) = diff_density/dtime
            
            ! The source rate for the cloud ice consists of two terms:
            ! 1.) the freezing
            ! 2.) the resublimation
            
            diag%phase_trans_rates(ji,jl,3) = state%rho(ji,jl,4)/dtime - diff_density/dtime
            
            ! the heat source rates
            diag%phase_trans_heating_rate(ji,jl) &
            ! freezing
            = -diag%phase_trans_rates(ji,jl,4)*phase_trans_heat(3,diag%temperature(ji,jl)) &
            ! resublimation
            - diff_density*phase_trans_heat(2,diag%temperature(ji,jl))/dtime
          endif
        endif
        
        ! Precipitation
        ! -------------
        diag%phase_trans_rates(ji,jl,1) = 0._wp
        diag%phase_trans_rates(ji,jl,2) = 0._wp
        ! snow
        if (diag%temperature(ji,jl)<t_0) then
          diag%phase_trans_rates(ji,jl,1) = max(state%rho(ji,jl,3) &
                                           - maximum_cloud_water_content*state%rho(ji,jl,5),0._wp)/1000._wp
          ! the snow creation comes at the cost of cloud ice particles
          diag%phase_trans_rates(ji,jl,3) = diag%phase_trans_rates(ji,jl,3) - diag%phase_trans_rates(ji,jl,1)
        ! rain
        elseif (diag%temperature(ji,jl)>=t_0) then
          diag%phase_trans_rates(ji,jl,2) = max(state%rho(ji,jl,4) &
                                            - maximum_cloud_water_content*state%rho(ji,jl,5),0._wp)/1000._wp
          ! the rain creation comes at the cost of cloud water particles
          diag%phase_trans_rates(ji,jl,4) = diag%phase_trans_rates(ji,jl,4) - diag%phase_trans_rates(ji,jl,2)
        endif
        
        ! turning of snow to rain
        if (diag%temperature(ji,jl)>=t_0 .and. state%rho(ji,jl,1)>0._wp) then
          diag%phase_trans_rates(ji,jl,1) = -state%rho(ji,jl,1)/dtime
          diag%phase_trans_rates(ji,jl,2) = diag%phase_trans_rates(ji,jl,2) - diag%phase_trans_rates(ji,jl,1)
          diag%phase_trans_heating_rate(ji,jl) = diag%phase_trans_heating_rate(ji,jl) &
                                              + diag%phase_trans_rates(ji,jl,1)*phase_trans_heat(3,diag%temperature(ji,jl))
        endif
        ! turning of rain to snow
        if (diag%temperature(ji,jl)<t_0 .and. state%rho(ji,jl,2)>0._wp) then
          diag%phase_trans_rates(ji,jl,2) = -state%rho(ji,jl,2)/dtime
          diag%phase_trans_rates(ji,jl,1) = diag%phase_trans_rates(ji,jl,1) - diag%phase_trans_rates(ji,jl,2)
          diag%phase_trans_heating_rate(ji,jl) = diag%phase_trans_heating_rate(ji,jl) - &
                                         diag%phase_trans_rates(ji,jl,2)*phase_trans_heat(3,diag%temperature(ji,jl))
        endif
        
        ! Surface effects
        ! ---------------
        if (jl==n_layers .and. lsfc_phase_trans) then
          
          ! evaporation and latent heat rates
          if (grid%land_fraction(ji)==0) then
            ! saturation pressure at surface temperature
            if (state%temperature_soil(ji,1)>=t_0) then
              saturation_pressure_sfc = saturation_pressure_over_water(state%temperature_soil(ji,1))
              saturation_pressure_sfc = enhancement_factor_over_water(air_pressure)*saturation_pressure_sfc
            else
              saturation_pressure_sfc = saturation_pressure_over_ice(state%temperature_soil(ji,1))
              saturation_pressure_sfc = enhancement_factor_over_ice(air_pressure)*saturation_pressure_sfc
            endif
            
            ! difference water vapour density between saturation at ground temperature and actual absolute humidity in the lowest model layer
            diff_density_sfc = saturation_pressure_sfc/(r_v*state%temperature_soil(ji,1)) &
            - state%rho(ji,jl,n_condensed_constituents+2)
            
            ! evporation, sublimation
            diag%phase_trans_rates(ji,jl,n_condensed_constituents+1) = &
            diag%phase_trans_rates(ji,jl,n_condensed_constituents+1) + &
            max(0._wp,diff_density_sfc/diag%scalar_flux_resistance(ji))/grid%layer_thickness(ji,jl)
            
            ! calculating the latent heat flux density affecting the surface
            if (state%temperature_soil(ji,1)>=t_0) then
              diag%power_flux_density_latent(ji) = -phase_trans_heat(1,state%temperature_soil(ji,1)) &
              *max(0._wp,diff_density_sfc/diag%scalar_flux_resistance(ji))
            else
              diag%power_flux_density_latent(ji) = -phase_trans_heat(2,state%temperature_soil(ji,1)) &
              *max(0._wp,diff_density_sfc/diag%scalar_flux_resistance(ji))
            endif
          endif
        endif
      enddo
    enddo
    
  end subroutine calc_h2otracers_source_rates

end module mo_phase_trans












