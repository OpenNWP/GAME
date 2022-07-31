! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module derived_quantities

  ! This module contains look-up functions for properties of the atmosphere.

  use definitions, only: wp
  use constants,   only: m_d,n_a,k_b,M_PI,t_0,r_v
  use dictionary,  only: saturation_pressure_over_water,saturation_pressure_over_ice
  
  implicit none
  
  private
  
  public :: rel_humidity
  public :: calc_diffusion_coeff
  
  contains
  
  function rel_humidity(abs_humidity,temperature) &
  bind(c,name = "rel_humidity")
    
    ! This function returns the relative humidity as a function of the absolute humidity in kg/m^3 and the temperature in K.
    
    real(wp), intent(in) :: abs_humidity,temperature
    real(wp)             :: rel_humidity
    
    ! local variables
    real(wp)             :: vapour_pressure     ! actual water vapour pressure
    real(wp)             :: saturation_pressure ! saturation water vapour pressure
    
    ! calculation of the water vapour pressure according to the equation of state
    vapour_pressure = abs_humidity*r_v*temperature
    
    if (temperature>t_0) then
      saturation_pressure = saturation_pressure_over_water(temperature)
    endif
    if (temperature<=t_0) then
      saturation_pressure = saturation_pressure_over_ice(temperature)
    endif
    
    rel_humidity = vapour_pressure/saturation_pressure
    
  end function rel_humidity
  
  function calc_diffusion_coeff(temperature,density) &
  bind(c,name = "calc_diffusion_coeff")
  
    ! This function calculates the molecular diffusion coefficient.
  
    real(wp), intent(in) :: temperature,density
    real(wp)             :: calc_diffusion_coeff
    
    ! local variables
    real(wp) :: particle_radius,particle_mass,thermal_velocity,particle_density,cross_section,mean_free_path

    ! these things are hardly ever modified
    particle_radius = 130e-12_wp
    particle_mass = m_d/n_a
    
    ! actual calculation
    thermal_velocity = sqrt(8.0_wp*k_b*temperature/(M_PI*particle_mass))
    particle_density = density/particle_mass
    cross_section = 4.0_wp*M_PI*particle_radius**2
    mean_free_path = 1.0_wp/(sqrt(2.0_wp)*particle_density*cross_section)
    calc_diffusion_coeff = 1.0_wp/3.0_wp*thermal_velocity*mean_free_path
    
  end function calc_diffusion_coeff

end module derived_quantities













