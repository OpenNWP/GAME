! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module derived_quantities

  ! This module contains look-up functions for properties of the atmosphere.

  use definitions, only: wp
  use constants,   only: m_d,n_a,k_b,M_PI
  
  implicit none
  
  private
  
  public :: calc_diffusion_coeff
  
  contains
  
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
