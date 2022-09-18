! This source file is part of the Limited-area GAME version (L-GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/L-GAME

module mo_derived

  ! This module contains look-up functions for properties of the atmosphere.

  use mo_definitions,      only: wp,t_grid,t_state,t_diag
  use mo_run_nml,          only: lmoist
  use mo_constants,        only: m_d,n_a,k_b,M_PI,t_0,r_v,c_d_v,c_v_v,r_d,m_v
  use mo_dictionary,       only: saturation_pressure_over_water,saturation_pressure_over_ice,c_p_cond
  use mo_grid_nml,         only: n_scalars
  use mo_constituents_nml, only: n_constituents,n_condensed_constituents
  
  implicit none
  
  contains
  
  function rel_humidity(abs_humidity,temperature)
    
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
  
  function calc_diffusion_coeff(temperature,density)
  
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
  
  function density_total(rho,ji)
  
    ! This function calculates the total density of the air at a certain gridpoint.
    
    ! input arguments
    real(wp),intent(in) :: rho(n_constituents*n_scalars)  ! density field
    integer, intent(in) :: ji   ! indices of the gridpoint
    ! output
    real(wp)            :: density_total ! the result
    
    ! local variables
    integer :: jc
    
    density_total = 0._wp
    
    do jc=0,n_constituents-1
      density_total = density_total + rho(1+ji+jc*n_scalars)
    enddo
    
  end function density_total

  subroutine temperature_diagnostics(state,diag,grid)
    
    ! This subroutine diagnoses the temperature of the gas phase.
    
    type(t_state), intent(in)    :: state
    type(t_diag),  intent(inout) :: diag
    type(t_grid),  intent(in)    :: grid
    
    ! local variables
    integer :: ji
    
    if (.not. lmoist) then
     !$omp parallel do private(ji)
      do ji=1,n_scalars
        diag%temperature(ji) = (grid%theta_v_bg(ji)+state%theta_v_pert(ji))*(grid%exner_bg(ji)+state%exner_pert(ji))
      enddo
      !$omp end parallel do
    endif
    if (lmoist) then
     !$omp parallel do private(ji)
      do ji=1,n_scalars
        diag%temperature(ji) = (grid%theta_v_bg(ji)+state%theta_v_pert(ji))*(grid%exner_bg(ji)+state%exner_pert(ji)) &
        /(1._wp+state%rho(ji,n_condensed_constituents+2) &
        /state%rho(ji,n_condensed_constituents+1)*(m_d/m_v-1._wp))
      enddo
      !$omp end parallel do
    endif
    
  end subroutine temperature_diagnostics

  function gas_constant_diagnostics(rho,grid_point_index)
    
    ! This function calculates the specific gas constant of the gas phase.
    
    real(wp), intent(in) :: rho(n_constituents*n_scalars)
    integer,  intent(in) :: grid_point_index
    real(wp)             :: gas_constant_diagnostics
    
    gas_constant_diagnostics = 0._wp
    if (.not. lmoist) then
      gas_constant_diagnostics = r_d
    endif
    if (lmoist) then
      gas_constant_diagnostics = (rho(n_condensed_constituents*n_scalars+1+grid_point_index) &
                                  - rho((n_condensed_constituents + 1)*n_scalars+1+grid_point_index))*r_d &
                                  + rho((n_condensed_constituents + 1)*n_scalars+1+grid_point_index)*r_v
      gas_constant_diagnostics = gas_constant_diagnostics/rho(n_condensed_constituents*n_scalars+1+grid_point_index)
    endif
  
  end function 

  function c_v_mass_weighted_air(rho,temperature,grid_point_index)
  
    ! This function calculates the mass-weighted c_v of the air.
    
    real(wp), intent(in) :: rho(n_constituents*n_scalars),temperature(n_scalars)
    integer,  intent(in) :: grid_point_index
    real(wp) :: c_v_mass_weighted_air
    
    ! local variables
    integer :: jc
    
    c_v_mass_weighted_air = 0._wp
    do jc=0,n_condensed_constituents-1
      ! It is correct to use c_p here because the compression of the condensates has almost no effect on the air pressure.
      c_v_mass_weighted_air = c_v_mass_weighted_air &
                              +rho(jc*n_scalars+1+grid_point_index) &
                              *c_p_cond(jc,temperature(1+grid_point_index))
    enddo
    c_v_mass_weighted_air = rho(n_condensed_constituents*n_scalars+1+grid_point_index)*c_d_v
    if (lmoist) then
      c_v_mass_weighted_air = c_v_mass_weighted_air &
                              +rho((n_condensed_constituents+1)*n_scalars+1+grid_point_index)*c_v_v
    endif
  
  end function c_v_mass_weighted_air

end module mo_derived













