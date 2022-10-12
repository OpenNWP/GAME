! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https:!github.com/OpenNWP/GAME

module mo_derived

  ! This file contains functions calculating derived thermodynamic quantities of the atmosphere.

  use mo_definitions,      only: wp,t_grid,t_state,t_diag
  use mo_constants,        only: m_d,n_a,k_b,M_PI,t_0,r_v,c_d_v,c_v_v,r_d,m_v
  use mo_dictionary,       only: saturation_pressure_over_water,saturation_pressure_over_ice,c_p_cond
  use mo_grid_nml,         only: n_cells,n_layers
  use mo_constituents_nml, only: lmoist,n_constituents,n_condensed_constituents
  
  implicit none
  
  contains

  subroutine temperature_diagnostics(state,diag,grid)
    
    ! This subroutine diagnoses the temperature of the gas phase.
    
    type(t_state), intent(in)    :: state ! state variables
    type(t_diag),  intent(inout) :: diag  ! diagnostic quantities
    type(t_grid),  intent(in)    :: grid  ! grid quantities (needed for the background state)
    
    if (.not. lmoist) then
      !$omp parallel workshare
      diag%temperature = (grid%theta_v_bg+state%theta_v_pert)*(grid%exner_bg+state%exner_pert)
      !$omp end parallel workshare
    endif
    if (lmoist) then
      !$omp parallel workshare
      diag%temperature = (grid%theta_v_bg+state%theta_v_pert)*(grid%exner_bg+state%exner_pert) &
      /(1._wp+state%rho(:,:,n_condensed_constituents+2) &
      /state%rho(:,:,n_condensed_constituents+1)*(m_d/m_v-1._wp))
      !$omp end parallel workshare
    endif
    
  end subroutine temperature_diagnostics
  
  function rel_humidity(abs_humidity,temperature)
    
    ! This function returns the relative humidity as a function of the absolute humidity in kg/m^3 and the temperature in K.
    
    real(wp), intent(in) :: abs_humidity ! absolute humidity (water vapour mass density (kg/m**3))
    real(wp), intent(in) :: temperature  ! temperature (K)
    real(wp)             :: rel_humidity ! result
    
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

  function gas_constant_diagnostics(rho,ji,jl)
    
    ! This function calculates the specific gas constant of the gas phase.
    
    real(wp), intent(in) :: rho(n_cells,n_layers,n_constituents)
    integer,  intent(in) :: ji,jl
    real(wp)             :: gas_constant_diagnostics
    
    gas_constant_diagnostics = 0._wp
    ! in the dry case this is just the individual gas constant of dry air
    if (.not. lmoist) then
      gas_constant_diagnostics = r_d
    endif
    if (lmoist) then
      gas_constant_diagnostics = (rho(ji,jl,n_condensed_constituents+1) &
                                  - rho(ji,jl,n_condensed_constituents+2))*r_d &
                                  + rho(ji,jl,n_condensed_constituents+2)*r_v
      gas_constant_diagnostics = gas_constant_diagnostics/rho(ji,jl,n_condensed_constituents+1)
    endif
  
  end function 

  function c_v_mass_weighted_air(rho,temperature,ji,jl)
  
    ! This function calculates the mass-weighted c_v of the air.
    
    real(wp), intent(in) :: rho(n_cells,n_layers,n_constituents),temperature(n_cells,n_layers)
    integer,  intent(in) :: ji,jl
    real(wp)             :: c_v_mass_weighted_air
    
    ! local variables
    integer :: jc ! constituent index
    
    c_v_mass_weighted_air = 0._wp
    do jc=1,n_condensed_constituents
      ! It is correct to use c_p here because the compression of the condensates has almost no effect on the air pressure.
      c_v_mass_weighted_air = c_v_mass_weighted_air + rho(ji,jl,jc)*c_p_cond(jc,temperature(ji,jl))
    enddo
    if (lmoist) then
      ! dry air
      c_v_mass_weighted_air = c_v_mass_weighted_air &
                              + (rho(ji,jl,n_condensed_constituents+1) - rho(ji,jl,n_condensed_constituents+2))*c_d_v
      ! water vapour
      c_v_mass_weighted_air = c_v_mass_weighted_air + rho(ji,jl,n_condensed_constituents+2)*c_v_v
    else
      c_v_mass_weighted_air = c_v_mass_weighted_air + rho(ji,jl,n_condensed_constituents+1)*c_d_v
    endif
  
  end function c_v_mass_weighted_air
  
  function calc_diffusion_coeff(temperature,density)
  
    ! This function calculates the molecular diffusion coefficient.
  
    real(wp), intent(in) :: temperature          ! temperature
    real(wp), intent(in) :: density              ! mass density
    real(wp)             :: calc_diffusion_coeff ! result
    
    ! local variables
    real(wp) :: particle_radius,particle_mass,thermal_velocity,particle_density,cross_section,mean_free_path

    ! these things are hardly ever modified
    particle_radius = 130e-12_wp
    particle_mass = m_d/n_a
    
    ! actual calculation
    thermal_velocity = sqrt(8._wp*k_b*temperature/(M_PI*particle_mass))
    particle_density = density/particle_mass
    cross_section = 4._wp*M_PI*particle_radius**2
    mean_free_path = 1._wp/(sqrt(2._wp)*particle_density*cross_section)
    calc_diffusion_coeff = 1._wp/3._wp*thermal_velocity*mean_free_path
    
  end function calc_diffusion_coeff

end module mo_derived













