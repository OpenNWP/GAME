! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_held_suarez
  
  ! This module calculates the Held-Suarez radiative forcing.
  
  use mo_definitions,      only: wp
  use mo_constants,        only: c_d_v,r_d,p_0
  use mo_constituents_nml, only: n_constituents,n_condensed_constituents
  use mo_grid_nml,         only: n_layers
  
  implicit none
  
  contains
  
  subroutine held_suar(latitude_scalar,mass_densities,temperature_gas,n_cells_rad,radiation_tendency)
    
    ! This subroutine calculates the Held-Suarez radiative forcing.
    
    integer,  intent(in)  :: n_cells_rad                                         ! number of columns of the given radiation slice
    real(wp), intent(in)  :: latitude_scalar(n_cells_rad)                        ! latitudes of the scalar gridpoints in this radiation block
    real(wp), intent(in)  :: mass_densities(n_cells_rad,n_layers,n_constituents) ! mass densities in this radiation block
    real(wp), intent(in)  :: temperature_gas(n_cells_rad,n_layers)               ! temperature of the gas phase in this radiation block
    real(wp), intent(out) :: radiation_tendency(n_cells_rad,n_layers)            ! resulting heating rate in W/m**3
    
    ! local variables
    integer  :: ji       ! horizontal index
    integer  :: jl       ! layer index
    real(wp) :: pressure ! air pressure
    
    !$omp parallel do private(ji,jl,pressure)
    do jl=1,n_layers
      do ji=1,n_cells_rad
        pressure = mass_densities(ji,jl,n_condensed_constituents+1)*r_d*temperature_gas(ji,jl)
        radiation_tendency(ji,jl) = -k_t(latitude_scalar(ji),pressure) &
                                    *(temperature_gas(ji,jl) - t_eq(latitude_scalar(ji),pressure))
        radiation_tendency(ji,jl) = c_d_v*mass_densities(ji,jl,n_condensed_constituents+1)*radiation_tendency(ji,jl)
      enddo
    enddo
    !$omp end parallel do
    
  end subroutine held_suar
  
  function t_eq(latitude,pressure)
    
    ! This function returns the equilibrium temperature at a certain point in the atmosphere.
    
    real(wp), intent(in) :: latitude ! geographical latitude of the point
    real(wp), intent(in) :: pressure ! pressure at the point (vertical coordinate)
    real(wp)             :: t_eq     ! equilibrium temperature in Kelvin (result)
    
    ! local variables
    real(wp) :: delta_t_y       ! for definition see Held-Suarez paper
    real(wp) :: delta_theta_v_z ! for definition see Held-Suarez paper
    real(wp) :: kappa           ! for definition see Held-Suarez paper
    
    delta_t_y = 60._wp
    delta_theta_v_z = 10._wp
    kappa = 2._wp/7._wp
    t_eq = 315._wp
    t_eq = t_eq - delta_t_y*sin(latitude)**2
    t_eq = t_eq - delta_theta_v_z*log(pressure/p_0)*cos(latitude)**2
    t_eq = t_eq*(pressure/p_0)**kappa
    t_eq = max(200._wp,t_eq)
    
  end function t_eq

  function k_t(latitude,pressure)
    
    ! This function returns the relaxation coefficient according to which the equilibrium temperature is approached.
    
    real(wp), intent(in) :: latitude ! geographical latitude of the point
    real(wp), intent(in) :: pressure ! pressure at the point (vertical coordinate)
    real(wp)             :: k_t      ! relaxation coefficient in 1/s (result)
    
    ! local variables
    real(wp) :: k_a     ! for definition see Held-Suarez paper
    real(wp) :: k_s     ! for definition see Held-Suarez paper
    real(wp) :: sigma_b ! for definition see Held-Suarez paper
    real(wp) :: sigma   ! for definition see Held-Suarez paper
    
    k_a = 1._wp/40._wp*1._wp/86400._wp
    k_s = 1._wp/4._wp*1._wp/86400._wp
    sigma_b = 0.7_wp
    sigma = pressure/p_0
    k_t = k_a + (k_s-k_a)*max(0._wp,(sigma - sigma_b)/(1._wp - sigma_b))*cos(latitude)**4
    
  end function k_t
  
end module mo_held_suarez













