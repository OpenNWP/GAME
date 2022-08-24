! This source file is part of the Geophysical Fluids Modeling Framework (GAME),which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_held_suarez

  ! This module calculates the Held-Suarez radiative forcing.
  
  use iso_c_binding
  use mo_definitions,      only: wp
  use mo_constants,        only: c_d_v,r_d,p_0
  use mo_constituents_nml, only: n_constituents,n_condensed_constituents
  use mo_rad_nml,          only: n_scals_rad,n_scals_rad_h
  
  implicit none
  
  contains

  subroutine held_suar(latitude_scalar,mass_densities,temperature_gas,radiation_tendency) &
  bind(c,name = "held_suar")
  
    real(wp), intent(in)  :: latitude_scalar(n_scals_rad_h), &
                             mass_densities(n_constituents*n_scals_rad),temperature_gas(n_scals_rad)
    real(wp), intent(out) :: radiation_tendency(n_scals_rad)
    
    ! local variables
    integer  :: ji,layer_index,h_index
    real(wp) :: pressure
  
    !$omp parallel do private(ji,layer_index,h_index,pressure)
    do ji=1,n_scals_rad
      layer_index = (ji-1)/n_scals_rad_h
      h_index = ji - layer_index*n_scals_rad_h
      pressure = mass_densities(n_condensed_constituents*n_scals_rad + ji)*r_d*temperature_gas(ji)
      radiation_tendency(ji) = -k_t(latitude_scalar(h_index),pressure) &
                               *(temperature_gas(ji) - t_eq(latitude_scalar(h_index),pressure))
      radiation_tendency(ji) = c_d_v*mass_densities(n_condensed_constituents*n_scals_rad_h + ji)*radiation_tendency(ji)
    enddo
    !$omp end parallel do
    
  end subroutine held_suar

  function t_eq(latitude,pressure)
  
    real(wp), intent(in) :: latitude,pressure
    real(wp)             :: t_eq
  
    ! local variables
    real(wp) :: delta_t_y,delta_theta_v_z,kappa
  
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
  
    real(wp), intent(in) :: latitude,pressure
    real(wp)             :: k_t
  
    ! local variables
    real(wp) :: k_a,k_s,sigma_b,sigma
  
    k_a = 1._wp/40._wp*1._wp/86400._wp
    k_s = 1._wp/4._wp*1._wp/86400._wp
    sigma_b = 0.7_wp
    sigma = pressure/p_0
    k_t = k_a + (k_s-k_a)*max(0._wp,(sigma - sigma_b)/(1._wp - sigma_b))*cos(latitude)**4
  
  end function k_t

end module mo_held_suarez








