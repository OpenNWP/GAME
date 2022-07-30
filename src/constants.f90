! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module constants

  use iso_c_binding
  use definitions, only: wp
  
  ! This is a collection of some quantities that are hardly ever changed.
  
  implicit none
  
  public
  
  ! physical constants
  ! ------------------
  real(c_double) :: r_e = 6371000.789927_wp ! Earth radius
  real(c_double) :: k_b = 1.380649e-23_wp   ! Boltzmann's constant
  real(c_double) :: n_a =  6.02214076e23_wp ! Avogadro's number
  real(c_double) :: t_0 = 273.15_wp         ! 273.15 K
  real(c_double) :: rho_h2o = 1024._wp      ! typical density of water
  real(c_double) :: p_0 = 100000._wp        ! reference pressure
  real(c_double) :: omega = 7.292115e-5     ! angular frequency of Earth rotation
  real(c_double) :: gravity = 9.80616_wp    ! average surface gravity value
  real(c_double) :: r_d = 287.057811_wp     ! specific gas constant of dry air
  real(c_double) :: c_d_p = 1005._wp        ! isobaric specific heat capacity of dry air
  
  ! non-physical constants
  ! ----------------------
  real(c_double) :: M_PI = 4._wp*atan(1._wp) ! pi
  real(c_double) :: EPSILON_SECURITY = 1e-10 ! security constant
  
  ! some properties of the standard atmosphere
  ! ------------------------------------------
  real(c_double) :: lapse_rate = 0.0065_wp    ! lapse_rate within the troposphere
  real(c_double) :: surface_temp = 288.15_wp  ! the temperature at the surface
  real(c_double) :: tropo_height = 11000._wp  ! the tropopause height
  real(c_double) :: inv_height = 20000._wp    ! height where the temperature inversion begins
  real(c_double) :: t_grad_inv = 0.001_wp     ! temperature gradient above the inversion
  real(c_double) :: p_0_standard = 101325._wp ! reference pressure of the standard atmosphere

end module constants








