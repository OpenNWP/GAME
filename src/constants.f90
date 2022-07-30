! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module constants

  use iso_c_binding
  
  ! This is a collection of some quantities that are hardly ever changed.
  
  implicit none
  
  public
  
  real(c_double) :: re = 6371000.789927      ! Earth radius
  real(c_double) :: k_B = 1.380649e-23       ! Boltzmann's constant
  real(c_double) :: N_A =  6.02214076e23     ! Avogadro's number
  real(c_double) :: T_0 = 273.15             ! 273.15 K
  real(c_double) :: density_water = 1024.    ! typical density of water
  real(c_double) :: p_0 = 100000.            ! reference pressure
  real(c_double) :: omega = 7.292115e-5      ! angular frequency of Earth rotation
  real(c_double) :: gravity = 9.80616_c_double        ! average surface gravity value
  
  ! non-physical constants
  ! ----------------------
  real(c_double) :: M_PI = 4._c_double*atan(1._c_double) ! pi
  real(c_double) :: EPSILON_SECURITY = 1e-10             ! security constant
  
  ! some properties of the standard atmosphere
  ! ------------------------------------------
  real(c_double) :: lapse_rate = 0.0065_c_double      ! lapse_rate within the troposphere
  real(c_double) :: surface_temp = 288.15_c_double    ! the temperature at the surface
  real(c_double) :: tropo_height = 11000._c_double    ! the tropopause height
  real(c_double) :: inv_height = 20000._c_double      ! height where the temperature inversion begins
  real(c_double) :: t_grad_inv = 0.001_c_double       ! temperature gradient above the inversion
  real(c_double) :: p_0_standard = 101325._c_double   ! reference pressure of the standard atmosphere

end module constants








