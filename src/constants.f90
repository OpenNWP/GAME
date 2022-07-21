! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module constants

  ! This is a collection of some quantities that are hardly ever changed.
  
  implicit none
  
  public
  
  real :: re = 6371000.789927      ! Earth radius
  real :: k_B = 1.380649e-23       ! Boltzmann's constant
  real :: N_A =  6.02214076e23     ! Avogadro's number
  real :: T_0 = 273.15             ! 273.15 K
  real :: density_water = 1024.    ! typical density of water
  real :: p_0 = 100000.            ! reference pressure
  real :: omega = 7.292115e-5      ! angular frequency of Earth rotation
  real :: gravity = 9.80616        ! average surface gravity value
  
  ! non-physical constants
  ! ----------------------
  real :: M_PI = 4.*atan(1.)       ! pi
  real :: EPSILON_SECURITY = 1e-10 ! security constant
  
  ! some properties of the standard atmosphere
  ! ------------------------------------------
  real :: lapse_rate = 0.0065      ! lapse_rate within the troposphere
  real :: surface_temp = 288.15    ! the temperature at the surface
  real :: tropo_height = 12000.    ! the tropopause height
  real :: inv_height = 20000.      ! height where the temperature inversion begins
  real :: t_grad_inv = 0.001       ! temperature gradient above the inversion
  real :: p_0_standard = 101325.   ! reference pressure of the standard atmosphere

end module constants








