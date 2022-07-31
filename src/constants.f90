! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module constants

  ! This is a collection of some quantities that are hardly ever changed.
  
  use definitions, only: wp
  
  implicit none
  
  public
  
  ! physical constants
  ! ------------------
  real(wp), parameter :: r_e = 6371000.789927_wp   ! Earth radius
  real(wp), parameter :: k_b = 1.380649e-23_wp     ! Boltzmann's constant
  real(wp), parameter :: n_a =  6.02214076e23_wp   ! Avogadro's number
  real(wp), parameter :: t_0 = 273.15_wp           ! 273.15 K
  real(wp), parameter :: rho_h2o = 1024._wp        ! typical density of water
  real(wp), parameter :: p_0 = 100000._wp          ! reference pressure
  real(wp), parameter :: omega = 7.292115e-5_wp    ! angular frequency of Earth rotation
  real(wp), parameter :: gravity = 9.80616_wp      ! average surface gravity value
  real(wp), parameter :: m_d = n_a*0.004810e-23_wp ! molar mass of dry air
  real(wp), parameter :: m_v = n_a*0.002991e-23_wp ! molar mass of water
  real(wp), parameter :: r_d = 287.057811_wp       ! specific gas constant of dry air
  real(wp), parameter :: r_v = 461.524879_wp       ! specific gas constant of water vapour
  real(wp), parameter :: c_d_v = 717.942189_wp     ! isochoric specific heat capacity of dry air
  real(wp), parameter :: c_v_v = 1396.475121_wp    ! isochoric specific heat capacity of water vapour
  real(wp), parameter :: c_d_p = 1005._wp          ! isobaric specific heat capacity of dry air
  real(wp), parameter :: c_v_p = 1858._wp          ! isobaric specific heat capacity of water vapour
  
  ! non-physical constants
  ! ----------------------
  real(wp), parameter :: M_PI = 4._wp*atan(1._wp)    ! pi
  real(wp), parameter :: EPSILON_SECURITY = 1e-10_wp ! security constant
  
  ! some properties of the standard atmosphere
  ! ------------------------------------------
  real(wp), parameter :: lapse_rate = 0.0065_wp    ! lapse_rate within the troposphere
  real(wp), parameter :: surface_temp = 288.15_wp  ! the temperature at the surface
  real(wp), parameter :: tropo_height = 11000._wp  ! the tropopause height
  real(wp), parameter :: inv_height = 20000._wp    ! height where the temperature inversion begins
  real(wp), parameter :: t_grad_inv = 0.001_wp     ! temperature gradient above the inversion
  real(wp), parameter :: p_0_standard = 101325._wp ! reference pressure of the standard atmosphere

end module constants







