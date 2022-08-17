! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module constituents_nml

  ! This namelist defines the constituents of the model atmosphere.

  use definitions, only: wp
  use run_nml,     only: lmoist

  implicit none
  
  integer  :: n_gaseous_constituents   ! number of constituents of the gas phase
  integer  :: n_condensed_constituents ! number of condensed constituents
  integer  :: n_constituents           ! the total number of constituents
  real(wp) :: snow_velocity            ! sedimentation velocity of snow
  real(wp) :: rain_velocity            ! sedimentation velocity of rain
  real(wp) :: cloud_droplets_velocity  ! sedimentation velocity of cloud droplets

  namelist /constituents/lmoist

  contains

  subroutine constituents_nml_setup() &
  bind(c,name = "constituents_nml_setup")
    
    n_condensed_constituents = 4
    n_gaseous_constituents = 2
    ! the dry case
    if (.not. lmoist) then
      n_condensed_constituents = 0
      n_gaseous_constituents = 1
    endif
    snow_velocity = 5._wp
    rain_velocity = 10._wp
    cloud_droplets_velocity = .01_wp
    n_constituents = n_condensed_constituents + n_gaseous_constituents
    
    if (lmoist) then
      write(*,*) "Moisture is turned on."
    else
      write(*,*) "Moisture is turned off."
    endif
    
  end subroutine constituents_nml_setup
  
end module constituents_nml






