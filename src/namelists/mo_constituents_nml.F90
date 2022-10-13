! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module mo_constituents_nml

  ! This namelist defines the constituents of the model atmosphere.

  use mo_definitions, only: wp

  implicit none
  
#ifndef COMPILE_TIME_CONFIG
  logical :: lmoist                   ! moisture switch
  integer :: n_gaseous_constituents   ! number of constituents of the gas phase
  integer :: n_condensed_constituents ! number of condensed constituents
  integer :: n_constituents           ! the total number of constituents
#endif
#ifdef COMPILE_TIME_CONFIG
  logical :: lmoist = .true.
  integer :: n_gaseous_constituents = 2
  integer :: n_condensed_constituents = 4
  integer :: n_constituents = 6
  logical :: dummy                        ! dummy variable with no effect
#endif
  real(wp) :: snow_velocity            ! sedimentation velocity of snow
  real(wp) :: rain_velocity            ! sedimentation velocity of rain
  real(wp) :: cloud_droplets_velocity  ! sedimentation velocity of cloud droplets

#ifndef COMPILE_TIME_CONFIG
  namelist /constituents/lmoist
#endif
#ifdef COMPILE_TIME_CONFIG
  namelist /constituents/dummy
#endif

  contains

  subroutine constituents_nml_setup()
  
    ! local variables
    integer :: fileunit
    
#ifndef COMPILE_TIME_CONFIG
    lmoist = .true.
    n_condensed_constituents = 4
    n_gaseous_constituents = 2
#endif
    snow_velocity = 5._wp
    rain_velocity = 10._wp
    cloud_droplets_velocity = .01_wp
    
    ! open and read namelist file
    open(action="read",file="namelist.nml",newunit=fileunit)
    read(nml=constituents,unit=fileunit)
        
    close(fileunit)
    
#ifndef COMPILE_TIME_CONFIG
    ! the dry case
    if (.not. lmoist) then
      n_condensed_constituents = 0
      n_gaseous_constituents = 1
    endif
    n_constituents = n_condensed_constituents + n_gaseous_constituents
#endif
    
    if (lmoist) then
      write(*,*) "Moisture is turned on."
    else
      write(*,*) "Moisture is turned off."
    endif
    
  end subroutine constituents_nml_setup
  
end module mo_constituents_nml






