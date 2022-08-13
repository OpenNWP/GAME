! This source file is part of the Geophysical Fluids Modeling Framework (GAME), which is released under the MIT license.
! Github repository: https://github.com/OpenNWP/GAME

module run_nml

  ! This is the namelist that configures the basic run properties of a model integration.
  
  use iso_c_binding
  use definitions, only: wp
  
  implicit none
  
  character(len=64) :: run_id ! ID of this run
  real(wp)          :: dtime  ! time step
  logical           :: lmoist ! moisture switch
  
  namelist /run/run_id,lmoist

  contains

  subroutine run_nml_setup() &
  bind(c,name = "run_nml_setup")
  
    run_id = "ideal"
    dtime = 360.312923_wp
    lmoist = .true.
  
  end subroutine run_nml_setup
  
end module run_nml












